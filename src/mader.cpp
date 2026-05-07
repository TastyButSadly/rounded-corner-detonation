#include "mader.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "grid.h"
#include "solver.h"

namespace {

constexpr double EPS = 1e-8;

struct FaceState {
    double rho = 0.0;
    double u = 0.0;
    double v = 0.0;
    double E = 0.0;
    double w = 0.0;
    bool blocked = true;
    bool physical = false;
    int i = -1;
    int j = -1;
};

bool is_physical_cell(int i, int j, int Nx, int Ny, int fict_x, int fict_y) {
    return i >= fict_x && i < Nx - fict_x && j >= fict_y && j < Ny - fict_y;
}

bool is_fluid_cell(int i, int j, int Nx, int Ny, int fict_x, int fict_y, const Mask2D& solid_mask) {
    return is_physical_cell(i, j, Nx, Ny, fict_x, fict_y) && solid_mask[i][j] == 0;
}

double clamp_unit(double value) {
    return std::max(0.0, std::min(1.0, value));
}

int physical_i_index(int i, int fict_x) {
    return i - fict_x + 1;
}

double radial_geom_factor(int donor_i, int fict_x, double alpha_abs, bool slab) {
    if (slab) {
        return 1.0;
    }
    const double ii = static_cast<double>(physical_i_index(donor_i, fict_x));
    const double denom = std::max(2.0 * ii - 1.0, EPS);
    return std::max(0.0, (2.0 * ii - alpha_abs) / denom);
}

double boundary_velocity(double value, int bc_code) {
    if (bc_code == 1) {
        return -value;
    }
    return value;
}

std::vector<double> make_cons(double rho, double u, double v, double p, double g) {
    return noncons_to_cons(std::vector<double>{rho, u, v, p}, g);
}

std::vector<double> boundary_cons(double rho, double u, double v, double p,
                                  int bc_code, int axis, double g) {
    return boundary(make_cons(rho, u, v, p, g), bc_code, g, axis);
}

double cell_internal_energy(double rho, double pressure, double g) {
    if (rho <= EPS) {
        return EPS;
    }
    return std::max(pressure / (rho * (g - 1.0)), EPS);
}

double corrected_transport_fraction(const ScalarField& w_field,
                                    int donor_i, int donor_j,
                                    int acceptor_i, int acceptor_j,
                                    int Nx, int Ny, int fict_x, int fict_y,
                                    const Mask2D& solid_mask,
                                    const MaderConfig& config) {
    const double donor_w = clamp_unit(w_field[donor_i][donor_j]);
    if (!config.shargatov_correction) {
        return donor_w;
    }
    if (donor_w <= config.gasw_threshold || donor_w >= 1.0 - config.gasw_threshold) {
        return donor_w;
    }

    auto classify = [&](int i, int j) -> int {
        if (!is_fluid_cell(i, j, Nx, Ny, fict_x, fict_y, solid_mask)) {
            return -1;
        }
        const double value = clamp_unit(w_field[i][j]);
        if (value <= config.gasw_threshold) {
            return 0;
        }
        if (value >= 1.0 - config.gasw_threshold) {
            return 1;
        }
        return -1;
    };

    std::vector<std::pair<int, int>> common_neighbors;
    if (donor_i != acceptor_i) {
        common_neighbors.push_back({std::min(donor_i, acceptor_i), donor_j - 1});
        common_neighbors.push_back({std::min(donor_i, acceptor_i), donor_j + 1});
    } else if (donor_j != acceptor_j) {
        common_neighbors.push_back({donor_i - 1, std::min(donor_j, acceptor_j)});
        common_neighbors.push_back({donor_i + 1, std::min(donor_j, acceptor_j)});
    }

    bool has_pure_zero = false;
    bool has_pure_one = false;
    constexpr int offsets[4][2] = {
        {-1, 0}, {1, 0}, {0, -1}, {0, 1}
    };

    for (const auto& offset : offsets) {
        const int ni = donor_i + offset[0];
        const int nj = donor_j + offset[1];
        const int cls = classify(ni, nj);
        has_pure_zero = has_pure_zero || cls == 0;
        has_pure_one = has_pure_one || cls == 1;
    }

    for (const auto& neighbor : common_neighbors) {
        const int cls = classify(neighbor.first, neighbor.second);
        has_pure_zero = has_pure_zero || cls == 0;
        has_pure_one = has_pure_one || cls == 1;
    }

    if (acceptor_i >= 0 && acceptor_j >= 0) {
        const int acceptor_cls = classify(acceptor_i, acceptor_j);
        has_pure_zero = has_pure_zero || acceptor_cls == 0;
        has_pure_one = has_pure_one || acceptor_cls == 1;
    }

    if (has_pure_zero && !has_pure_one) {
        return 1.0;
    }
    if (has_pure_one && !has_pure_zero) {
        return 0.0;
    }
    return donor_w;
}

double right_boundary_alpha(const ScalarField& u_tilde,
                            int donor_i,
                            int j,
                            int fict_x,
                            double dt,
                            double dx) {
    const int inner_i = std::max(donor_i - 1, fict_x);
    return 0.5 * (3.0 * u_tilde[donor_i][j] - u_tilde[inner_i][j]) * dt / dx;
}

double top_boundary_beta(const ScalarField& v_tilde,
                         int i,
                         int donor_j,
                         double dt,
                         double dy) {
    return v_tilde[i][donor_j] * dt / dy;
}

double transport_fraction_x(const FaceState& left,
                            const FaceState& right,
                            const ScalarField& u_tilde,
                            int j,
                            int fict_x,
                            int right_bc_code,
                            double dt,
                            double dx) {
    if (right_bc_code == 2 && left.physical && !right.physical && right.i < 0) {
        return right_boundary_alpha(u_tilde, left.i, j, fict_x, dt, dx);
    }
    const double left_velocity = left.physical ? u_tilde[left.i][j] : left.u;
    const double right_velocity = right.physical ? u_tilde[right.i][j] : right.u;
    double denom = 1.0 + (left_velocity - right_velocity) * dt / dx;
    if (std::abs(denom) < EPS) {
        denom = (denom >= 0.0) ? EPS : -EPS;
    }
    return 0.5 * (left_velocity + right_velocity) * dt / dx / denom;
}

double transport_fraction_y(const FaceState& down,
                            const FaceState& up,
                            const ScalarField& v_tilde,
                            int i,
                            int up_bc_code,
                            double dt,
                            double dy) {
    if (up_bc_code == 2 && down.physical && !up.physical && up.j < 0) {
        return top_boundary_beta(v_tilde, i, down.j, dt, dy);
    }
    const double down_velocity = down.physical ? v_tilde[i][down.j] : down.v;
    const double up_velocity = up.physical ? v_tilde[i][up.j] : up.v;
    double denom = 1.0 + (down_velocity - up_velocity) * dt / dy;
    if (std::abs(denom) < EPS) {
        denom = (denom >= 0.0) ? EPS : -EPS;
    }
    return 0.5 * (down_velocity + up_velocity) * dt / dy / denom;
}

double sample_pressure(double rho, double u, double v, double p,
                       double neighbor_rho, double neighbor_u, double neighbor_v, double neighbor_p,
                       int i, int j, int ni, int nj,
                       int Nx, int Ny, int fict_x, int fict_y,
                       int left_bc_code, int right_bc_code, int up_bc_code, int down_bc_code,
                       const Mask2D& solid_mask,
                       double g) {
    if (ni < fict_x) {
        return cons_to_noncons(boundary_cons(rho, u, v, p, left_bc_code, 0, g), g)[P];
    }
    if (ni >= Nx - fict_x) {
        return cons_to_noncons(boundary_cons(rho, u, v, p, right_bc_code, 0, g), g)[P];
    }
    if (nj < fict_y) {
        return cons_to_noncons(boundary_cons(rho, u, v, p, down_bc_code, 1, g), g)[P];
    }
    if (nj >= Ny - fict_y) {
        return cons_to_noncons(boundary_cons(rho, u, v, p, up_bc_code, 1, g), g)[P];
    }
    if (solid_mask[ni][nj]) {
        const int axis = (ni != i) ? 0 : 1;
        return cons_to_noncons(boundary_cons(rho, u, v, p, 1, axis, g), g)[P];
    }
    (void)neighbor_rho;
    (void)neighbor_u;
    (void)neighbor_v;
    return neighbor_p;
}

double sample_normal_velocity(const ScalarField& velocity,
                              int i, int j, int ni, int nj,
                              int Nx, int Ny, int fict_x, int fict_y,
                              int left_bc_code, int right_bc_code, int up_bc_code, int down_bc_code,
                              const Mask2D& solid_mask,
                              int axis) {
    const double current = velocity[i][j];
    if (ni < fict_x) {
        return boundary_velocity(current, left_bc_code);
    }
    if (ni >= Nx - fict_x) {
        return boundary_velocity(current, right_bc_code);
    }
    if (nj < fict_y) {
        return boundary_velocity(current, down_bc_code);
    }
    if (nj >= Ny - fict_y) {
        return boundary_velocity(current, up_bc_code);
    }
    if (solid_mask[ni][nj]) {
        return -current;
    }
    (void)axis;
    return velocity[ni][nj];
}

FaceState make_vertical_side(int cell_i, int j,
                             int Nx, int Ny, int fict_x, int fict_y,
                             int left_bc_code, int right_bc_code,
                             const Mask2D& solid_mask,
                             const ScalarField& rho,
                             const ScalarField& u,
                             const ScalarField& v,
                             const ScalarField& E,
                             const ScalarField& w) {
    FaceState state;
    const int first_phys = fict_x;
    const int last_phys = Nx - fict_x - 1;

    if (cell_i < first_phys) {
        if (left_bc_code == 1) {
            return state;
        }
        const int ref_i = first_phys;
        if (!is_fluid_cell(ref_i, j, Nx, Ny, fict_x, fict_y, solid_mask)) {
            return state;
        }
        state.rho = rho[ref_i][j];
        state.u = u[ref_i][j];
        state.v = v[ref_i][j];
        state.E = E[ref_i][j];
        state.w = clamp_unit(w[ref_i][j]);
        state.blocked = false;
        state.physical = false;
        return state;
    }

    if (cell_i > last_phys) {
        if (right_bc_code == 1) {
            return state;
        }
        const int ref_i = last_phys;
        if (!is_fluid_cell(ref_i, j, Nx, Ny, fict_x, fict_y, solid_mask)) {
            return state;
        }
        state.rho = rho[ref_i][j];
        state.u = u[ref_i][j];
        state.v = v[ref_i][j];
        state.E = E[ref_i][j];
        state.w = clamp_unit(w[ref_i][j]);
        state.blocked = false;
        state.physical = false;
        return state;
    }

    if (!is_fluid_cell(cell_i, j, Nx, Ny, fict_x, fict_y, solid_mask)) {
        return state;
    }

    state.rho = rho[cell_i][j];
    state.u = u[cell_i][j];
    state.v = v[cell_i][j];
    state.E = E[cell_i][j];
    state.w = clamp_unit(w[cell_i][j]);
    state.blocked = false;
    state.physical = true;
    state.i = cell_i;
    state.j = j;
    return state;
}

FaceState make_horizontal_side(int i, int cell_j,
                               int Nx, int Ny, int fict_x, int fict_y,
                               int up_bc_code, int down_bc_code,
                               const Mask2D& solid_mask,
                               const ScalarField& rho,
                               const ScalarField& u,
                               const ScalarField& v,
                               const ScalarField& E,
                               const ScalarField& w) {
    FaceState state;
    const int first_phys = fict_y;
    const int last_phys = Ny - fict_y - 1;

    if (cell_j < first_phys) {
        if (down_bc_code == 1) {
            return state;
        }
        const int ref_j = first_phys;
        if (!is_fluid_cell(i, ref_j, Nx, Ny, fict_x, fict_y, solid_mask)) {
            return state;
        }
        state.rho = rho[i][ref_j];
        state.u = u[i][ref_j];
        state.v = v[i][ref_j];
        state.E = E[i][ref_j];
        state.w = clamp_unit(w[i][ref_j]);
        state.blocked = false;
        state.physical = false;
        return state;
    }

    if (cell_j > last_phys) {
        if (up_bc_code == 1) {
            return state;
        }
        const int ref_j = last_phys;
        if (!is_fluid_cell(i, ref_j, Nx, Ny, fict_x, fict_y, solid_mask)) {
            return state;
        }
        state.rho = rho[i][ref_j];
        state.u = u[i][ref_j];
        state.v = v[i][ref_j];
        state.E = E[i][ref_j];
        state.w = clamp_unit(w[i][ref_j]);
        state.blocked = false;
        state.physical = false;
        return state;
    }

    if (!is_fluid_cell(i, cell_j, Nx, Ny, fict_x, fict_y, solid_mask)) {
        return state;
    }

    state.rho = rho[i][cell_j];
    state.u = u[i][cell_j];
    state.v = v[i][cell_j];
    state.E = E[i][cell_j];
    state.w = clamp_unit(w[i][cell_j]);
    state.blocked = false;
    state.physical = true;
    state.i = i;
    state.j = cell_j;
    return state;
}

void accumulate_transport(ScalarField& DM,
                          ScalarField& DE,
                          ScalarField& DW,
                          ScalarField& DPU,
                          ScalarField& DPV,
                          const FaceState& donor,
                          const FaceState& acceptor,
                          double amount,
                          double transported_w) {
    if (amount <= 0.0) {
        return;
    }

    if (donor.physical) {
        DM[donor.i][donor.j] -= amount;
        DE[donor.i][donor.j] -= donor.E * amount;
        DW[donor.i][donor.j] -= transported_w * amount;
        DPU[donor.i][donor.j] -= donor.u * amount;
        DPV[donor.i][donor.j] -= donor.v * amount;
    }
    if (acceptor.physical) {
        DM[acceptor.i][acceptor.j] += amount;
        DE[acceptor.i][acceptor.j] += donor.E * amount;
        DW[acceptor.i][acceptor.j] += transported_w * amount;
        DPU[acceptor.i][acceptor.j] += donor.u * amount;
        DPV[acceptor.i][acceptor.j] += donor.v * amount;
    }
}

} // namespace

ScalarField make_scalar_field(int Nx, int Ny, double value) {
    return ScalarField(Nx, std::vector<double>(Ny, value));
}

Mask2D build_step_mask(double x_min, double y_min, double dx, double dy,
                       int Nx, int Ny, int fict_x, int fict_y,
                       const MaderConfig& config) {
    Mask2D solid_mask(Nx, std::vector<unsigned char>(Ny, 0));
    if (!config.step_enabled) {
        return solid_mask;
    }

    const double radius = config.corner_radius;
    const bool rounded_corner = radius > EPS;
    const double corner_center_x = config.step_x_end - radius;
    const double corner_center_y = config.step_y_end - radius;
    const double radius_sq = radius * radius;

    for (int i = fict_x; i < Nx - fict_x; ++i) {
        const double x = x_min + (i - fict_x + 0.5) * dx;
        for (int j = fict_y; j < Ny - fict_y; ++j) {
            const double y = y_min + (j - fict_y + 0.5) * dy;
            if (x <= config.step_x_end + EPS && y <= config.step_y_end + EPS) {
                bool is_solid = true;
                if (rounded_corner && x > corner_center_x + EPS && y > corner_center_y + EPS) {
                    const double dx_corner = x - corner_center_x;
                    const double dy_corner = y - corner_center_y;
                    is_solid = dx_corner * dx_corner + dy_corner * dy_corner <= radius_sq + EPS;
                }
                solid_mask[i][j] = is_solid ? 1 : 0;
            }
        }
    }

    return solid_mask;
}

void apply_step_mask(Field3& u,
                     const Mask2D& solid_mask,
                     const std::vector<double>& fill_cons) {
    for (size_t i = 0; i < u.size(); ++i) {
        for (size_t j = 0; j < u[i].size(); ++j) {
            if (solid_mask[i][j]) {
                u[i][j] = fill_cons;
            }
        }
    }
}

void apply_step_mask_scalar(ScalarField& field,
                            const Mask2D& solid_mask,
                            double fill_value) {
    for (size_t i = 0; i < field.size(); ++i) {
        for (size_t j = 0; j < field[i].size(); ++j) {
            if (solid_mask[i][j]) {
                field[i][j] = fill_value;
            }
        }
    }
}

double mader_temperature(double rho, double pressure, const MaderConfig& config) {
    const double gas_constant = std::max(config.gas_constant, EPS);
    return pressure / std::max(rho * gas_constant, EPS);
}

void exchange_halos_scalar(ScalarField& field,
                           MPI_Comm cart_comm,
                           int halo_x,
                           int halo_y,
                           int i_start,
                           int i_end,
                           int j_start,
                           int j_end,
                           int left_rank,
                           int right_rank,
                           int down_rank,
                           int up_rank) {
    MPI_Status status;

    const int y_len = j_end - j_start;
    const int x_len = (i_end - i_start)
        + (left_rank != MPI_PROC_NULL ? halo_x : 0)
        + (right_rank != MPI_PROC_NULL ? halo_x : 0);

    std::vector<double> send_left(static_cast<size_t>(halo_x * y_len), 0.0);
    std::vector<double> recv_left(static_cast<size_t>(halo_x * y_len), 0.0);
    std::vector<double> send_right(static_cast<size_t>(halo_x * y_len), 0.0);
    std::vector<double> recv_right(static_cast<size_t>(halo_x * y_len), 0.0);

    if (left_rank != MPI_PROC_NULL) {
        int idx = 0;
        for (int i = 0; i < halo_x; ++i) {
            for (int j = j_start; j < j_end; ++j) {
                send_left[idx++] = field[i_start + i][j];
            }
        }
    }

    if (right_rank != MPI_PROC_NULL) {
        int idx = 0;
        for (int i = 0; i < halo_x; ++i) {
            for (int j = j_start; j < j_end; ++j) {
                send_right[idx++] = field[i_end - halo_x + i][j];
            }
        }
    }

    MPI_Sendrecv(send_left.data(), static_cast<int>(send_left.size()), MPI_DOUBLE, left_rank, 10,
                 recv_right.data(), static_cast<int>(recv_right.size()), MPI_DOUBLE, right_rank, 10,
                 cart_comm, &status);
    MPI_Sendrecv(send_right.data(), static_cast<int>(send_right.size()), MPI_DOUBLE, right_rank, 11,
                 recv_left.data(), static_cast<int>(recv_left.size()), MPI_DOUBLE, left_rank, 11,
                 cart_comm, &status);

    if (left_rank != MPI_PROC_NULL) {
        int idx = 0;
        for (int i = 0; i < halo_x; ++i) {
            for (int j = j_start; j < j_end; ++j) {
                field[i_start - halo_x + i][j] = recv_left[idx++];
            }
        }
    }

    if (right_rank != MPI_PROC_NULL) {
        int idx = 0;
        for (int i = 0; i < halo_x; ++i) {
            for (int j = j_start; j < j_end; ++j) {
                field[i_end + i][j] = recv_right[idx++];
            }
        }
    }

    const int x_begin = (left_rank != MPI_PROC_NULL) ? i_start - halo_x : i_start;
    const int x_end = (right_rank != MPI_PROC_NULL) ? i_end + halo_x : i_end;
    std::vector<double> send_down(static_cast<size_t>(halo_y * x_len), 0.0);
    std::vector<double> recv_down(static_cast<size_t>(halo_y * x_len), 0.0);
    std::vector<double> send_up(static_cast<size_t>(halo_y * x_len), 0.0);
    std::vector<double> recv_up(static_cast<size_t>(halo_y * x_len), 0.0);

    if (down_rank != MPI_PROC_NULL) {
        int idx = 0;
        for (int i = x_begin; i < x_end; ++i) {
            for (int j = 0; j < halo_y; ++j) {
                send_down[idx++] = field[i][j_start + j];
            }
        }
    }

    if (up_rank != MPI_PROC_NULL) {
        int idx = 0;
        for (int i = x_begin; i < x_end; ++i) {
            for (int j = 0; j < halo_y; ++j) {
                send_up[idx++] = field[i][j_end - halo_y + j];
            }
        }
    }

    MPI_Sendrecv(send_down.data(), static_cast<int>(send_down.size()), MPI_DOUBLE, down_rank, 12,
                 recv_up.data(), static_cast<int>(recv_up.size()), MPI_DOUBLE, up_rank, 12,
                 cart_comm, &status);
    MPI_Sendrecv(send_up.data(), static_cast<int>(send_up.size()), MPI_DOUBLE, up_rank, 13,
                 recv_down.data(), static_cast<int>(recv_down.size()), MPI_DOUBLE, down_rank, 13,
                 cart_comm, &status);

    if (down_rank != MPI_PROC_NULL) {
        int idx = 0;
        for (int i = x_begin; i < x_end; ++i) {
            for (int j = 0; j < halo_y; ++j) {
                field[i][j_start - halo_y + j] = recv_down[idx++];
            }
        }
    }

    if (up_rank != MPI_PROC_NULL) {
        int idx = 0;
        for (int i = x_begin; i < x_end; ++i) {
            for (int j = 0; j < halo_y; ++j) {
                field[i][j_end + j] = recv_up[idx++];
            }
        }
    }
}

void gather_scalar_to_root(ScalarField& field,
                           int Nx,
                           int Ny,
                           int fict_x,
                           int fict_y,
                           int i_start,
                           int i_end,
                           int j_start,
                           int j_end,
                           int rank,
                           int size,
                           MPI_Comm cart_comm,
                           int* dims) {
    const int local_nx = i_end - i_start;
    const int local_ny = j_end - j_start;
    std::vector<double> send_buf(static_cast<size_t>(local_nx * local_ny), 0.0);

    int idx = 0;
    for (int i = i_start; i < i_end; ++i) {
        for (int j = j_start; j < j_end; ++j) {
            send_buf[idx++] = field[i][j];
        }
    }

    if (rank == 0) {
        for (int p = 1; p < size; ++p) {
            int coords[2] = {0, 0};
            MPI_Cart_coords(cart_comm, p, 2, coords);

            int p_i_start = 0;
            int p_i_end = 0;
            int p_j_start = 0;
            int p_j_end = 0;
            get_subdomain_bounds(Nx, dims[0], coords[0], p_i_start, p_i_end);
            get_subdomain_bounds(Ny, dims[1], coords[1], p_j_start, p_j_end);
            p_i_start += fict_x;
            p_i_end += fict_x;
            p_j_start += fict_y;
            p_j_end += fict_y;

            std::vector<double> recv_buf(static_cast<size_t>((p_i_end - p_i_start) * (p_j_end - p_j_start)), 0.0);
            MPI_Recv(recv_buf.data(), static_cast<int>(recv_buf.size()), MPI_DOUBLE, p, 20, cart_comm, MPI_STATUS_IGNORE);

            int recv_idx = 0;
            for (int i = p_i_start; i < p_i_end; ++i) {
                for (int j = p_j_start; j < p_j_end; ++j) {
                    field[i][j] = recv_buf[recv_idx++];
                }
            }
        }
    } else {
        MPI_Send(send_buf.data(), static_cast<int>(send_buf.size()), MPI_DOUBLE, 0, 20, cart_comm);
    }
}

double calc_time_step_masked(const Field3& v_cons,
                             const Mask2D& solid_mask,
                             double dx,
                             double dy,
                             double cfl,
                             double g,
                             int fict_x,
                             int fict_y) {
    double max_speed = 0.0;
    const int Nx = static_cast<int>(v_cons.size());
    const int Ny = static_cast<int>(v_cons[0].size());

    for (int i = fict_x; i < Nx - fict_x; ++i) {
        for (int j = fict_y; j < Ny - fict_y; ++j) {
            if (solid_mask[i][j]) {
                continue;
            }
            const std::vector<double> prim = cons_to_noncons(v_cons[i][j], g);
            const double c = calc_sound_speed(prim, g);
            const double flow_speed = std::sqrt(prim[U] * prim[U] + prim[V] * prim[V]);
            max_speed = std::max(max_speed, c + flow_speed);
        }
    }

    if (max_speed <= EPS) {
        return cfl * std::min(dx, dy);
    }
    return cfl * std::min(dx, dy) / max_speed;
}

double calc_time_step_masked_local(const Field3& v_cons,
                                   const Mask2D& solid_mask,
                                   double dx,
                                   double dy,
                                   double cfl,
                                   double g,
                                   int fict_x,
                                   int fict_y,
                                   int i_start,
                                   int i_end,
                                   int j_start,
                                   int j_end) {
    double max_speed = 0.0;
    bool has_fluid_cells = false;

    for (int i = i_start; i < i_end; ++i) {
        for (int j = j_start; j < j_end; ++j) {
            if (!is_fluid_cell(i, j, static_cast<int>(v_cons.size()), static_cast<int>(v_cons[0].size()),
                               fict_x, fict_y, solid_mask)) {
                continue;
            }
            has_fluid_cells = true;
            const std::vector<double> prim = cons_to_noncons(v_cons[i][j], g);
            const double c = calc_sound_speed(prim, g);
            const double flow_speed = std::sqrt(prim[U] * prim[U] + prim[V] * prim[V]);
            max_speed = std::max(max_speed, c + flow_speed);
        }
    }

    if (!has_fluid_cells || max_speed <= EPS) {
        return std::numeric_limits<double>::max();
    }
    return cfl * std::min(dx, dy) / max_speed;
}

void mader_step_2d(const Field3& u_prev,
                   Field3& u_next,
                   const ScalarField& w_prev,
                   ScalarField& w_next,
                   const ScalarField& rho_ref,
                   double dt,
                   double dx,
                   double dy,
                   double g,
                   int Nx,
                   int Ny,
                   int fict_x,
                   int fict_y,
                   int left_bc_code,
                   int right_bc_code,
                   int up_bc_code,
                   int down_bc_code,
                   const Mask2D& solid_mask,
                   const MaderConfig& config,
                   int step_index) {
    mader_step_2d_local(u_prev, u_next, w_prev, w_next, rho_ref,
                        dt, dx, dy, g, Nx, Ny, fict_x, fict_y,
                        left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                        solid_mask, config, step_index,
                        fict_x, Nx - fict_x, fict_y, Ny - fict_y);
}

void mader_step_2d_local(const Field3& u_prev,
                         Field3& u_next,
                         const ScalarField& w_prev,
                         ScalarField& w_next,
                         const ScalarField& rho_ref,
                         double dt,
                         double dx,
                         double dy,
                         double g,
                         int Nx,
                         int Ny,
                         int fict_x,
                         int fict_y,
                         int left_bc_code,
                         int right_bc_code,
                         int up_bc_code,
                         int down_bc_code,
                         const Mask2D& solid_mask,
                         const MaderConfig& config,
                         int step_index,
                         int i_start,
                         int i_end,
                         int j_start,
                         int j_end) {
    u_next = u_prev;
    w_next = w_prev;

    const int phys_i_begin = fict_x;
    const int phys_i_end = Nx - fict_x;
    const int phys_j_begin = fict_y;
    const int phys_j_end = Ny - fict_y;

    const int phase_i_begin = std::max(phys_i_begin, i_start - 1);
    const int phase_i_end = std::min(phys_i_end, i_end + 1);
    const int phase_j_begin = std::max(phys_j_begin, j_start - 1);
    const int phase_j_end = std::min(phys_j_end, j_end + 1);

    const int work_i_begin = std::max(phys_i_begin, phase_i_begin - 1);
    const int work_i_end = std::min(phys_i_end, phase_i_end + 1);
    const int work_j_begin = std::max(phys_j_begin, phase_j_begin - 1);
    const int work_j_end = std::min(phys_j_end, phase_j_end + 1);

    ScalarField rho = make_scalar_field(Nx, Ny, 0.0);
    ScalarField u = make_scalar_field(Nx, Ny, 0.0);
    ScalarField v = make_scalar_field(Nx, Ny, 0.0);
    ScalarField p = make_scalar_field(Nx, Ny, 0.0);
    ScalarField I = make_scalar_field(Nx, Ny, 0.0);
    ScalarField temperature = make_scalar_field(Nx, Ny, 0.0);
    ScalarField w_reacted = w_prev;

    for (int i = work_i_begin; i < work_i_end; ++i) {
        for (int j = work_j_begin; j < work_j_end; ++j) {
            if (!is_fluid_cell(i, j, Nx, Ny, fict_x, fict_y, solid_mask)) {
                continue;
            }

            const std::vector<double> prim = cons_to_noncons(u_prev[i][j], g);
            rho[i][j] = prim[RHO];
            u[i][j] = prim[U];
            v[i][j] = prim[V];
            p[i][j] = prim[P];
            I[i][j] = cell_internal_energy(prim[RHO], prim[P], g);
            temperature[i][j] = mader_temperature(rho[i][j], p[i][j], config);

            if (config.reaction_enabled
                && step_index >= config.reaction_delay_steps
                && temperature[i][j] > config.min_temperature
                && w_reacted[i][j] > config.gasw_threshold) {
                const double exponent = -config.reaction_activation_energy
                    / std::max(config.gas_constant * temperature[i][j], EPS);
                const double reaction_rate = config.reaction_rate * w_reacted[i][j] * std::exp(exponent);
                const double delta_w = std::min(w_reacted[i][j], dt * reaction_rate);
                w_reacted[i][j] = clamp_unit(w_reacted[i][j] - delta_w);
                if (w_reacted[i][j] < config.gasw_threshold) {
                    w_reacted[i][j] = 0.0;
                }

                I[i][j] += config.reaction_heat_release * delta_w;
                p[i][j] = std::max((g - 1.0) * rho[i][j] * I[i][j], EPS);
                temperature[i][j] = mader_temperature(rho[i][j], p[i][j], config);
            }
        }
    }

    ScalarField qx = make_scalar_field(Nx, Ny, 0.0);
    ScalarField qy = make_scalar_field(Nx, Ny, 0.0);

    for (int i = phase_i_begin; i < std::min(phase_i_end, phys_i_end - 1); ++i) {
        for (int j = phase_j_begin; j < phase_j_end; ++j) {
            if (!is_fluid_cell(i, j, Nx, Ny, fict_x, fict_y, solid_mask)
                || !is_fluid_cell(i + 1, j, Nx, Ny, fict_x, fict_y, solid_mask)) {
                continue;
            }
            const double du = u[i][j] - u[i + 1][j];
            if (du > 0.0) {
                qx[i][j] = config.visc * rho[i][j] * du;
            }
        }
    }

    for (int i = phase_i_begin; i < phase_i_end; ++i) {
        for (int j = phase_j_begin; j < std::min(phase_j_end, phys_j_end - 1); ++j) {
            if (!is_fluid_cell(i, j, Nx, Ny, fict_x, fict_y, solid_mask)
                || !is_fluid_cell(i, j + 1, Nx, Ny, fict_x, fict_y, solid_mask)) {
                continue;
            }
            const double dv = v[i][j] - v[i][j + 1];
            if (dv > 0.0) {
                qy[i][j] = config.visc * rho[i][j] * dv;
            }
        }
    }

    ScalarField u_tilde = u;
    ScalarField v_tilde = v;

    for (int i = phase_i_begin; i < phase_i_end; ++i) {
        for (int j = phase_j_begin; j < phase_j_end; ++j) {
            if (!is_fluid_cell(i, j, Nx, Ny, fict_x, fict_y, solid_mask)) {
                continue;
            }

            const double p_left = sample_pressure(rho[i][j], u[i][j], v[i][j], p[i][j],
                                                  rho[i - 1][j], u[i - 1][j], v[i - 1][j], p[i - 1][j],
                                                  i, j, i - 1, j,
                                                  Nx, Ny, fict_x, fict_y,
                                                  left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                                                  solid_mask, g);
            const double p_right = sample_pressure(rho[i][j], u[i][j], v[i][j], p[i][j],
                                                   rho[i + 1][j], u[i + 1][j], v[i + 1][j], p[i + 1][j],
                                                   i, j, i + 1, j,
                                                   Nx, Ny, fict_x, fict_y,
                                                   left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                                                   solid_mask, g);
            const double p_down = sample_pressure(rho[i][j], u[i][j], v[i][j], p[i][j],
                                                  rho[i][j - 1], u[i][j - 1], v[i][j - 1], p[i][j - 1],
                                                  i, j, i, j - 1,
                                                  Nx, Ny, fict_x, fict_y,
                                                  left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                                                  solid_mask, g);
            const double p_up = sample_pressure(rho[i][j], u[i][j], v[i][j], p[i][j],
                                                rho[i][j + 1], u[i][j + 1], v[i][j + 1], p[i][j + 1],
                                                i, j, i, j + 1,
                                                Nx, Ny, fict_x, fict_y,
                                                left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                                                solid_mask, g);

            const double q_left = (i > fict_x && is_fluid_cell(i - 1, j, Nx, Ny, fict_x, fict_y, solid_mask)) ? qx[i - 1][j] : 0.0;
            const double q_right = (i < Nx - fict_x - 1 && is_fluid_cell(i + 1, j, Nx, Ny, fict_x, fict_y, solid_mask)) ? qx[i][j] : 0.0;
            const double q_down = (j > fict_y && is_fluid_cell(i, j - 1, Nx, Ny, fict_x, fict_y, solid_mask)) ? qy[i][j - 1] : 0.0;
            const double q_up = (j < Ny - fict_y - 1 && is_fluid_cell(i, j + 1, Nx, Ny, fict_x, fict_y, solid_mask)) ? qy[i][j] : 0.0;
            const double p_face_left = 0.5 * (p[i][j] + p_left);
            const double p_face_right = 0.5 * (p[i][j] + p_right);
            const double p_face_down = 0.5 * (p[i][j] + p_down);
            const double p_face_up = 0.5 * (p[i][j] + p_up);
            u_tilde[i][j] = u[i][j] - (dt / std::max(rho[i][j], EPS) / dx)
                * ((p_face_right + q_right) - (p_face_left + q_left));
            v_tilde[i][j] = v[i][j] - (dt / std::max(rho[i][j], EPS) / dy)
                * ((p_face_up + q_up) - (p_face_down + q_down));
        }
    }

    ScalarField rho_tilde = rho;
    ScalarField I_tilde = I;
    ScalarField E = make_scalar_field(Nx, Ny, 0.0);

    for (int i = phase_i_begin; i < phase_i_end; ++i) {
        for (int j = phase_j_begin; j < phase_j_end; ++j) {
            if (!is_fluid_cell(i, j, Nx, Ny, fict_x, fict_y, solid_mask)) {
                continue;
            }

            const double u_left_old = sample_normal_velocity(u, i, j, i - 1, j,
                                                             Nx, Ny, fict_x, fict_y,
                                                             left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                                                             solid_mask, 0);
            const double u_right_old = sample_normal_velocity(u, i, j, i + 1, j,
                                                              Nx, Ny, fict_x, fict_y,
                                                              left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                                                              solid_mask, 0);
            const double u_left_new = sample_normal_velocity(u_tilde, i, j, i - 1, j,
                                                             Nx, Ny, fict_x, fict_y,
                                                             left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                                                             solid_mask, 0);
            const double u_right_new = sample_normal_velocity(u_tilde, i, j, i + 1, j,
                                                              Nx, Ny, fict_x, fict_y,
                                                              left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                                                              solid_mask, 0);
            const double v_down_old = sample_normal_velocity(v, i, j, i, j - 1,
                                                             Nx, Ny, fict_x, fict_y,
                                                             left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                                                             solid_mask, 1);
            const double v_up_old = sample_normal_velocity(v, i, j, i, j + 1,
                                                           Nx, Ny, fict_x, fict_y,
                                                           left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                                                           solid_mask, 1);
            const double v_down_new = sample_normal_velocity(v_tilde, i, j, i, j - 1,
                                                             Nx, Ny, fict_x, fict_y,
                                                             left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                                                             solid_mask, 1);
            const double v_up_new = sample_normal_velocity(v_tilde, i, j, i, j + 1,
                                                           Nx, Ny, fict_x, fict_y,
                                                           left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                                                           solid_mask, 1);

            const double U1 = 0.5 * (u_left_old + u_left_new);
            const double U2 = 0.5 * (u_right_old + u_right_new);
            const double V1 = 0.5 * (v_down_old + v_down_new);
            const double V2 = 0.5 * (v_up_old + v_up_new);
            const double T3 = 0.5 * (u[i][j] + u_tilde[i][j]);
            const double T1 = 0.5 * (v[i][j] + v_tilde[i][j]);

            const double q_left = (i > fict_x && is_fluid_cell(i - 1, j, Nx, Ny, fict_x, fict_y, solid_mask)) ? qx[i - 1][j] : 0.0;
            const double q_right = (i < Nx - fict_x - 1 && is_fluid_cell(i + 1, j, Nx, Ny, fict_x, fict_y, solid_mask)) ? qx[i][j] : 0.0;
            const double q_down = (j > fict_y && is_fluid_cell(i, j - 1, Nx, Ny, fict_x, fict_y, solid_mask)) ? qy[i][j - 1] : 0.0;
            const double q_up = (j < Ny - fict_y - 1 && is_fluid_cell(i, j + 1, Nx, Ny, fict_x, fict_y, solid_mask)) ? qy[i][j] : 0.0;

            rho_tilde[i][j] = std::max(rho[i][j] - rho[i][j] * dt
                * ((U2 - U1) / (2.0 * dx) + (V2 - V1) / (2.0 * dy)), EPS);

            const double zip_term
                = (p[i][j] / dx) * (U2 - U1)
                + (q_right / dx) * (U2 - T3)
                + (q_left / dx) * (T3 - U1)
                + (p[i][j] / dy) * (V2 - V1)
                + (q_up / dy) * (V2 - T1)
                + (q_down / dy) * (T1 - V1);

            I_tilde[i][j] = std::max(I[i][j] - (dt / (4.0 * std::max(rho[i][j], EPS))) * zip_term, EPS);
            E[i][j] = I_tilde[i][j] + 0.5 * (u_tilde[i][j] * u_tilde[i][j] + v_tilde[i][j] * v_tilde[i][j]);
        }
    }

    ScalarField DM = make_scalar_field(Nx, Ny, 0.0);
    ScalarField DE = make_scalar_field(Nx, Ny, 0.0);
    ScalarField DW = make_scalar_field(Nx, Ny, 0.0);
    ScalarField DPU = make_scalar_field(Nx, Ny, 0.0);
    ScalarField DPV = make_scalar_field(Nx, Ny, 0.0);

    for (int iface = i_start; iface <= i_end; ++iface) {
        for (int j = j_start; j < j_end; ++j) {
            if (j < phys_j_begin || j >= phys_j_end) {
                continue;
            }
            const int il = iface - 1;
            const int ir = iface;
            FaceState left = make_vertical_side(il, j, Nx, Ny, fict_x, fict_y,
                                                left_bc_code, right_bc_code,
                                                solid_mask, rho_tilde, u_tilde, v_tilde, E, w_reacted);
            FaceState right = make_vertical_side(ir, j, Nx, Ny, fict_x, fict_y,
                                                 left_bc_code, right_bc_code,
                                                 solid_mask, rho_tilde, u_tilde, v_tilde, E, w_reacted);
            if (left.blocked || right.blocked) {
                continue;
            }

            const double alpha = transport_fraction_x(left, right, u_tilde, j, fict_x, right_bc_code, dt, dx);
            const FaceState& donor = (alpha >= 0.0) ? left : right;
            const FaceState& acceptor = (alpha >= 0.0) ? right : left;
            if (!donor.physical) {
                continue;
            }

            const double alpha_abs = std::abs(alpha);
            const double amount = donor.rho * alpha_abs * radial_geom_factor(donor.i, fict_x, alpha_abs, config.slab);
            const double transported_w = corrected_transport_fraction(w_reacted,
                                                                      donor.i, donor.j,
                                                                      acceptor.i, acceptor.j,
                                                                      Nx, Ny, fict_x, fict_y,
                                                                      solid_mask, config);
            accumulate_transport(DM, DE, DW, DPU, DPV, donor, acceptor, amount, transported_w);
        }
    }

    for (int i = i_start; i < i_end; ++i) {
        for (int jface = j_start; jface <= j_end; ++jface) {
            if (i < phys_i_begin || i >= phys_i_end) {
                continue;
            }
            const int jd = jface - 1;
            const int ju = jface;
            FaceState down = make_horizontal_side(i, jd, Nx, Ny, fict_x, fict_y,
                                                  up_bc_code, down_bc_code,
                                                  solid_mask, rho_tilde, u_tilde, v_tilde, E, w_reacted);
            FaceState up = make_horizontal_side(i, ju, Nx, Ny, fict_x, fict_y,
                                                up_bc_code, down_bc_code,
                                                solid_mask, rho_tilde, u_tilde, v_tilde, E, w_reacted);
            if (down.blocked || up.blocked) {
                continue;
            }

            const double beta = transport_fraction_y(down, up, v_tilde, i, up_bc_code, dt, dy);
            const FaceState& donor = (beta >= 0.0) ? down : up;
            const FaceState& acceptor = (beta >= 0.0) ? up : down;
            if (!donor.physical) {
                continue;
            }

            const double amount = donor.rho * std::abs(beta);
            const double transported_w = corrected_transport_fraction(w_reacted,
                                                                      donor.i, donor.j,
                                                                      acceptor.i, acceptor.j,
                                                                      Nx, Ny, fict_x, fict_y,
                                                                      solid_mask, config);
            accumulate_transport(DM, DE, DW, DPU, DPV, donor, acceptor, amount, transported_w);
        }
    }

    for (int i = i_start; i < i_end; ++i) {
        for (int j = j_start; j < j_end; ++j) {
            if (!is_physical_cell(i, j, Nx, Ny, fict_x, fict_y)) {
                continue;
            }
            if (solid_mask[i][j]) {
                w_next[i][j] = 0.0;
                continue;
            }

            double rho_floor = EPS;
            if (rho_ref[i][j] > EPS && config.min_density_factor > 0.0) {
                rho_floor = std::max(rho_floor, config.min_density_factor * rho_ref[i][j]);
            }

            const double rho_base = rho[i][j];
            const double rho_new = std::max(rho_base + DM[i][j], rho_floor);
            const double mom_x = rho_base * u_tilde[i][j] + DPU[i][j];
            const double mom_y = rho_base * v_tilde[i][j] + DPV[i][j];
            const double u_new = mom_x / rho_new;
            const double v_new = mom_y / rho_new;
            const double E_new = (rho_base * E[i][j] + DE[i][j]) / rho_new;
            double I_new = E_new - 0.5 * (u_new * u_new + v_new * v_new);
            if (I_new < EPS) {
                I_new = EPS;
            }
            const double p_new = std::max((g - 1.0) * rho_new * I_new, EPS);
            const double w_new = clamp_unit((rho_base * w_reacted[i][j] + DW[i][j]) / rho_new);

            u_next[i][j] = make_cons(rho_new, u_new, v_new, p_new, g);
            enforce_physical_state(u_next[i][j], g);
            w_next[i][j] = w_new;
        }
    }
}
