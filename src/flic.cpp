#include "flic.h"

#include <algorithm>
#include <vector>

#include "grid.h"
#include "solver.h"

namespace {

using Field3 = std::vector<std::vector<std::vector<double>>>;

void apply_boundaries(Field3& u, int Nx, int Ny, int fict_x, int fict_y,
                      int left_bc_code, int right_bc_code, int up_bc_code, int down_bc_code, double g) {
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < fict_x; ++k) {
            u[k][j] = boundary(u[fict_x][j], left_bc_code, g, 0);
            u[Nx - 1 - k][j] = boundary(u[Nx - fict_x - 1][j], right_bc_code, g, 0);
        }
    }
    for (int i = 0; i < Nx; ++i) {
        for (int k = 0; k < fict_y; ++k) {
            u[i][k] = boundary(u[i][fict_y], down_bc_code, g, 1);
            u[i][Ny - 1 - k] = boundary(u[i][Ny - fict_y - 1], up_bc_code, g, 1);
        }
    }
}

void apply_boundaries_local(Field3& u, int Nx, int Ny, int fict_x, int fict_y,
                            int left_bc_code, int right_bc_code, int up_bc_code, int down_bc_code,
                            int i_start, int i_end, int j_start, int j_end,
                            int left_rank, int right_rank, int down_rank, int up_rank,
                            double g) {
    apply_physical_boundaries_local(u,
                                    Nx, Ny, fict_x, fict_y,
                                    i_start, i_end, j_start, j_end,
                                    left_rank, right_rank, down_rank, up_rank,
                                    left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                                    g);
}

void lagrangian_pressure_x(const Field3& in, Field3& out, double dt, double dx, double g,
                           int Nx, int Ny, int fict_x, int fict_y) {
    out = in;
    for (int i = fict_x; i < Nx - fict_x; ++i) {
        for (int j = fict_y; j < Ny - fict_y; ++j) {
            const std::vector<double> prim_l = cons_to_noncons(in[i - 1][j], g);
            const std::vector<double> prim_c = cons_to_noncons(in[i][j], g);
            const std::vector<double> prim_r = cons_to_noncons(in[i + 1][j], g);

            const double rho = std::max(prim_c[RHO], 1e-8);
            const double u = prim_c[U];
            const double v = prim_c[V];
            const double E_spec = in[i][j][e] / rho;

            const double p_ip = 0.5 * (prim_c[P] + prim_r[P]);
            const double p_im = 0.5 * (prim_l[P] + prim_c[P]);
            const double u_ip = 0.5 * (prim_c[U] + prim_r[U]);
            const double u_im = 0.5 * (prim_l[U] + prim_c[U]);

            const double u_tilde = u - (dt / rho) * (p_ip - p_im) / dx;
            const double E_tilde_spec = E_spec - (dt / rho) * ((p_ip * u_ip) - (p_im * u_im)) / dx;

            out[i][j][r] = rho;
            out[i][j][ru] = rho * u_tilde;
            out[i][j][rv] = rho * v;
            out[i][j][e] = rho * E_tilde_spec;
            enforce_physical_state(out[i][j], g);
        }
    }
}

void lagrangian_pressure_x_local(const Field3& in, Field3& out, double dt, double dx, double g,
                                 int i_start, int i_end, int j_start, int j_end) {
    out = in;
    for (int i = i_start; i < i_end; ++i) {
        for (int j = j_start; j < j_end; ++j) {
            const std::vector<double> prim_l = cons_to_noncons(in[i - 1][j], g);
            const std::vector<double> prim_c = cons_to_noncons(in[i][j], g);
            const std::vector<double> prim_r = cons_to_noncons(in[i + 1][j], g);

            const double rho = std::max(prim_c[RHO], 1e-8);
            const double u = prim_c[U];
            const double v = prim_c[V];
            const double E_spec = in[i][j][e] / rho;

            const double p_ip = 0.5 * (prim_c[P] + prim_r[P]);
            const double p_im = 0.5 * (prim_l[P] + prim_c[P]);
            const double u_ip = 0.5 * (prim_c[U] + prim_r[U]);
            const double u_im = 0.5 * (prim_l[U] + prim_c[U]);

            const double u_tilde = u - (dt / rho) * (p_ip - p_im) / dx;
            const double E_tilde_spec = E_spec - (dt / rho) * ((p_ip * u_ip) - (p_im * u_im)) / dx;

            out[i][j][r] = rho;
            out[i][j][ru] = rho * u_tilde;
            out[i][j][rv] = rho * v;
            out[i][j][e] = rho * E_tilde_spec;
            enforce_physical_state(out[i][j], g);
        }
    }
}

void lagrangian_pressure_y(const Field3& in, Field3& out, double dt, double dy, double g,
                           int Nx, int Ny, int fict_x, int fict_y) {
    out = in;
    for (int i = fict_x; i < Nx - fict_x; ++i) {
        for (int j = fict_y; j < Ny - fict_y; ++j) {
            const std::vector<double> prim_d = cons_to_noncons(in[i][j - 1], g);
            const std::vector<double> prim_c = cons_to_noncons(in[i][j], g);
            const std::vector<double> prim_u = cons_to_noncons(in[i][j + 1], g);

            const double rho = std::max(prim_c[RHO], 1e-8);
            const double u = prim_c[U];
            const double v = prim_c[V];
            const double E_spec = in[i][j][e] / rho;

            const double p_jp = 0.5 * (prim_c[P] + prim_u[P]);
            const double p_jm = 0.5 * (prim_d[P] + prim_c[P]);
            const double v_jp = 0.5 * (prim_c[V] + prim_u[V]);
            const double v_jm = 0.5 * (prim_d[V] + prim_c[V]);

            const double v_tilde = v - (dt / rho) * (p_jp - p_jm) / dy;
            const double E_tilde_spec = E_spec - (dt / rho) * ((p_jp * v_jp) - (p_jm * v_jm)) / dy;

            out[i][j][r] = rho;
            out[i][j][ru] = rho * u;
            out[i][j][rv] = rho * v_tilde;
            out[i][j][e] = rho * E_tilde_spec;
            enforce_physical_state(out[i][j], g);
        }
    }
}

void lagrangian_pressure_y_local(const Field3& in, Field3& out, double dt, double dy, double g,
                                 int i_start, int i_end, int j_start, int j_end) {
    out = in;
    for (int i = i_start; i < i_end; ++i) {
        for (int j = j_start; j < j_end; ++j) {
            const std::vector<double> prim_d = cons_to_noncons(in[i][j - 1], g);
            const std::vector<double> prim_c = cons_to_noncons(in[i][j], g);
            const std::vector<double> prim_u = cons_to_noncons(in[i][j + 1], g);

            const double rho = std::max(prim_c[RHO], 1e-8);
            const double u = prim_c[U];
            const double v = prim_c[V];
            const double E_spec = in[i][j][e] / rho;

            const double p_jp = 0.5 * (prim_c[P] + prim_u[P]);
            const double p_jm = 0.5 * (prim_d[P] + prim_c[P]);
            const double v_jp = 0.5 * (prim_c[V] + prim_u[V]);
            const double v_jm = 0.5 * (prim_d[V] + prim_c[V]);

            const double v_tilde = v - (dt / rho) * (p_jp - p_jm) / dy;
            const double E_tilde_spec = E_spec - (dt / rho) * ((p_jp * v_jp) - (p_jm * v_jm)) / dy;

            out[i][j][r] = rho;
            out[i][j][ru] = rho * u;
            out[i][j][rv] = rho * v_tilde;
            out[i][j][e] = rho * E_tilde_spec;
            enforce_physical_state(out[i][j], g);
        }
    }
}

void advection_x(const Field3& in, Field3& out, double dt, double dx, double g,
                 int Nx, int Ny, int fict_x, int fict_y) {
    out = in;
    for (int i = fict_x; i < Nx - fict_x; ++i) {
        for (int j = fict_y; j < Ny - fict_y; ++j) {
            std::vector<double> flux_left(M, 0.0);
            std::vector<double> flux_right(M, 0.0);

            {
                const std::vector<double> pl = cons_to_noncons(in[i - 1][j], g);
                const std::vector<double> pc = cons_to_noncons(in[i][j], g);
                const double u_face = 0.5 * (pl[U] + pc[U]);
                const std::vector<double>& donor = (u_face > 0.0) ? in[i - 1][j] : in[i][j];
                for (int k = 0; k < M; ++k) {
                    flux_left[k] = u_face * donor[k];
                }
            }
            {
                const std::vector<double> pc = cons_to_noncons(in[i][j], g);
                const std::vector<double> pr = cons_to_noncons(in[i + 1][j], g);
                const double u_face = 0.5 * (pc[U] + pr[U]);
                const std::vector<double>& donor = (u_face > 0.0) ? in[i][j] : in[i + 1][j];
                for (int k = 0; k < M; ++k) {
                    flux_right[k] = u_face * donor[k];
                }
            }

            for (int k = 0; k < M; ++k) {
                out[i][j][k] = in[i][j][k] - (dt / dx) * (flux_right[k] - flux_left[k]);
            }
            enforce_physical_state(out[i][j], g);
        }
    }
}

void advection_x_local(const Field3& in, Field3& out, double dt, double dx, double g,
                       int i_start, int i_end, int j_start, int j_end) {
    out = in;
    for (int i = i_start; i < i_end; ++i) {
        for (int j = j_start; j < j_end; ++j) {
            std::vector<double> flux_left(M, 0.0);
            std::vector<double> flux_right(M, 0.0);

            {
                const std::vector<double> pl = cons_to_noncons(in[i - 1][j], g);
                const std::vector<double> pc = cons_to_noncons(in[i][j], g);
                const double u_face = 0.5 * (pl[U] + pc[U]);
                const std::vector<double>& donor = (u_face > 0.0) ? in[i - 1][j] : in[i][j];
                for (int k = 0; k < M; ++k) {
                    flux_left[k] = u_face * donor[k];
                }
            }
            {
                const std::vector<double> pc = cons_to_noncons(in[i][j], g);
                const std::vector<double> pr = cons_to_noncons(in[i + 1][j], g);
                const double u_face = 0.5 * (pc[U] + pr[U]);
                const std::vector<double>& donor = (u_face > 0.0) ? in[i][j] : in[i + 1][j];
                for (int k = 0; k < M; ++k) {
                    flux_right[k] = u_face * donor[k];
                }
            }

            for (int k = 0; k < M; ++k) {
                out[i][j][k] = in[i][j][k] - (dt / dx) * (flux_right[k] - flux_left[k]);
            }
            enforce_physical_state(out[i][j], g);
        }
    }
}

void advection_y(const Field3& in, Field3& out, double dt, double dy, double g,
                 int Nx, int Ny, int fict_x, int fict_y) {
    out = in;
    for (int i = fict_x; i < Nx - fict_x; ++i) {
        for (int j = fict_y; j < Ny - fict_y; ++j) {
            std::vector<double> flux_down(M, 0.0);
            std::vector<double> flux_up(M, 0.0);

            {
                const std::vector<double> pd = cons_to_noncons(in[i][j - 1], g);
                const std::vector<double> pc = cons_to_noncons(in[i][j], g);
                const double v_face = 0.5 * (pd[V] + pc[V]);
                const std::vector<double>& donor = (v_face > 0.0) ? in[i][j - 1] : in[i][j];
                for (int k = 0; k < M; ++k) {
                    flux_down[k] = v_face * donor[k];
                }
            }
            {
                const std::vector<double> pc = cons_to_noncons(in[i][j], g);
                const std::vector<double> pu = cons_to_noncons(in[i][j + 1], g);
                const double v_face = 0.5 * (pc[V] + pu[V]);
                const std::vector<double>& donor = (v_face > 0.0) ? in[i][j] : in[i][j + 1];
                for (int k = 0; k < M; ++k) {
                    flux_up[k] = v_face * donor[k];
                }
            }

            for (int k = 0; k < M; ++k) {
                out[i][j][k] = in[i][j][k] - (dt / dy) * (flux_up[k] - flux_down[k]);
            }
            enforce_physical_state(out[i][j], g);
        }
    }
}

void advection_y_local(const Field3& in, Field3& out, double dt, double dy, double g,
                       int i_start, int i_end, int j_start, int j_end) {
    out = in;
    for (int i = i_start; i < i_end; ++i) {
        for (int j = j_start; j < j_end; ++j) {
            std::vector<double> flux_down(M, 0.0);
            std::vector<double> flux_up(M, 0.0);

            {
                const std::vector<double> pd = cons_to_noncons(in[i][j - 1], g);
                const std::vector<double> pc = cons_to_noncons(in[i][j], g);
                const double v_face = 0.5 * (pd[V] + pc[V]);
                const std::vector<double>& donor = (v_face > 0.0) ? in[i][j - 1] : in[i][j];
                for (int k = 0; k < M; ++k) {
                    flux_down[k] = v_face * donor[k];
                }
            }
            {
                const std::vector<double> pc = cons_to_noncons(in[i][j], g);
                const std::vector<double> pu = cons_to_noncons(in[i][j + 1], g);
                const double v_face = 0.5 * (pc[V] + pu[V]);
                const std::vector<double>& donor = (v_face > 0.0) ? in[i][j] : in[i][j + 1];
                for (int k = 0; k < M; ++k) {
                    flux_up[k] = v_face * donor[k];
                }
            }

            for (int k = 0; k < M; ++k) {
                out[i][j][k] = in[i][j][k] - (dt / dy) * (flux_up[k] - flux_down[k]);
            }
            enforce_physical_state(out[i][j], g);
        }
    }
}

} // namespace

void flic_step_2d(const Field3& u_prev, Field3& u_next, double dt, double dx, double dy, double g,
                  int Nx, int Ny, int fict_x, int fict_y,
                  int left_bc_code, int right_bc_code, int up_bc_code, int down_bc_code) {
    Field3 stage0 = u_prev;
    Field3 stage1 = u_prev;
    Field3 stage2 = u_prev;
    Field3 stage3 = u_prev;
    Field3 stage4 = u_prev;

    const double dt_half = 0.5 * dt;

    apply_boundaries(stage0, Nx, Ny, fict_x, fict_y, left_bc_code, right_bc_code, up_bc_code, down_bc_code, g);
    lagrangian_pressure_x(stage0, stage1, dt_half, dx, g, Nx, Ny, fict_x, fict_y);
    apply_boundaries(stage1, Nx, Ny, fict_x, fict_y, left_bc_code, right_bc_code, up_bc_code, down_bc_code, g);
    advection_x(stage1, stage2, dt_half, dx, g, Nx, Ny, fict_x, fict_y);
    apply_boundaries(stage2, Nx, Ny, fict_x, fict_y, left_bc_code, right_bc_code, up_bc_code, down_bc_code, g);

    lagrangian_pressure_y(stage2, stage3, dt, dy, g, Nx, Ny, fict_x, fict_y);
    apply_boundaries(stage3, Nx, Ny, fict_x, fict_y, left_bc_code, right_bc_code, up_bc_code, down_bc_code, g);
    advection_y(stage3, stage4, dt, dy, g, Nx, Ny, fict_x, fict_y);
    apply_boundaries(stage4, Nx, Ny, fict_x, fict_y, left_bc_code, right_bc_code, up_bc_code, down_bc_code, g);

    lagrangian_pressure_x(stage4, stage1, dt_half, dx, g, Nx, Ny, fict_x, fict_y);
    apply_boundaries(stage1, Nx, Ny, fict_x, fict_y, left_bc_code, right_bc_code, up_bc_code, down_bc_code, g);
    advection_x(stage1, u_next, dt_half, dx, g, Nx, Ny, fict_x, fict_y);
    apply_boundaries(u_next, Nx, Ny, fict_x, fict_y, left_bc_code, right_bc_code, up_bc_code, down_bc_code, g);
}

void flic_step_2d_local(const Field3& u_prev,
                        Field3& u_next,
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
                        MPI_Comm cart_comm,
                        int i_start,
                        int i_end,
                        int j_start,
                        int j_end,
                        int left_rank,
                        int right_rank,
                        int down_rank,
                        int up_rank) {
    Field3 stage0 = u_prev;
    Field3 stage1 = u_prev;
    Field3 stage2 = u_prev;
    Field3 stage3 = u_prev;
    Field3 stage4 = u_prev;

    const double dt_half = 0.5 * dt;

    exchange_halos_global(stage0, cart_comm, fict_x, fict_y,
                          i_start, i_end, j_start, j_end,
                          left_rank, right_rank, down_rank, up_rank);
    apply_boundaries_local(stage0, Nx, Ny, fict_x, fict_y,
                           left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                           i_start, i_end, j_start, j_end,
                           left_rank, right_rank, down_rank, up_rank,
                           g);
    lagrangian_pressure_x_local(stage0, stage1, dt_half, dx, g, i_start, i_end, j_start, j_end);

    exchange_halos_global(stage1, cart_comm, fict_x, fict_y,
                          i_start, i_end, j_start, j_end,
                          left_rank, right_rank, down_rank, up_rank);
    apply_boundaries_local(stage1, Nx, Ny, fict_x, fict_y,
                           left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                           i_start, i_end, j_start, j_end,
                           left_rank, right_rank, down_rank, up_rank,
                           g);
    advection_x_local(stage1, stage2, dt_half, dx, g, i_start, i_end, j_start, j_end);

    exchange_halos_global(stage2, cart_comm, fict_x, fict_y,
                          i_start, i_end, j_start, j_end,
                          left_rank, right_rank, down_rank, up_rank);
    apply_boundaries_local(stage2, Nx, Ny, fict_x, fict_y,
                           left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                           i_start, i_end, j_start, j_end,
                           left_rank, right_rank, down_rank, up_rank,
                           g);
    lagrangian_pressure_y_local(stage2, stage3, dt, dy, g, i_start, i_end, j_start, j_end);

    exchange_halos_global(stage3, cart_comm, fict_x, fict_y,
                          i_start, i_end, j_start, j_end,
                          left_rank, right_rank, down_rank, up_rank);
    apply_boundaries_local(stage3, Nx, Ny, fict_x, fict_y,
                           left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                           i_start, i_end, j_start, j_end,
                           left_rank, right_rank, down_rank, up_rank,
                           g);
    advection_y_local(stage3, stage4, dt, dy, g, i_start, i_end, j_start, j_end);

    exchange_halos_global(stage4, cart_comm, fict_x, fict_y,
                          i_start, i_end, j_start, j_end,
                          left_rank, right_rank, down_rank, up_rank);
    apply_boundaries_local(stage4, Nx, Ny, fict_x, fict_y,
                           left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                           i_start, i_end, j_start, j_end,
                           left_rank, right_rank, down_rank, up_rank,
                           g);
    lagrangian_pressure_x_local(stage4, stage1, dt_half, dx, g, i_start, i_end, j_start, j_end);

    exchange_halos_global(stage1, cart_comm, fict_x, fict_y,
                          i_start, i_end, j_start, j_end,
                          left_rank, right_rank, down_rank, up_rank);
    apply_boundaries_local(stage1, Nx, Ny, fict_x, fict_y,
                           left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                           i_start, i_end, j_start, j_end,
                           left_rank, right_rank, down_rank, up_rank,
                           g);
    advection_x_local(stage1, u_next, dt_half, dx, g, i_start, i_end, j_start, j_end);

    apply_boundaries_local(u_next, Nx, Ny, fict_x, fict_y,
                           left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                           i_start, i_end, j_start, j_end,
                           left_rank, right_rank, down_rank, up_rank,
                           g);
}
