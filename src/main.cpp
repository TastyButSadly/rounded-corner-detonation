#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "mpi_compat.h"

#include "./parser.h"
#include "./utils.h"
#include "./solver.h"
#include "./flic.h"
#include "./mader.h"
#include "./point.h"
#include "./grid.h"
#include <script.h>

int main(int argc, char** argv) {
// Изменено на единый файл input.ini
    MPI_Init(&argc, &argv);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    auto abort_with_message = [&](const std::string& message) -> int {
        if (rank == 0) {
            std::cerr << message << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    };

    auto parse_bool = [](const std::string& value, bool default_value = false) -> bool {
        if (value == "1" || value == "true" || value == "True" || value == "TRUE") {
            return true;
        }
        if (value == "0" || value == "false" || value == "False" || value == "FALSE") {
            return false;
        }
        return default_value;
    };

    std::string system_ini;
    if (argc > 1) {
        system_ini = argv[1];
    } else if (std::filesystem::exists("configs/input.ini")) {
        system_ini = "configs/input.ini";
    } else {
        system_ini = "../configs/input.ini";
    }

    // чтение системных параметров из объединенного файла
    IniParser sys(system_ini);
    std::string case_name = "mader_2de_case";
    try { case_name = sys.getString("case_name"); } catch (...) {}
    double x_min = sys.getDouble("x_min");
    double x_max = sys.getDouble("x_max");
    double y_min = sys.getDouble("y_min");
    double y_max = sys.getDouble("y_max");
    int Nx = sys.getInt("Nx");
    int Ny = sys.getInt("Ny");
    int fict_x = sys.getInt("fict_x");
    int fict_y = sys.getInt("fict_y");
    double tmax = sys.getDouble("tmax");
    double cfl = sys.getDouble("cfl");
    double g = sys.getDouble("g");
    double dt_out = sys.getDouble("dt_out");
    int step_num = sys.getInt("step_num");
    int px = 0;
    int py = 0;
    try { px = sys.getInt("px"); } catch (...) {}
    try { py = sys.getInt("py"); } catch (...) {}
    std::string boundary_type_left = sys.getString("boundary_type_left");
    std::string boundary_type_right = sys.getString("boundary_type_right");
    std::string boundary_type_up = sys.getString("boundary_type_up");
    std::string boundary_type_down = sys.getString("boundary_type_down");
    int equation_type = sys.getInt("equation_type");
    std::string analytic_axis = "x";
    int analytic_profile_index = -1;
    try { analytic_axis = sys.getString("analytic_axis"); } catch (...) {}
    try { analytic_profile_index = sys.getInt("analytic_profile_index"); } catch (...) {}
    double default_w = 0.0;
    try { default_w = sys.getDouble("default_w"); } catch (...) {}
    MaderConfig mader_config;
    try { mader_config.visc = sys.getDouble("mader_visc"); } catch (...) {}
    try {
        std::string step_enabled = sys.getString("step_enabled");
        mader_config.step_enabled = parse_bool(step_enabled, false);
    } catch (...) {}
    try {
        std::string slab_flag = sys.getString("mader_slab");
        mader_config.slab = parse_bool(slab_flag, true);
    } catch (...) {}
    try {
        std::string reactive_flag = sys.getString("reaction_enabled");
        mader_config.reaction_enabled = parse_bool(reactive_flag, false);
    } catch (...) {}
    try {
        std::string shargatov_flag = sys.getString("shargatov_correction");
        mader_config.shargatov_correction = parse_bool(shargatov_flag, false);
    } catch (...) {}
    try { mader_config.step_x_end = sys.getDouble("step_x_end"); } catch (...) {}
    try { mader_config.step_y_end = sys.getDouble("step_y_end"); } catch (...) {}
    try { mader_config.corner_radius = sys.getDouble("corner_radius"); } catch (...) {}
    try { mader_config.gas_constant = sys.getDouble("gas_constant"); } catch (...) {}
    try { mader_config.reaction_rate = sys.getDouble("reaction_rate"); } catch (...) {}
    try { mader_config.reaction_activation_energy = sys.getDouble("reaction_activation_energy"); } catch (...) {}
    try { mader_config.reaction_heat_release = sys.getDouble("reaction_heat_release"); } catch (...) {}
    try { mader_config.min_temperature = sys.getDouble("min_temperature"); } catch (...) {}
    try { mader_config.gasw_threshold = sys.getDouble("gasw_threshold"); } catch (...) {}
    try { mader_config.reaction_delay_steps = sys.getInt("reaction_delay_steps"); } catch (...) {}
    try { mader_config.min_density_factor = sys.getDouble("min_density_factor"); } catch (...) {}
    const int analytic_axis_code = (analytic_axis == "y" || analytic_axis == "Y") ? 1 : 0;

    if (px < 0 || py < 0) {
        return abort_with_message("px and py must be non-negative.");
    }
    if (mader_config.step_enabled) {
        if (!std::isfinite(mader_config.corner_radius) || mader_config.corner_radius < 0.0) {
            return abort_with_message("corner_radius must be non-negative.");
        }
        const double max_corner_radius = std::min(mader_config.step_x_end - x_min,
                                                  mader_config.step_y_end - y_min);
        if (mader_config.corner_radius > max_corner_radius + 1e-12) {
            return abort_with_message("corner_radius is too large for the solid step geometry.");
        }
    }
    if (px > 0 && py > 0 && px * py != size) {
        return abort_with_message("Configured px * py must match the MPI world size.");
    }

    int dims[2] = {px, py};
    MPI_Dims_create(size, 2, dims);
    if ((px > 0 && dims[0] != px) || (py > 0 && dims[1] != py) || dims[0] * dims[1] != size) {
        return abort_with_message("MPI_Dims_create could not satisfy the requested px/py decomposition.");
    }

    int periods[2] = {0, 0};
    MPI_Comm cart_comm = MPI_COMM_NULL;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
    if (cart_comm == MPI_COMM_NULL) {
        return abort_with_message("Failed to create MPI Cartesian communicator.");
    }

    int coords[2] = {0, 0};
    MPI_Cart_coords(cart_comm, rank, 2, coords);

    int left_rank = MPI_PROC_NULL;
    int right_rank = MPI_PROC_NULL;
    int down_rank = MPI_PROC_NULL;
    int up_rank = MPI_PROC_NULL;
    MPI_Cart_shift(cart_comm, 0, 1, &left_rank, &right_rank);
    MPI_Cart_shift(cart_comm, 1, 1, &down_rank, &up_rank);

    int i_start_phys = 0;
    int i_end_phys = 0;
    int j_start_phys = 0;
    int j_end_phys = 0;
    get_subdomain_bounds(Nx, dims[0], coords[0], i_start_phys, i_end_phys);
    get_subdomain_bounds(Ny, dims[1], coords[1], j_start_phys, j_end_phys);

    int i_start = i_start_phys + fict_x;
    int i_end = i_end_phys + fict_x;
    int j_start = j_start_phys + fict_y;
    int j_end = j_end_phys + fict_y;

    double dx = (x_max - x_min) / Nx;
    double dy = (y_max - y_min) / Ny;
    int Nx_with_fict_cells = Nx + 2 * fict_x;
    int Ny_with_fict_cells = Ny + 2 * fict_y;

    // инициализация сетки
    std::vector<std::vector<std::vector<double>>> u_prev(
        Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> u_next(
        Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    ScalarField w_prev = make_scalar_field(Nx_with_fict_cells, Ny_with_fict_cells, 0.0);
    ScalarField w_next = make_scalar_field(Nx_with_fict_cells, Ny_with_fict_cells, 0.0);
    ScalarField rho_ref = make_scalar_field(Nx_with_fict_cells, Ny_with_fict_cells, 0.0);
    if (equation_type == 9) {
        parseAndInitializeMader(u_prev, w_prev, system_ini, fict_x, fict_y, g, default_w);
        parseAndInitializeMader(u_next, w_next, system_ini, fict_x, fict_y, g, default_w);
    } else {
        parseAndInitialize(u_prev, system_ini, fict_x, fict_y, g);
        parseAndInitialize(u_next, system_ini, fict_x, fict_y, g);
    }

    Mask2D solid_mask(Nx_with_fict_cells, std::vector<unsigned char>(Ny_with_fict_cells, 0));
    if (equation_type == 9) {
        solid_mask = build_step_mask(x_min, y_min, dx, dy,
                                     Nx_with_fict_cells, Ny_with_fict_cells,
                                     fict_x, fict_y, mader_config);
        apply_step_mask_scalar(w_prev, solid_mask, 0.0);
        apply_step_mask_scalar(w_next, solid_mask, 0.0);
        for (int i = fict_x; i < Nx_with_fict_cells - fict_x; ++i) {
            for (int j = fict_y; j < Ny_with_fict_cells - fict_y; ++j) {
                if (solid_mask[i][j]) {
                    continue;
                }
                rho_ref[i][j] = cons_to_noncons(u_prev[i][j], g)[RHO];
            }
        }
    }

    const auto u_initial = u_prev;
    //set_sod_initial_conditions(u_prev, u_next, Nx, x_min, x_max, Ny, y_min, y_max,
      //                        rho_R, u_R, v_R, p_R, rho_L, u_L, v_L, p_L, g, fict_x, fict_y);

    auto parse_boundary_code = [](const std::string& value) -> int {
        try {
            int code = std::stoi(value);
            return (code <= 0) ? 2 : code;
        } catch (...) {
            return 2; // default to outflow
        }
    };

    const int left_bc_code = parse_boundary_code(boundary_type_left);
    const int right_bc_code = parse_boundary_code(boundary_type_right);
    const int up_bc_code = parse_boundary_code(boundary_type_up);
    const int down_bc_code = parse_boundary_code(boundary_type_down);

    std::string steps_dir;
    if (argc > 2) {
        steps_dir = argv[2];
    } else {
        try {
            steps_dir = sys.getString("output_dir");
        } catch (...) {
            std::filesystem::path default_dir = std::filesystem::path("output") / case_name;
            steps_dir = default_dir.string();
        }
    }
    if (rank == 0) {
        std::filesystem::create_directories(steps_dir);
    }
    MPI_Barrier(cart_comm);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    auto write_state_csv = [&](const std::string& path,
                               const auto& field,
                               const ScalarField* mader_w,
                               double time,
                               bool with_exact) {
        std::vector<std::vector<std::vector<double>>> analytic_solution;
        if (with_exact) {
            analytic_solution = compute_analytic_solution_2d(
                x_min, x_max, y_min, y_max, Nx, Ny, fict_x, fict_y,
                time, u_initial, analytic_axis_code, analytic_profile_index, g
            );
        }

        std::ofstream fout(path);
        fout << "x,y,rho,u,v,p,w,temperature,rho_exact,u_exact,v_exact,p_exact,is_solid\n";
        for (int i = fict_x; i < Nx_with_fict_cells - fict_x; ++i) {
            const double x = x_min + (i - fict_x + 0.5) * dx;
            for (int j = fict_y; j < Ny_with_fict_cells - fict_y; ++j) {
                const double y = y_min + (j - fict_y + 0.5) * dy;
                const bool is_solid = solid_mask[i][j] != 0;

                if (is_solid) {
                    fout << x << "," << y << "," << nan << "," << nan << "," << nan << "," << nan << ","
                         << nan << "," << nan << ","
                         << nan << "," << nan << "," << nan << "," << nan << ",1\n";
                    continue;
                }

                const std::vector<double> prim = cons_to_noncons(field[i][j], g);
                const double w_value = (mader_w != nullptr) ? (*mader_w)[i][j] : nan;
                const double temperature = (mader_w != nullptr) ? mader_temperature(prim[RHO], prim[P], mader_config) : nan;
                double rho_exact = nan;
                double u_exact = nan;
                double v_exact = nan;
                double p_exact = nan;

                if (with_exact) {
                    const std::vector<double>& prim_exact = analytic_solution[i][j];
                    rho_exact = prim_exact[0];
                    u_exact = prim_exact[1];
                    v_exact = prim_exact[2];
                    p_exact = prim_exact[3];
                }

                fout << x << "," << y << "," << prim[RHO] << "," << prim[U] << "," << prim[V] << "," << prim[P] << ","
                     << w_value << "," << temperature << ","
                     << rho_exact << "," << u_exact << "," << v_exact << "," << p_exact << ",0\n";
            }
        }
    };

    // сохранение начальных условий
    gather_to_root(u_prev, Nx, Ny, fict_x, fict_y, i_start, i_end, j_start, j_end, rank, size, cart_comm, dims);
    if (equation_type == 9) {
        gather_scalar_to_root(w_prev, Nx, Ny, fict_x, fict_y, i_start, i_end, j_start, j_end, rank, size, cart_comm, dims);
    }
    if (rank == 0) {
        write_state_csv((std::filesystem::path(steps_dir) / "step_0_initial.csv").string(),
                        u_prev,
                        equation_type == 9 ? &w_prev : nullptr,
                        0.0,
                        equation_type != 9);
    }

    // основной цикл по времени
    double curr_time = 0;
    int step = 0;

// Массивы для хранения реконструкции
    // X-direction structures
    std::vector<std::vector<std::vector<double>>> left_face_x(
        Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> right_face_x(
        Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> delta_x(
        Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));

    // Y-direction structures (NEW)
    // left_face_y соответствует "нижней" грани (j-1/2), right_face_y - "верхней" (j+1/2) в локальных координатах
    std::vector<std::vector<std::vector<double>>> left_face_y(
        Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> right_face_y(
        Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> delta_y(
        Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));

    // Массивы реконструкции для Родионова (Corrector step arrays)
    std::vector<std::vector<std::vector<double>>> rec_x_L(
        Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> rec_x_R(
        Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> rec_y_L(
        Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> rec_y_R(
        Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> u_half(
        Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));

    double next_output_time = (dt_out > 0.0) ? dt_out : std::numeric_limits<double>::infinity();
    auto start = std::chrono::high_resolution_clock::now();
    while (curr_time < tmax) {
        exchange_halos_global(u_prev, cart_comm, fict_x, fict_y, i_start, i_end, j_start, j_end,
                              left_rank, right_rank, down_rank, up_rank);
        if (equation_type == 9) {
            exchange_halos_scalar(w_prev, cart_comm, fict_x, fict_y, i_start, i_end, j_start, j_end,
                                  left_rank, right_rank, down_rank, up_rank);
        }

        // 1. Установка граничных условий на текущем слое
        // 1. Граничные условия (Ghost Cells)
        if (equation_type != 9) {
            apply_physical_boundaries_local(u_prev, Nx_with_fict_cells, Ny_with_fict_cells, fict_x, fict_y,
                                            i_start, i_end, j_start, j_end,
                                            left_rank, right_rank, down_rank, up_rank,
                                            left_bc_code, right_bc_code, up_bc_code, down_bc_code, g);
        }

        double local_dt = (equation_type == 9)
            ? calc_time_step_masked_local(u_prev, solid_mask, dx, dy, cfl, g, fict_x, fict_y,
                                          i_start, i_end, j_start, j_end)
            : calc_time_step_local(u_prev, dx, dy, cfl, g, i_start, i_end, j_start, j_end);

        double dt = local_dt;
        MPI_Allreduce(&local_dt, &dt, 1, MPI_DOUBLE, MPI_MIN, cart_comm);
        if ((curr_time + dt) > tmax) {
            dt = tmax - curr_time; // Исправлена логика обрезки шага
        }
        if (!std::isfinite(dt) || dt <= 0.0) {
            return abort_with_message("Failed to compute a positive time step.");
        }

        // --- ПРЕДВАРИТЕЛЬНЫЕ РАСЧЕТЫ (Вне цикла по ячейкам) ---
        const int rec_i_start = std::max(0, i_start - 1);
        const int rec_i_end = std::min(Nx_with_fict_cells - 1, i_end);
        const int rec_j_start = std::max(0, j_start - 1);
        const int rec_j_end = std::min(Ny_with_fict_cells - 1, j_end);

        // TYPE 1: KOLGAN (2D Extension)
        if (equation_type == 1) {
            // 1. Reconstruct along X (for each row j)
            for (int j = rec_j_start; j <= rec_j_end; ++j) {
                for (int i = rec_i_start; i <= rec_i_end; ++i) {
                    if (i == 0) {
                        Kolgan(u_prev[i][j], u_prev[i][j], u_prev[i + 1][j], left_face_x[i][j], right_face_x[i][j], g);
                    } else if (i == Nx_with_fict_cells - 1) {
                        Kolgan(u_prev[i - 1][j], u_prev[i][j],
                               u_prev[i][j], left_face_x[i][j],
                               right_face_x[i][j], g);
                    } else {
                        Kolgan(u_prev[i - 1][j], u_prev[i][j], u_prev[i + 1][j],
                               left_face_x[i][j], right_face_x[i][j], g);
                    }
                }
            }
            // 2. Reconstruct along Y (for each column i)
            for (int i = rec_i_start; i <= rec_i_end; ++i) {
                for (int j = rec_j_start; j <= rec_j_end; ++j) {
                    if (j == 0) {
                        Kolgan(u_prev[i][j], u_prev[i][j], u_prev[i][j + 1], left_face_y[i][j], right_face_y[i][j], g);
                    } else if (j == Ny_with_fict_cells - 1) {
                        Kolgan(u_prev[i][j - 1],
                               u_prev[i][j],
                               u_prev[i][j], left_face_y[i][j], right_face_y[i][j], g);
                    } else {
                        Kolgan(u_prev[i][j - 1], u_prev[i][j], u_prev[i][j + 1], left_face_y[i][j], right_face_y[i][j], g);
                    }
                }
            }
        }

        // TYPE 2: RODIONOV (2D Extension)
        if (equation_type == 2) {
            // 1. Расчет дельт (наклонов) по X
            for (int j = rec_j_start; j <= rec_j_end; ++j) {
                for (int i = rec_i_start; i <= rec_i_end; ++i) {
                    if (i == 0) {
                        Kolgan_for_Rodionov(u_prev[0][j], u_prev[0][j], u_prev[1][j], delta_x[i][j], g);
                    } else if (i == Nx_with_fict_cells - 1) {
                        Kolgan_for_Rodionov(u_prev[Nx_with_fict_cells - 2][j],
                                            u_prev[Nx_with_fict_cells - 1][j], u_prev[Nx_with_fict_cells - 1][j],
                                            delta_x[i][j], g);
                    } else {
                        Kolgan_for_Rodionov(u_prev[i - 1][j], u_prev[i][j], u_prev[i + 1][j], delta_x[i][j], g);
                    }
                }
            }
            // 2. Расчет дельт (наклонов) по Y
            for (int i = rec_i_start; i <= rec_i_end; ++i) {
                for (int j = rec_j_start; j <= rec_j_end; ++j) {
                    if (j == 0) {
                        Kolgan_for_Rodionov(u_prev[i][0], u_prev[i][0], u_prev[i][1], delta_y[i][j], g);
                    } else if (j == Ny_with_fict_cells - 1) {
                        Kolgan_for_Rodionov(u_prev[i][Ny_with_fict_cells - 2],
                                            u_prev[i][Ny_with_fict_cells - 1],
                                            u_prev[i][Ny_with_fict_cells - 1], delta_y[i][j], g);
                    } else {
                        Kolgan_for_Rodionov(u_prev[i][j - 1], u_prev[i][j], u_prev[i][j + 1], delta_y[i][j], g);
                    }
                }
            }

            // 3. Шаг Предиктор (расчет u_half) - теперь локальный MPI-совместимый
            u_half = u_prev;
            Rodionov_local(u_prev, u_half, delta_x, delta_y, dt, dx, dy, g,
                           i_start, i_end, j_start, j_end);
            exchange_halos_global(u_half, cart_comm, fict_x, fict_y, i_start, i_end, j_start, j_end,
                                  left_rank, right_rank, down_rank, up_rank);
            apply_physical_boundaries_local(u_half, Nx_with_fict_cells, Ny_with_fict_cells, fict_x, fict_y,
                                            i_start, i_end, j_start, j_end,
                                            left_rank, right_rank, down_rank, up_rank,
                                            left_bc_code, right_bc_code, up_bc_code, down_bc_code, g);

            // 4. Шаг Корректор: Реконструкция на гранях
            for (int i = rec_i_start; i <= rec_i_end; ++i) {
                for (int j = rec_j_start; j <= rec_j_end; ++j) {
                    for (int var = 0; var < M; ++var) {
                        double u_avg = 0.5 * (u_prev[i][j][var] + u_half[i][j][var]);

                        // По X: rec_x_L это "левое" значение внутри ячейки (для потока слева это right state)
                        rec_x_L[i][j][var] = u_avg - 0.5 * delta_x[i][j][var];
                        rec_x_R[i][j][var] = u_avg + 0.5 * delta_x[i][j][var];

                        // По Y:
                        rec_y_L[i][j][var] = u_avg - 0.5 * delta_y[i][j][var];
                        rec_y_R[i][j][var] = u_avg + 0.5 * delta_y[i][j][var];
                    }
                    enforce_physical_state(rec_x_L[i][j], g);
                    enforce_physical_state(rec_x_R[i][j], g);
                    enforce_physical_state(rec_y_L[i][j], g);
                    enforce_physical_state(rec_y_R[i][j], g);
                }
            }
        }

        // --- ОСНОВНОЙ ШАГ ПО ВРЕМЕНИ ---
        if (equation_type == 8) { // FLIC
            flic_step_2d_local(u_prev, u_next, dt, dx, dy, g,
                               Nx_with_fict_cells, Ny_with_fict_cells, fict_x, fict_y,
                               left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                               cart_comm, i_start, i_end, j_start, j_end,
                               left_rank, right_rank, down_rank, up_rank);
        } else if (equation_type == 9) { // Mader 2DE
            mader_step_2d_local(u_prev, u_next,
                                w_prev, w_next, rho_ref,
                                dt, dx, dy, g,
                                Nx_with_fict_cells, Ny_with_fict_cells, fict_x, fict_y,
                                left_bc_code, right_bc_code, up_bc_code, down_bc_code,
                                solid_mask, mader_config, step,
                                i_start, i_end, j_start, j_end);
        } else {
            u_next = u_prev;
            // --- ОСНОВНОЙ ЦИКЛ ПО ЯЧЕЙКАМ (Расчет потоков и обновление) ---
            for (int i = i_start; i < i_end; ++i) {
                for (int j = j_start; j < j_end; ++j) {
                    std::vector<double> left_flux(M), right_flux(M);
                    std::vector<double> up_flux(M), down_flux(M);

                    // Выбор схемы для расчета потоков
                    if (equation_type == 0) { // Godunov
                        left_flux = godunov_flux_x(u_prev[i - 1][j], u_prev[i][j], g);
                        right_flux = godunov_flux_x(u_prev[i][j], u_prev[i + 1][j], g);
                        up_flux = godunov_flux_y(u_prev[i][j], u_prev[i][j + 1], g);
                        down_flux = godunov_flux_y(u_prev[i][j - 1], u_prev[i][j], g);
                    } else if (equation_type == 1) { // Kolgan
                        // Используем заранее рассчитанные грани
                        left_flux = godunov_flux_x(right_face_x[i - 1][j], left_face_x[i][j], g);
                        right_flux = godunov_flux_x(right_face_x[i][j], left_face_x[i + 1][j], g);
                        up_flux = godunov_flux_y(right_face_y[i][j], left_face_y[i][j + 1], g);
                        down_flux = godunov_flux_y(right_face_y[i][j - 1], left_face_y[i][j], g);
                    } else if (equation_type == 2) { // Rodionov
                        // Используем заранее рассчитанные left/right
                        left_flux = godunov_flux_x(rec_x_R[i - 1][j], rec_x_L[i][j], g);
                        right_flux = godunov_flux_x(rec_x_R[i][j], rec_x_L[i + 1][j], g);
                        up_flux = godunov_flux_y(rec_y_R[i][j], rec_y_L[i][j + 1], g);
                        down_flux = godunov_flux_y(rec_y_R[i][j - 1], rec_y_L[i][j], g);
                    } else if (equation_type == 3) { // HLL
                        left_flux = hll_flux_new(u_prev[i - 1][j], u_prev[i][j], g, 0);
                        right_flux = hll_flux_new(u_prev[i][j], u_prev[i + 1][j], g, 0);
                        up_flux = hll_flux_new(u_prev[i][j], u_prev[i][j + 1], g, 1);
                        down_flux = hll_flux_new(u_prev[i][j - 1], u_prev[i][j], g, 1);
                    } else if (equation_type == 4) { // HLLC
                        left_flux = hllc_flux_new(u_prev[i - 1][j], u_prev[i][j], g, 0);
                        right_flux = hllc_flux_new(u_prev[i][j], u_prev[i + 1][j], g, 0);
                        up_flux = hllc_flux_new(u_prev[i][j], u_prev[i][j + 1], g, 1);
                        down_flux = hllc_flux_new(u_prev[i][j - 1], u_prev[i][j], g, 1);
                    } else if (equation_type == 5) { // Rusanov
                        left_flux = rusanov_2d(u_prev[i - 1][j], u_prev[i][j], g, 0);
                        right_flux = rusanov_2d(u_prev[i][j], u_prev[i + 1][j], g, 0);
                        up_flux = rusanov_2d(u_prev[i][j], u_prev[i][j + 1], g, 1);
                        down_flux = rusanov_2d(u_prev[i][j - 1], u_prev[i][j], g, 1);
                    } else if (equation_type == 6) { // Osher
                        left_flux = osher_flux_2d(u_prev[i - 1][j], u_prev[i][j], g, 0);
                        right_flux = osher_flux_2d(u_prev[i][j], u_prev[i + 1][j], g, 0);
                        up_flux = osher_flux_2d(u_prev[i][j], u_prev[i][j + 1], g, 1);
                        down_flux = osher_flux_2d(u_prev[i][j - 1], u_prev[i][j], g, 1);
                    } else if (equation_type == 7) { // Roe
                        left_flux = roe_flux_2d(u_prev[i - 1][j], u_prev[i][j], g, 0);
                        right_flux = roe_flux_2d(u_prev[i][j], u_prev[i + 1][j], g, 0);
                        up_flux = roe_flux_2d(u_prev[i][j], u_prev[i][j + 1], g, 1);
                        down_flux = roe_flux_2d(u_prev[i][j - 1], u_prev[i][j], g, 1);
                    } else {
                        throw std::runtime_error("Unsupported equation_type: " + std::to_string(equation_type));
                    }

                    // Обновление решения u_next
                    for (int k = 0; k < M; ++k) {
                        // Стандартный конечно-объемный апдейт
                        u_next[i][j][k] = u_prev[i][j][k] - (dt / dx) * (right_flux[k] - left_flux[k])
                                                            - (dt / dy) * (up_flux[k] - down_flux[k]);
                    }

                    enforce_physical_state(u_next[i][j], g);
                }
            }
        }
            /*
        // Установка граничных условий для нового слоя
        for(int i = 0; i < fict_x; ++i) {
            u_next[i] = boundary(u_next[fict_x], left_bc_code, g);
            u_next[Nx_with_fict_cells - 1 - i] = boundary(u_next[Nx_with_fict_cells - fict_x - 1], right_bc_code, g);
        }
            */

        std::swap(u_prev, u_next);
        if (equation_type == 9) {
            std::swap(w_prev, w_next);
        }
        curr_time += dt;
        step++;

        // Вывод и сохранение
        bool should_output = false;
        if (step_num > 0 && step % step_num == 0) {
            should_output = true;
        }
        if (dt_out > 0.0 && curr_time + 1e-12 >= next_output_time) {
            should_output = true;
            while (next_output_time <= curr_time + 1e-12) {
                next_output_time += dt_out;
            }
        }
        if (should_output) {
            gather_to_root(u_prev, Nx, Ny, fict_x, fict_y, i_start, i_end, j_start, j_end, rank, size, cart_comm, dims);
            if (equation_type == 9) {
                gather_scalar_to_root(w_prev, Nx, Ny, fict_x, fict_y, i_start, i_end, j_start, j_end, rank, size, cart_comm, dims);
            }
            if (rank == 0) {
                std::cout << "Step: " << step << ", Time: " << curr_time << ", dt: " << dt << std::endl;
                std::ostringstream basename;
                basename << "step_" << step << "_time_" << std::fixed << std::setprecision(6) << curr_time << ".csv";
                write_state_csv((std::filesystem::path(steps_dir) / basename.str()).string(),
                                u_prev,
                                equation_type == 9 ? &w_prev : nullptr,
                                curr_time,
                                equation_type != 9);
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start);

    // Сохранение финальных результатов
    gather_to_root(u_prev, Nx, Ny, fict_x, fict_y, i_start, i_end, j_start, j_end, rank, size, cart_comm, dims);
    if (equation_type == 9) {
        gather_scalar_to_root(w_prev, Nx, Ny, fict_x, fict_y, i_start, i_end, j_start, j_end, rank, size, cart_comm, dims);
    }
    if (rank == 0) {
        std::filesystem::create_directories(steps_dir);
        write_state_csv((std::filesystem::path(steps_dir) / "final_results.csv").string(),
                        u_prev,
                        equation_type == 9 ? &w_prev : nullptr,
                        curr_time,
                        equation_type != 9);

        // Вывод времени выполнения
        double seconds = duration.count() / 1000.0;
        std::cout << "==========================================" << std::endl;
        std::cout << "Calculation completed! Final time: " << curr_time << std::endl;
        std::cout << "Total execution time: " << std::fixed << std::setprecision(3)
                  << seconds << " seconds (" << duration.count() << " milliseconds)" << std::endl;
        std::cout << "==========================================" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
