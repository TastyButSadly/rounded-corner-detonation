#ifndef SOLVER
#define SOLVER

#include "mpi_compat.h"
#include <vector>
#include <array>
#include "./grid.h"
#include "grid.h"





std::vector<double> godunov_flux_x(const std::vector<double>& left_param, const std::vector<double>& right_param, double g);

std::vector<double> godunov_flux_y(const std::vector<double>& up_param, const std::vector<double>& down_param, double g);

std::vector<double> cons_to_noncons(const std::vector<double>& v_cons, double g);

std::vector<double> noncons_to_cons(const std::vector<double>& v_noncons, double g);

double calc_time_step(const std::vector<std::vector<std::vector<double>>>& v_cons, 
                      double dx, double dy, double cfl, double g, int fict_cells);

double calc_time_step_local(const std::vector<std::vector<std::vector<double>>>& v_cons,
                            double dx, double dy, double cfl, double g,
                            int i_start, int i_end, int j_start, int j_end);

double calc_sound_speed(const std::vector<double>& v_noncons, double g);

std::vector<double> boundary(const std::vector<double>& v_cons, int boundary_type, double g, int axis);

void enforce_physical_state(std::vector<double>& v_cons, double g);

void set_sod_initial_conditions(std::vector<std::vector<std::vector<double>>>& u_prev, 
                               std::vector<std::vector<std::vector<double>>>& u_next,
                               int Nx, double x_min, double x_max, int Ny, double y_min, 
                               double y_max, double rho_r, double u_r, double v_r, double p_r, 
                                double rho_l, double u_l, double v_l, double p_l, double g, int fict_x, int fict_y); 

std::vector<double> find_solution_x(const std::vector<double>& v_noncons_l, const std::vector<double>& v_noncons_r,
                                  double v_cont, double p_cont, double c_l, double c_r, double S, double g);

std::vector<double> find_solution_y(const std::vector<double>& v_noncons_up, const std::vector<double>& v_noncons_down,
                                  double v_cont, double p_cont, double c_up, double c_down, double S, double g);

std::tuple<double, double> calc_contact_pressure_velocity_x(
    const std::vector<double>& v_ncons_l, 
    const std::vector<double>& v_ncons_r,
    double cl, double cr, double g);

std::tuple<double, double> calc_contact_pressure_velocity_y(
    const std::vector<double>& v_ncons_up, 
    const std::vector<double>& v_ncons_down,
    double cup, double cdown, double g);

std::tuple<double, double> calc_F_and_DF_x(double curr_press, const std::vector<double>& v_ncons, double c, double g);  
std::tuple<double, double> calc_F_and_DF_y(double curr_press, const std::vector<double>& v_ncons, double c, double g);
double pressure_initial_guess_x(const std::vector<double>& v_ncons_l, const std::vector<double>& v_ncons_r, double cl, double cr, double g);
double pressure_initial_guess_y(const std::vector<double>& v_ncons_up, const std::vector<double>& v_ncons_down, double cup, double cdown, double g);
std::vector<double> diff_flux_ncons_x(const std::vector<double>& v_ncons, double g);
std::vector<double> diff_flux_ncons_y(const std::vector<double>& v_ncons, double g);

std::vector<std::vector<std::vector<double>>> compute_analytic_solution_2d(
    double x_min, double x_max, double y_min, double y_max,
    int Nx, int Ny, int fict_x, int fict_y, 
    double t,
    const std::vector<std::vector<std::vector<double>>>& initial_cons,
    int analytic_axis,
    int analytic_profile_index,
    double g);

double minmod(double a, double b);



void Kolgan(const std::vector<double>& v_cons_left, 
           const std::vector<double>& v_cons_center, 
           const std::vector<double>& v_cons_right,
           std::vector<double>& v_cons_left_face, 
           std::vector<double>& v_cons_right_face, 
           double g);




void Kolgan_for_Rodionov(const std::vector<double>& v_cons_left, 
           const std::vector<double>& v_cons_center, 
           const std::vector<double>& v_cons_right,
           std::vector<double>& delta, 
           double g);

// ОБНОВЛЕННАЯ СИГНАТУРА: принимаем delta_x и delta_y
void Rodionov(const std::vector<std::vector<std::vector<double>>>& u_prev,
              std::vector<std::vector<std::vector<double>>>& u_half,
              const std::vector<std::vector<std::vector<double>>>& delta_x,
              const std::vector<std::vector<std::vector<double>>>& delta_y, // Добавлен аргумент
              double dt, double dx, double dy, double g, int Nx, int Ny, int fict_x, int fict_y,
              int left_bc_code, int right_bc_code, int up_bc_code, int down_bc_code);

void Rodionov_local(const std::vector<std::vector<std::vector<double>>>& u_prev,
                    std::vector<std::vector<std::vector<double>>>& u_half,
                    const std::vector<std::vector<std::vector<double>>>& delta_x,
                    const std::vector<std::vector<std::vector<double>>>& delta_y,
                    double dt, double dx, double dy, double g,
                    int i_start, int i_end, int j_start, int j_end);


std::vector<double> hll_flux(const std::vector<double>& left_cons, 
                            const std::vector<double>& right_cons, 
                            double g);


std::vector<double> hllc_flux(const std::vector<double>& left_cons, 
                             const std::vector<double>& right_cons, 
                             double g);


std::vector<double> Rusanov(const std::vector<double>& v_right, const std::vector<double>& v_left, 
                            double dt, double dx, double g);

std::vector<double> osher_flux(const std::vector<double>& left_cons, 
                               const std::vector<double>& right_cons, 
                               double g);

std::vector<double> roe_flux(const std::vector<double>& left_cons, 
                            const std::vector<double>& right_cons, double g);
std::vector<double> hllc_flux_new(const std::vector<double>& left_cons, 
                              const std::vector<double>& right_cons, 
                              double g, int axis);

std::vector<double> rusanov_flux(const std::vector<double>& left_cons, 
                                 const std::vector<double>& right_cons, 
                                 double g);

std::vector<double> hll_flux_new(const std::vector<double>& left_cons, 
                             const std::vector<double>& right_cons, 
                             double g, int axis);


                             std::vector<double> roe_flux_new(const std::vector<double>& left_cons, 
                             const std::vector<double>& right_cons, 
                             double g);

std::vector<double> osher_flux(const std::vector<double>& left_cons, 
                               const std::vector<double>& right_cons, 
                               double g);

// Осher интермедиаты в примитивах (используется внутри osher_flux)
struct OsherStates {
    std::array<double,3> U_1_3;
    std::array<double,3> U_2_3;
    std::array<double,3> U_s0;
    std::array<double,3> U_s1;
    double u_star;
    double p_star;
    double c_1_3;
    double c_2_3;
};

OsherStates calc_osher_states(const double prim_l[3], 
                             const double prim_r[3], 
                             double g);


std::vector<double> osher_2d(const std::vector<double>& UL, const std::vector<double>& UR, double g, int dir);
std::vector<double> roe_2d(const std::vector<double>& UL, const std::vector<double>& UR, double g, int dir);
std::vector<double> rusanov_2d(const std::vector<double>& UL, const std::vector<double>& UR, double g, int dir);
std::vector<double> osher_flux_2d(const std::vector<double>& left_cons, 
                                  const std::vector<double>& right_cons, 
                                  double g, int axis);
                                  // --- ОПТИМИЗИРОВАННЫЙ ROE 2D ---
std::vector<double> roe_flux_2d(const std::vector<double>& left_cons, 
                               const std::vector<double>& right_cons, 
                               double g, int axis);

void get_subdomain_bounds(int N, int parts, int rank, int& start, int& end);

void exchange_halos_global_width(std::vector<std::vector<std::vector<double>>>& u,
                                 MPI_Comm cart_comm, int halo_x, int halo_y,
                                 int i_start, int i_end, int j_start, int j_end,
                                 int left_rank, int right_rank, int down_rank, int up_rank);

void exchange_halos_global(std::vector<std::vector<std::vector<double>>>& u,
                           MPI_Comm cart_comm, int fict_x, int fict_y,
                           int i_start, int i_end, int j_start, int j_end,
                           int left_rank, int right_rank, int down_rank, int up_rank);

void apply_physical_boundaries_local(std::vector<std::vector<std::vector<double>>>& u,
                                     int Nx, int Ny, int fict_x, int fict_y,
                                     int i_start, int i_end, int j_start, int j_end,
                                     int left_rank, int right_rank, int down_rank, int up_rank,
                                     int left_bc_code, int right_bc_code, int up_bc_code, int down_bc_code,
                                     double g);

void gather_to_root(std::vector<std::vector<std::vector<double>>>& u,
                    int Nx, int Ny, int fict_x, int fict_y,
                    int i_start, int i_end, int j_start, int j_end,
                    int rank, int size, MPI_Comm cart_comm, int* dims);
#endif
