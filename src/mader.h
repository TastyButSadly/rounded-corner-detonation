#ifndef MADER_H
#define MADER_H

#include <vector>

#include "mpi_compat.h"

using Field3 = std::vector<std::vector<std::vector<double>>>;
using ScalarField = std::vector<std::vector<double>>;
using Mask2D = std::vector<std::vector<unsigned char>>;

struct MaderConfig {
    double visc = 0.5;
    bool slab = true;
    bool step_enabled = false;
    double step_x_end = 0.0;
    double step_y_end = 0.0;
    double corner_radius = 0.0;
    bool reaction_enabled = false;
    bool shargatov_correction = false;
    double gas_constant = 1.0;
    double reaction_rate = 0.0;
    double reaction_activation_energy = 0.0;
    double reaction_heat_release = 0.0;
    double min_temperature = 1000.0;
    double gasw_threshold = 0.02;
    int reaction_delay_steps = 25;
    double min_density_factor = 0.0;
};

ScalarField make_scalar_field(int Nx, int Ny, double value = 0.0);

Mask2D build_step_mask(double x_min, double y_min, double dx, double dy,
                       int Nx, int Ny, int fict_x, int fict_y,
                       const MaderConfig& config);

void apply_step_mask(Field3& u,
                     const Mask2D& solid_mask,
                     const std::vector<double>& fill_cons);

void apply_step_mask_scalar(ScalarField& field,
                            const Mask2D& solid_mask,
                            double fill_value);

double calc_time_step_masked(const Field3& v_cons,
                             const Mask2D& solid_mask,
                             double dx, double dy, double cfl, double g,
                             int fict_x, int fict_y);

double calc_time_step_masked_local(const Field3& v_cons,
                                   const Mask2D& solid_mask,
                                   double dx, double dy, double cfl, double g,
                                   int fict_x, int fict_y,
                                   int i_start, int i_end, int j_start, int j_end);

double mader_temperature(double rho, double pressure, const MaderConfig& config);

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
                           int up_rank);

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
                           int* dims);

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
                   int step_index);

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
                         int j_end);

#endif
