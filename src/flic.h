#ifndef FLIC_H
#define FLIC_H

#include "mpi_compat.h"
#include <vector>

// One full 2D FLIC step with Strang splitting:
// X(dt/2) -> Y(dt) -> X(dt/2).
void flic_step_2d(
    const std::vector<std::vector<std::vector<double>>>& u_prev,
    std::vector<std::vector<std::vector<double>>>& u_next,
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
    int down_bc_code);

void flic_step_2d_local(
    const std::vector<std::vector<std::vector<double>>>& u_prev,
    std::vector<std::vector<std::vector<double>>>& u_next,
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
    int up_rank);

#endif
