#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#pragma optimize("gt", on)           // Агрессивная оптимизация скорости
#endif



#include "solver.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <limits>

//RimanSolver1D::RimanSolver1D(Grid* grid) {

//}



double calc_sound_speed(const std::vector<double>& v_noncons, double g){
    return std::sqrt(g * v_noncons[P] / v_noncons[RHO]);
}

// Обновленный расчет шага по времени для 2D (CFL condition)
double calc_time_step(const std::vector<std::vector<std::vector<double>>>& v_cons, 
                      double dx, double dy, double cfl, double g, int fict_cells) {
    double max_sig_x = 0.0;
    double max_sig_y = 0.0;
    
    int Nx = v_cons.size();
    int Ny = v_cons[0].size();

    // Проход по физической области
    for(int i = fict_cells; i < Nx - fict_cells; ++i) {
        for(int j = fict_cells; j < Ny - fict_cells; ++j) {
            auto prim = cons_to_noncons(v_cons[i][j], g);
            double c = calc_sound_speed(prim, g);
            
            double sig_x = std::abs(prim[U]) + c;
            double sig_y = std::abs(prim[V]) + c;
            
            if (sig_x > max_sig_x) max_sig_x = sig_x;
            if (sig_y > max_sig_y) max_sig_y = sig_y;
        }
    }
    
    
    
    // dt = CFL / (max_lambda_x/dx + max_lambda_y/dy)
    double inv_dt = (max_sig_x / dx) + (max_sig_y / dy);
    return cfl / inv_dt;
}

double calc_time_step_local(const std::vector<std::vector<std::vector<double>>>& v_cons,
                            double dx, double dy, double cfl, double g,
                            int i_start, int i_end, int j_start, int j_end) {
    double local_dt = std::numeric_limits<double>::max();

    for (int i = i_start; i < i_end; ++i) {
        for (int j = j_start; j < j_end; ++j) {
            const std::vector<double> prim = cons_to_noncons(v_cons[i][j], g);
            const double c = calc_sound_speed(prim, g);
            const double sig_x = std::abs(prim[U]) + c;
            const double sig_y = std::abs(prim[V]) + c;
            const double inv_dt = (sig_x / dx) + (sig_y / dy);
            if (inv_dt > 0.0) {
                local_dt = std::min(local_dt, cfl / inv_dt);
            }
        }
    }

    if (local_dt == std::numeric_limits<double>::max()) {
        return cfl * std::min(dx, dy);
    }
    return local_dt;
}


namespace {
constexpr double RHO_FLOOR = 1e-8;
constexpr double P_FLOOR = 1e-8;
constexpr double EPS = 1e-6;
constexpr int MAX_ITER_NUM = 20;
constexpr double P_MAX_RATIO = 2.0;

struct HaloExchangeScratch {
    std::vector<double> send_left;
    std::vector<double> recv_left;
    std::vector<double> send_right;
    std::vector<double> recv_right;
    std::vector<double> send_down;
    std::vector<double> recv_down;
    std::vector<double> send_up;
    std::vector<double> recv_up;
};

HaloExchangeScratch& halo_exchange_scratch() {
    static thread_local HaloExchangeScratch scratch;
    return scratch;
}
}

std::vector<double> cons_to_noncons(const std::vector<double>& v_cons, double g) {
    double rho = v_cons[r];
    if (rho <= RHO_FLOOR) rho = RHO_FLOOR;

    double u = v_cons[ru] / rho;
    double v = v_cons[rv] / rho; // 2D component
    
    // E = p/(g-1) + 0.5*rho*(u^2 + v^2)
    // p = (g-1) * (E - 0.5*rho*(u^2 + v^2))
    double kin_energy = 0.5 * rho * (u * u + v * v);
    double p = (g - 1.0) * (v_cons[e] - kin_energy);

    if (p < P_FLOOR) p = P_FLOOR;
    
    return std::vector<double>{rho, u, v, p};
}

std::vector<double> noncons_to_cons(const std::vector<double>& v_noncons, double g) {
    double rho = v_noncons[RHO];
    double u = v_noncons[U];
    double v = v_noncons[V];
    double p = v_noncons[P];

    double E = p / (g - 1.0) + 0.5 * rho * (u * u + v * v);
    
    return std::vector<double>{rho, rho * u, rho * v, E};
}

void enforce_physical_state(std::vector<double>& v_cons, double g) {
    double rho = v_cons[r];
    if (rho < RHO_FLOOR) rho = RHO_FLOOR;

    double u = v_cons[ru] / rho;
    double v = v_cons[rv] / rho;
    
    double kinetic = 0.5 * rho * (u * u + v * v);
    double internal = v_cons[e] - kinetic;
    double p = (g - 1.0) * internal;

    if (p < P_FLOOR) {
        p = P_FLOOR;
        internal = p / (g - 1.0);
        v_cons[e] = internal + kinetic;
    }
    
    v_cons[r] = rho;
    v_cons[ru] = rho * u;
    v_cons[rv] = rho * v;
    // e уже обновлено
}



std::vector<double> godunov_flux_x(const std::vector<double>& left_param, const std::vector<double>& right_param, double g){
    std::vector<double> v_noncons_l = cons_to_noncons(left_param, g);
    std::vector<double> v_noncons_r = cons_to_noncons(right_param, g);

    double r_l = v_noncons_l[RHO];
    double u_l = v_noncons_l[U];
    double p_l = v_noncons_l[P];

    double r_r = v_noncons_r[RHO];
    double u_r = v_noncons_r[U];
    double p_r = v_noncons_r[P];

    double c_l = calc_sound_speed(v_noncons_l, g);
    double c_r = calc_sound_speed(v_noncons_r, g);

    //double p_cont, v_cont;

    if ( 2. * (c_l + c_r) / (g - 1) <= u_r - u_l){
        std::cout << "error vacuum in godunov flux" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto [v_cont, p_cont] = calc_contact_pressure_velocity_x(v_noncons_l, v_noncons_r, c_l, c_r, g);

    
    std::vector<double> v_ncons = find_solution_x(v_noncons_l, v_noncons_r, v_cont, p_cont, c_l, c_r, 0.0, g);


    std::vector<double> flux = diff_flux_ncons_x(v_ncons, g);

    return flux;
    
}

std::vector<double> godunov_flux_y(const std::vector<double>& down_param, const std::vector<double>& up_param, double g){
    std::vector<double> v_noncons_up = cons_to_noncons(up_param, g);
    std::vector<double> v_noncons_down = cons_to_noncons(down_param, g);

    double r_up = v_noncons_up[RHO];
    double u_up = v_noncons_up[V];
    double p_up = v_noncons_up[P];

    double r_down = v_noncons_down[RHO];
    double u_down = v_noncons_down[V];
    double p_down = v_noncons_down[P];

    double c_up = calc_sound_speed(v_noncons_up, g);
    double c_down = calc_sound_speed(v_noncons_down, g);

    //double p_cont, v_cont;

    if ( 2. * (c_up + c_down) / (g - 1) <= u_up - u_down){
        std::cout << "error vacuum in godunov flux" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto [v_cont, p_cont] = calc_contact_pressure_velocity_y(v_noncons_down, v_noncons_up, c_down, c_up, g);

    
    std::vector<double> v_ncons = find_solution_y(v_noncons_down, v_noncons_up, v_cont, p_cont, c_down, c_up, 0.0, g);


    std::vector<double> flux = diff_flux_ncons_y(v_ncons, g);

    return flux;
    
}









std::vector<double> diff_flux_ncons_x(const std::vector<double>& v_ncons, double g) {
    double r_u = v_ncons[RHO] * v_ncons[U];
    double p_ru_2 = v_ncons[P] + (v_ncons[RHO] * v_ncons[U] * v_ncons[U]);
    double r_u_v = v_ncons[RHO] * v_ncons[U] * v_ncons[V];
    double u_e_p = v_ncons[U] * (v_ncons[P] + 0.5 * v_ncons[RHO] * (v_ncons[U] * v_ncons[U] + v_ncons[V] * v_ncons[V]) +( v_ncons[P] / (g - 1.0)));
    std::vector<double> flux{r_u, p_ru_2, r_u_v,u_e_p};

    return flux;
}

std::vector<double> diff_flux_ncons_y(const std::vector<double>& v_ncons, double g) {
    double r_v = v_ncons[RHO] * v_ncons[V];
    double r_v_u = v_ncons[RHO] * v_ncons[V] * v_ncons[U];
    double p_rv_2 = v_ncons[P] + (v_ncons[RHO] * fabs(v_ncons[V]) * fabs(v_ncons[V]));
    double v_e_p = v_ncons[V] * (v_ncons[P] + 0.5 * v_ncons[RHO] * (v_ncons[V] * v_ncons[V] + v_ncons[U] * v_ncons[U]) + v_ncons[P] / (g - 1.0));
    std::vector<double> flux{r_v, r_v_u, p_rv_2,v_e_p};

    return flux;
}
 
std::vector<double> boundary(const std::vector<double>& v_cons, int boundary_type, double g, int axis) {
    std::vector<double> bound_ncons = cons_to_noncons(v_cons, g);
    
    // axis == 0: граница по X (отражение u)
    // axis == 1: граница по Y (отражение v)

    if (boundary_type == 1) { // Стенка (Slip wall)
        if (axis == 0) bound_ncons[U] = -bound_ncons[U]; // Отражение нормальной скорости X
        if (axis == 1) bound_ncons[V] = -bound_ncons[V]; // Отражение нормальной скорости Y
    }
    else if (boundary_type == 2) { // Свободный выход (extrapolation)
        // Ничего не меняем, градиент 0
    }
    else {
        // Default outflow
    }
    
    std::vector<double> ret = noncons_to_cons(bound_ncons, g);
    enforce_physical_state(ret, g);
    return ret;
}



void set_sod_initial_conditions(std::vector<std::vector<std::vector<double>>>& u_prev, 
                               std::vector<std::vector<std::vector<double>>>& u_next,
                               int Nx, double x_min, double x_max, int Ny, double y_min, 
                               double y_max, double rho_r, double u_r, double v_r, double p_r, 
                                double rho_l, double u_l, double v_l, double p_l, double g, int fict_x, int fict_y) {
    

    
    
    
    
    int Nx_physical = Nx + 2 * fict_x;
    int Ny_physical = Ny + 2 * fict_y;
    
    double dx = (x_max - x_min) / Nx;
    double dy = (y_max - y_min) / Ny;
    
    //u_prev.resize(Nx, std::vector<double>(M, 0.0));
    //u_next.resize(Nx, std::vector<double>(M, 0.0));
   
    
    double x_mid = (x_min + x_max) / 2.0;
    double y_mid = (y_min + y_max) / 2.0;
    
    for (int i = fict_x; i < Nx_physical - fict_x; i++) {
        for (int j = fict_y; j < Ny_physical - fict_y; j++){
            double x_cell = x_min + (i - fict_x + 0.5) * dx;
        
            if (x_cell < x_mid) {
            
                u_prev[i][j][r] = rho_l;  
                u_prev[i][j][ru] = rho_l * u_l;
                u_prev[i][j][rv] = rho_l * v_l;  
                u_prev[i][j][e] = 0.5 * rho_l * (u_l * u_l + v_l * v_l) + p_l / (g - 1.0); 
            } else {
            
                u_prev[i][j][r] = rho_r; 
                u_prev[i][j][ru] = rho_r * u_r;
                u_prev[i][j][rv] = rho_r * v_r;     
                u_prev[i][j][e] = 0.5 * rho_r * (u_r * u_r + v_r * v_r) + p_r / (g - 1.0);   
            }
        
    
        u_next[i][j] = u_prev[i][j];
       
        }      
    }
}

std::vector<double> find_solution_x(const std::vector<double>& v_noncons_l, const std::vector<double>& v_noncons_r,
                                  double v_cont, double p_cont, double c_l, double c_r, double S, double g){
    const double rl = v_noncons_l[RHO];
    const double vl = v_noncons_l[U];
    const double ul = v_noncons_l[V];
    const double pl = v_noncons_l[P];

    const double rr = v_noncons_r[RHO];
    const double vr = v_noncons_r[U];
    const double ur = v_noncons_r[V];
    const double pr = v_noncons_r[P];

    const double g1 = 0.5 * (g - 1.0) / g;
    const double g2 = 0.5 * (g + 1.0) / g;
    const double g3 = 2.0 * g / (g - 1.0);
    const double g4 = 2.0 / (g - 1.0);
    const double g5 = 2.0 / (g + 1.0);
    const double g6 = (g - 1.0) / (g + 1.0);
    const double g7 = 0.5 * (g - 1.0);

    double rho_star = 0.0;
    double vel_star = 0.0;
    double press_star = 0.0;

    if (S <= v_cont) {
        if (p_cont <= pl) {
            const double shl = vl - c_l;
            if (S <= shl) {
                rho_star = rl;
                vel_star = vl;
                press_star = pl;
            } else {
                const double p_ratio = p_cont / pl;
                const double cml = c_l * std::pow(p_ratio, g1);
                const double stl = v_cont - cml;
                if (S > stl) {
                    rho_star = rl * std::pow(p_ratio, 1.0 / g);
                    vel_star = v_cont;
                    press_star = p_cont;
                } else {
                    const double vel = g5 * (c_l + g7 * vl + S);
                    const double c = g5 * (c_l + g7 * (vl - S));
                    rho_star = rl * std::pow(c / c_l, g4);
                    press_star = pl * std::pow(c / c_l, g3);
                    vel_star = vel;
                }
            }
        } else {
            const double p_ratio = p_cont / pl;
            const double sl = vl - c_l * std::sqrt(g2 * p_ratio + g1);
            if (S <= sl) {
                rho_star = rl;
                vel_star = vl;
                press_star = pl;
            } else {
                rho_star = rl * (p_ratio + g6) / (p_ratio * g6 + 1.0);
                vel_star = v_cont;
                press_star = p_cont;
            }
        }
    } else {
        if (p_cont > pr) {
            const double p_ratio = p_cont / pr;
            const double sr = vr + c_r * std::sqrt(g2 * p_ratio + g1);
            if (S >= sr) {
                rho_star = rr;
                vel_star = vr;
                press_star = pr;
            } else {
                rho_star = rr * (p_ratio + g6) / (p_ratio * g6 + 1.0);
                vel_star = v_cont;
                press_star = p_cont;
            }
        } else {
            const double shr = vr + c_r;
            if (S >= shr) {
                rho_star = rr;
                vel_star = vr;
                press_star = pr;
            } else {
                const double p_ratio = p_cont / pr;
                const double cmr = c_r * std::pow(p_ratio, g1);
                const double str = v_cont + cmr;
                if (S <= str) {
                    rho_star = rr * std::pow(p_ratio, 1.0 / g);
                    vel_star = v_cont;
                    press_star = p_cont;
                } else {
                    const double vel = g5 * (-c_r + g7 * vr + S);
                    const double c = g5 * (c_r - g7 * (vr - S));
                    rho_star = rr * std::pow(c / c_r, g4);
                    press_star = pr * std::pow(c / c_r, g3);
                    vel_star = vel;
                }
            }
        }
    }

    std::vector<double> v_noncons(M);
    v_noncons[RHO] = rho_star;
    v_noncons[U] = vel_star;
    v_noncons[V] = (S <= v_cont) ? ul : ur;
    v_noncons[P] = press_star;
    return v_noncons;
}




std::vector<double> find_solution_y(const std::vector<double>& v_noncons_l, const std::vector<double>& v_noncons_r,
                                  double v_cont, double p_cont, double c_l, double c_r, double S, double g){
    const double rl = v_noncons_l[RHO];
    const double vl = v_noncons_l[U];
    const double ul = v_noncons_l[V];
    const double pl = v_noncons_l[P];

    const double rr = v_noncons_r[RHO];
    const double vr = v_noncons_r[U];
    const double ur = v_noncons_r[V];
    const double pr = v_noncons_r[P];

    const double g1 = 0.5 * (g - 1.0) / g;
    const double g2 = 0.5 * (g + 1.0) / g;
    const double g3 = 2.0 * g / (g - 1.0);
    const double g4 = 2.0 / (g - 1.0);
    const double g5 = 2.0 / (g + 1.0);
    const double g6 = (g - 1.0) / (g + 1.0);
    const double g7 = 0.5 * (g - 1.0);

    double rho_star = 0.0;
    double vel_star = 0.0;
    double press_star = 0.0;

    if (S <= v_cont) {
        if (p_cont <= pl) {
            const double shl = ul - c_l;
            if (S <= shl) {
                rho_star = rl;
                vel_star = ul;
                press_star = pl;
            } else {
                const double p_ratio = p_cont / pl;
                const double cml = c_l * std::pow(p_ratio, g1);
                const double stl = v_cont - cml;
                if (S > stl) {
                    rho_star = rl * std::pow(p_ratio, 1.0 / g);
                    vel_star = v_cont;
                    press_star = p_cont;
                } else {
                    const double vel = g5 * (c_l + g7 * ul + S);
                    const double c = g5 * (c_l + g7 * (ul - S));
                    rho_star = rl * std::pow(c / c_l, g4);
                    press_star = pl * std::pow(c / c_l, g3);
                    vel_star = vel;
                }
            }
        } else {
            const double p_ratio = p_cont / pl;
            const double sl = ul - c_l * std::sqrt(g2 * p_ratio + g1);
            if (S <= sl) {
                rho_star = rl;
                vel_star = ul;
                press_star = pl;
            } else {
                rho_star = rl * (p_ratio + g6) / (p_ratio * g6 + 1.0);
                vel_star = v_cont;
                press_star = p_cont;
            }
        }
    } else {
        if (p_cont > pr) {
            const double p_ratio = p_cont / pr;
            const double sr = ur + c_r * std::sqrt(g2 * p_ratio + g1);
            if (S >= sr) {
                rho_star = rr;
                vel_star = ur;
                press_star = pr;
            } else {
                rho_star = rr * (p_ratio + g6) / (p_ratio * g6 + 1.0);
                vel_star = v_cont;
                press_star = p_cont;
            }
        } else {
            const double shr = ur + c_r;
            if (S >= shr) {
                rho_star = rr;
                vel_star = ur;
                press_star = pr;
            } else {
                const double p_ratio = p_cont / pr;
                const double cmr = c_r * std::pow(p_ratio, g1);
                const double str = v_cont + cmr;
                if (S <= str) {
                    rho_star = rr * std::pow(p_ratio, 1.0 / g);
                    vel_star = v_cont;
                    press_star = p_cont;
                } else {
                    const double vel = g5 * (-c_r + g7 * ur + S);
                    const double c = g5 * (c_r - g7 * (ur - S));
                    rho_star = rr * std::pow(c / c_r, g4);
                    press_star = pr * std::pow(c / c_r, g3);
                    vel_star = vel;
                }
            }
        }
    }

    std::vector<double> v_noncons(M);
    v_noncons[RHO] = rho_star;
    v_noncons[U] = (S <= v_cont) ? vl : vr;
    v_noncons[V] = vel_star;
    v_noncons[P] = press_star;
    return v_noncons;
}






double pressure_initial_guess_x(const std::vector<double>& v_ncons_l, const std::vector<double>& v_ncons_r, double cl, double cr, double g) {

    double rl, vl, pl;                  
    double rr, vr, pr;                  
                        
  
      
    double p_lin;
    double p_min, p_max;                
    double p_ratio;                    
    double p1, p2, g1, g2;              
    
    rl = v_ncons_l[RHO];
    vl = v_ncons_l[U];
    pl = v_ncons_l[P];
    rr = v_ncons_r[RHO];
    vr = v_ncons_r[U];
    pr = v_ncons_r[P];
    


    p_lin = std::max( 0.0, 0.5 * ( pl + pr ) - 0.125 * ( vr - vl ) * ( rl + rr ) * ( cl + cr ) );
    p_min = std::min( pl, pr );
    p_max = std::max( pl, pr );
    p_ratio = p_max / p_min;

    if ( ( p_ratio <= P_MAX_RATIO ) &&
       ( ( p_min < p_lin && p_lin < p_max ) || ( std::fabs( p_min - p_lin ) < EPS || std::fabs( p_max - p_lin ) < EPS ) ) ) {
        return p_lin;
    }

    if ( p_lin < p_min ) {
        g1 = 0.5 * ( g - 1.0 ) / g;
        return std::pow( ( ( cl + cr - 0.5 * ( g - 1.0 ) * ( vr - vl ) ) / ( cl / std::pow( pl, g1 ) + cr / std::pow( pr, g1 ) ) ), 1.0 / g1 );
    }

    g1 = 2.0 / ( g + 1.0 );
    g2 = ( g - 1.0 ) / ( g + 1.0 );
    p1 = std::sqrt( g1 / rl / ( g2 * pl + p_lin ) );
    p2 = std::sqrt( g1 / rr / ( g2 * pr + p_lin ) );
    return ( p1 * pl + p2 * pr - ( vr - vl ) ) / ( p1 + p2 );
    
}



double pressure_initial_guess_y(const std::vector<double>& v_ncons_l, const std::vector<double>& v_ncons_r, double cl, double cr, double g) {

    double rl, vl, pl;                  
    double rr, vr, pr;                  
                        
  
      
    double p_lin;
    double p_min, p_max;                
    double p_ratio;                    
    double p1, p2, g1, g2;              
    
    rl = v_ncons_l[RHO];
    vl = v_ncons_l[V];
    pl = v_ncons_l[P];
    rr = v_ncons_r[RHO];
    vr = v_ncons_r[V];
    pr = v_ncons_r[P];
    


    p_lin = std::max( 0.0, 0.5 * ( pl + pr ) - 0.125 * ( vr - vl ) * ( rl + rr ) * ( cl + cr ) );
    p_min = std::min( pl, pr );
    p_max = std::max( pl, pr );
    p_ratio = p_max / p_min;

    if ( ( p_ratio <= P_MAX_RATIO ) &&
       ( ( p_min < p_lin && p_lin < p_max ) || ( std::fabs( p_min - p_lin ) < EPS || std::fabs( p_max - p_lin ) < EPS ) ) ) {
        return p_lin;
    }

    if ( p_lin < p_min ) {
        g1 = 0.5 * ( g - 1.0 ) / g;
        return std::pow( ( ( cl + cr - 0.5 * ( g - 1.0 ) * ( vr - vl ) ) / ( cl / std::pow( pl, g1 ) + cr / std::pow( pr, g1 ) ) ), 1.0 / g1 );
    }

    g1 = 2.0 / ( g + 1.0 );
    g2 = ( g - 1.0 ) / ( g + 1.0 );
    p1 = std::sqrt( g1 / rl / ( g2 * pl + p_lin ) );
    p2 = std::sqrt( g1 / rr / ( g2 * pr + p_lin ) );
    return ( p1 * pl + p2 * pr - ( vr - vl ) ) / ( p1 + p2 );
    
}





std::tuple<double, double> calc_F_and_DF_x(double curr_press, const std::vector<double>& v_ncons, double c, double g) {
    double rho, v, p;
    double p_ratio, fg, q;
    double F, DF;

    rho = v_ncons[RHO];
    v = v_ncons[U];
    p = v_ncons[P];
    
    p_ratio = curr_press / p;
    if (curr_press <= p) {
        
        fg = 2.0 / (g - 1.0);
        F = fg * c * (std::pow(p_ratio, 1.0 / fg / g) - 1.0);
        DF = (1.0 / (rho * c)) * std::pow(p_ratio, -0.5 * (g + 1.0) / g);
    } else {
        
        q = std::sqrt(0.5 * (g + 1.0) / g * p_ratio + 0.5 * (g - 1.0) / g);
        F = (curr_press - p) / (c * rho * q);
        DF = 0.25 * ((g + 1.0) * p_ratio + 3 * g - 1.0) / (g * rho * c * std::pow(q, 3.0));
    }

    return std::make_tuple(F, DF);
}


std::tuple<double, double> calc_F_and_DF_y(double curr_press, const std::vector<double>& v_ncons, double c, double g) {
    double rho, v, p;
    double p_ratio, fg, q;
    double F, DF;

    rho = v_ncons[RHO];
    v = v_ncons[V];
    p = v_ncons[P];
    
    p_ratio = curr_press / p;
    if (curr_press <= p) {
        
        fg = 2.0 / (g - 1.0);
        F = fg * c * (std::pow(p_ratio, 1.0 / fg / g) - 1.0);
        DF = (1.0 / (rho * c)) * std::pow(p_ratio, -0.5 * (g + 1.0) / g);
    } else {
        
        q = std::sqrt(0.5 * (g + 1.0) / g * p_ratio + 0.5 * (g - 1.0) / g);
        F = (curr_press - p) / (c * rho * q);
        DF = 0.25 * ((g + 1.0) * p_ratio + 3 * g - 1.0) / (g * rho * c * std::pow(q, 3.0));
    }

    return std::make_tuple(F, DF);
}




std::tuple<double, double> calc_contact_pressure_velocity_x(
    const std::vector<double>& v_ncons_l, 
    const std::vector<double>& v_ncons_r,
    double cl, double cr, double g) 
{
    double vl, vr;
    double p_old;
    double fl, fr;
    double fld, frd;
    int iter_num = 0;
    double criteria;
    double p_cont, v_cont;

    vl = v_ncons_l[U];
    vr = v_ncons_r[U];
    
    if (2.0 * (cl + cr) / (g - 1.0) <= vr - vl) {
        std::cout << "\ncalc_contact_pressure_velocity_x -> vacuum is generated\n" << std::endl;
        exit(EXIT_FAILURE);
    }

    
    p_old = pressure_initial_guess_x(v_ncons_l, v_ncons_r, cl, cr, g);
    if (p_old < 0.0) {
        std::cout << "\ncalc_contact_pressure_velocity_x -> initial pressure guess is negative " << std::endl;
        exit(EXIT_FAILURE);
    }
    
    criteria = EPS + 1.0;
    while (criteria > EPS) {
        
        auto [fl_temp, fld_temp] = calc_F_and_DF_x(p_old, v_ncons_l, cl, g);
        auto [fr_temp, frd_temp] = calc_F_and_DF_x(p_old, v_ncons_r, cr, g);
        
        fl = fl_temp;
        fld = fld_temp;
        fr = fr_temp;
        frd = frd_temp;
        
        p_cont = p_old - (fl + fr + vr - vl) / (fld + frd);
        criteria = 2.0 * std::fabs((p_cont - p_old) / (p_cont + p_old));
        iter_num++;
        
        if (iter_num > MAX_ITER_NUM) {
            std::cout << ("\ncalc_contact_pressure_velocity_x -> number of iterations exceeds the maximum value.\n")<<std::endl;
            exit(EXIT_FAILURE);
        }
        if (p_cont < 0.0) {
            std::cout<<("\ncalc_contact_pressure_velocity_x -> pressure is negative.\n")<<std::endl;
            exit(EXIT_FAILURE);            
        }
        p_old = p_cont;
    }

    
    v_cont = 0.5 * (vl + vr + fr - fl);

    return std::make_tuple(v_cont, p_cont);
}





std::tuple<double, double> calc_contact_pressure_velocity_y(
    const std::vector<double>& v_ncons_l, 
    const std::vector<double>& v_ncons_r,
    double cl, double cr, double g) 
{
    double vl, vr;
    double p_old;
    double fl, fr;
    double fld, frd;
    int iter_num = 0;
    double criteria;
    double p_cont, v_cont;

    vl = v_ncons_l[V];
    vr = v_ncons_r[V];
    
    if (2.0 * (cl + cr) / (g - 1.0) <= vr - vl) {
        std::cout << "\ncalc_contact_pressure_velocity_y -> vacuum is generated\n" << std::endl;
        exit(EXIT_FAILURE);
    }

    
    p_old = pressure_initial_guess_y(v_ncons_l, v_ncons_r, cl, cr, g);
    if (p_old < 0.0) {
        std::cout << "\ncalc_contact_pressure_velocity_y -> initial pressure guess is negative " << std::endl;
        exit(EXIT_FAILURE);
    }
    
    criteria = EPS + 1.0;
    while (criteria > EPS) {
        
        auto [fl_temp, fld_temp] = calc_F_and_DF_y(p_old, v_ncons_l, cl, g);
        auto [fr_temp, frd_temp] = calc_F_and_DF_y(p_old, v_ncons_r, cr, g);
        
        fl = fl_temp;
        fld = fld_temp;
        fr = fr_temp;
        frd = frd_temp;
        
        p_cont = p_old - (fl + fr + vr - vl) / (fld + frd);
        criteria = 2.0 * std::fabs((p_cont - p_old) / (p_cont + p_old));
        iter_num++;
        
        if (iter_num > MAX_ITER_NUM) {
            std::cout << ("\ncalc_contact_pressure_velocity_y -> number of iterations exceeds the maximum value.\n")<<std::endl;
            exit(EXIT_FAILURE);
        }
        if (p_cont < 0.0) {
            std::cout<<("\ncalc_contact_pressure_velocity_y -> pressure is negative.\n")<<std::endl;
            std::cout << "\n p_cont = " << p_cont << "\n" << std::endl;
            exit(EXIT_FAILURE);            
        }
        p_old = p_cont;
    }

    
    v_cont = 0.5 * (vl + vr + fr - fl);

    return std::make_tuple(v_cont, p_cont);
}









std::vector<std::vector<std::vector<double>>> compute_analytic_solution_2d(
    double x_min, double x_max, double y_min, double y_max,
    int Nx, int Ny, int fict_x, int fict_y, 
    double t,
    const std::vector<std::vector<std::vector<double>>>& initial_cons,
    int analytic_axis,
    int analytic_profile_index,
    double g) 
{
    // Расчет параметров сетки
    double dx = (x_max - x_min) / Nx;
    double dy = (y_max - y_min) / Ny;
    int Nx_total = Nx + 2 * fict_x;
    int Ny_total = Ny + 2 * fict_y;
    
    // Результирующий массив: [Nx][Ny][4 величины: rho, u, v, p]
    std::vector<std::vector<std::vector<double>>> analytic_solution(
        Nx_total, std::vector<std::vector<double>>(Ny_total, std::vector<double>(M))
    );
    
    // Константы для удобства (предполагаем стандартный порядок)
    

    auto clamp_index = [](int idx, int lo, int hi) {
        if (idx < lo) return lo;
        if (idx > hi) return hi;
        return idx;
    };

    // Ищем интерфейс разрыва на выбранном профиле по максимальному скачку в примитивах.
    auto jump_score = [g](const std::vector<double>& a_cons, const std::vector<double>& b_cons) {
        std::vector<double> a = cons_to_noncons(a_cons, g);
        std::vector<double> b = cons_to_noncons(b_cons, g);
        return std::abs(a[RHO] - b[RHO]) + std::abs(a[U] - b[U]) + std::abs(a[V] - b[V]) + std::abs(a[P] - b[P]);
    };

    const bool solve_along_x = (analytic_axis == 0);
    if (solve_along_x) {
        const int j_profile_phys = (analytic_profile_index >= 0) ? analytic_profile_index : Ny / 4;
        const int j = fict_y + clamp_index(j_profile_phys, 0, Ny - 1);

        int i_left = fict_x + (Nx / 2) - 1;
        i_left = clamp_index(i_left, fict_x, Nx_total - fict_x - 2);
        double best = -std::numeric_limits<double>::infinity();
        for (int i = fict_x; i < Nx_total - fict_x - 1; ++i) {
            double s = jump_score(initial_cons[i][j], initial_cons[i + 1][j]);
            if (s > best) {
                best = s;
                i_left = i;
            }
        }
        const int i_right = i_left + 1;

        std::vector<double> left_prim = cons_to_noncons(initial_cons[i_left][j], g);
        std::vector<double> right_prim = cons_to_noncons(initial_cons[i_right][j], g);
        double cl = calc_sound_speed(left_prim, g);
        double cr = calc_sound_speed(right_prim, g);
        auto [v_cont, p_cont] = calc_contact_pressure_velocity_x(left_prim, right_prim, cl, cr, g);

        const double x_discont = x_min + (i_left - fict_x + 1) * dx;
        for (int i = fict_x; i < Nx_total - fict_x; ++i) {
            double x = x_min + (i - fict_x + 0.5) * dx;
            double xi = (t > 1e-12) ? (x - x_discont) / t : (x < x_discont ? -1e10 : 1e10);
            std::vector<double> res_1d = find_solution_x(left_prim, right_prim, v_cont, p_cont, cl, cr, xi, g);
            for (int jj = fict_y; jj < Ny_total - fict_y; ++jj) {
                analytic_solution[i][jj] = res_1d;
            }
        }
    } else {
        const int i_profile_phys = (analytic_profile_index >= 0) ? analytic_profile_index : Nx / 4;
        const int i = fict_x + clamp_index(i_profile_phys, 0, Nx - 1);

        int j_down = fict_y + (Ny / 2) - 1;
        j_down = clamp_index(j_down, fict_y, Ny_total - fict_y - 2);
        double best = -std::numeric_limits<double>::infinity();
        for (int j = fict_y; j < Ny_total - fict_y - 1; ++j) {
            double s = jump_score(initial_cons[i][j], initial_cons[i][j + 1]);
            if (s > best) {
                best = s;
                j_down = j;
            }
        }
        const int j_up = j_down + 1;

        std::vector<double> down_prim = cons_to_noncons(initial_cons[i][j_down], g);
        std::vector<double> up_prim = cons_to_noncons(initial_cons[i][j_up], g);
        double c_down = calc_sound_speed(down_prim, g);
        double c_up = calc_sound_speed(up_prim, g);
        auto [v_cont, p_cont] = calc_contact_pressure_velocity_y(down_prim, up_prim, c_down, c_up, g);

        const double y_discont = y_min + (j_down - fict_y + 1) * dy;
        for (int j = fict_y; j < Ny_total - fict_y; ++j) {
            double y = y_min + (j - fict_y + 0.5) * dy;
            double eta = (t > 1e-12) ? (y - y_discont) / t : (y < y_discont ? -1e10 : 1e10);
            std::vector<double> res_1d = find_solution_y(down_prim, up_prim, v_cont, p_cont, c_down, c_up, eta, g);
            for (int ii = fict_x; ii < Nx_total - fict_x; ++ii) {
                analytic_solution[ii][j] = res_1d;
            }
        }
    }

    return analytic_solution;
}


double minmod(double a, double b){
    if (a * b <= 0.0) return 0.0;
    return (fabs(a) < fabs(b)) ? a : b;
}


void Kolgan(const std::vector<double>& v_cons_left, 
           const std::vector<double>& v_cons_center, 
           const std::vector<double>& v_cons_right,
           std::vector<double>& v_cons_left_face, 
           std::vector<double>& v_cons_right_face, 
           double g){
    
    // Р’С‹С‡РёСЃР»СЏРµРј РЅР°РєР»РѕРЅС‹ РґР»СЏ РєР°Р¶РґРѕР№ РєРѕРЅСЃРµСЂРІР°С‚РёРІРЅРѕР№ РїРµСЂРµРјРµРЅРЅРѕР№
    std::vector<double> v_noncons_left = cons_to_noncons(v_cons_left, g);
    std::vector<double> v_noncons_center = cons_to_noncons(v_cons_center, g);
    std::vector<double> v_noncons_right = cons_to_noncons(v_cons_right, g);
    for (int var = 0; var < M; ++var) {
        double a = (v_noncons_center[var] - v_noncons_left[var]);    // РЅР°РєР»РѕРЅ СЃР»РµРІР°
        double b = (v_noncons_right[var] - v_noncons_center[var]);   // РЅР°РєР»РѕРЅ СЃРїСЂР°РІР°
        
        // РњРѕРґРёС„РёС†РёСЂРѕРІР°РЅРЅС‹Р№ РѕРіСЂР°РЅРёС‡РёС‚РµР»СЊ minmod
        //double delta = minmod(a, b);
        double delta = minmod((a + b) / 2.0, 2.0 * minmod(a, b));
        // Р РµРєРѕРЅСЃС‚СЂСѓРєС†РёСЏ Р·РЅР°С‡РµРЅРёР№ РЅР° РіСЂР°РЅСЏС…
        v_noncons_left[var] = v_noncons_center[var] - 0.5 * delta;
        v_noncons_right[var] = v_noncons_center[var] + 0.5 * delta;
    }
    v_cons_left_face = noncons_to_cons(v_noncons_left, g);
    v_cons_right_face = noncons_to_cons(v_noncons_right, g);

}
double mc(double rr) {
    if (rr <= 0.0) return 0.0;
    
    double val1 = 2.0 * rr;
    double val2 = (1.0 + rr) / 2.0;
    
    // Берем максимум из 0, а затем минимум из трех значений
    return std::min({val1, val2, 2.0});
}



void Kolgan_for_Rodionov(const std::vector<double>& v_cons_left, 
           const std::vector<double>& v_cons_center, 
           const std::vector<double>& v_cons_right,
           std::vector<double>& delta, 
           double g){
    
    // Вычисляем наклоны для примитивных переменных
    std::vector<double> v_noncons_left = cons_to_noncons(v_cons_left, g);
    std::vector<double> v_noncons_center = cons_to_noncons(v_cons_center, g);
    std::vector<double> v_noncons_right = cons_to_noncons(v_cons_right, g);
    
    std::vector<double> delta_prim(M, 0.0);
    
    for (int var = 0; var < M; var++) {
        double a = (v_noncons_center[var] - v_noncons_left[var]);
        double b = (v_noncons_right[var] - v_noncons_center[var]);
        
        // Модифицированный ограничитель minmod
        delta_prim[var] = minmod((a + b) / 2.0, 2.0 * minmod(a, b));
        //delta_prim[var] = minmod(a, b);
        
    }
    
    // Вычисляем дельту в консервативных переменных
    std::vector<double> left_prim(M), right_prim(M);
    for (int var = 0; var < M; var++) {
        left_prim[var] = v_noncons_center[var] - 0.5 * delta_prim[var];
        right_prim[var] = v_noncons_center[var] + 0.5 * delta_prim[var];
    }
    
    std::vector<double> left_cons = noncons_to_cons(left_prim, g);
    std::vector<double> right_cons = noncons_to_cons(right_prim, g);
    
    for (int var = 0; var < M; var++) {
        delta[var] =  (right_cons[var] - left_cons[var]);
    }
}


// Обновленная функция Родионова для 2D предиктора
void Rodionov(const std::vector<std::vector<std::vector<double>>>& u_prev,
              std::vector<std::vector<std::vector<double>>>& u_half,
              const std::vector<std::vector<std::vector<double>>>& delta_x,
              const std::vector<std::vector<std::vector<double>>>& delta_y,
              double dt, double dx, double dy, double g, 
              int Nx, int Ny, int fict_x, int fict_y,
              int left_bc_code, int right_bc_code, int up_bc_code, int down_bc_code) {
    
    // Предиктор: u^{n+1/2} = u^n - dt/2 * (div F(u))
    // Вычисляем дивергенцию, используя линейную реконструкцию внутри ячейки
    
    // Проход по внутренней области
    for(int i = fict_x; i < Nx - fict_x; ++i) {
        for(int j = fict_y; j < Ny - fict_y; ++j) {
            
            // --- X DIRECTION ---
            std::vector<double> x_left_val(M), x_right_val(M);
            for (int var = 0; var < M; var++) {
                x_left_val[var]  = u_prev[i][j][var] - 0.5 * delta_x[i][j][var];
                x_right_val[var] = u_prev[i][j][var] + 0.5 * delta_x[i][j][var];
            }
            enforce_physical_state(x_left_val, g);
            enforce_physical_state(x_right_val, g);

            std::vector<double> flux_x_left  = diff_flux_ncons_x(cons_to_noncons(x_left_val, g), g);
            std::vector<double> flux_x_right = diff_flux_ncons_x(cons_to_noncons(x_right_val, g), g);

            // --- Y DIRECTION ---
            std::vector<double> y_down_val(M), y_up_val(M);
            for (int var = 0; var < M; var++) {
                y_down_val[var] = u_prev[i][j][var] - 0.5 * delta_y[i][j][var];
                y_up_val[var]   = u_prev[i][j][var] + 0.5 * delta_y[i][j][var];
            }
            enforce_physical_state(y_down_val, g);
            enforce_physical_state(y_up_val, g);

            std::vector<double> flux_y_down = diff_flux_ncons_y(cons_to_noncons(y_down_val, g), g);
            std::vector<double> flux_y_up   = diff_flux_ncons_y(cons_to_noncons(y_up_val, g), g);
    
            // Обновление на полшага
            for(int k = 0; k < M; ++k) {
                u_half[i][j][k] = u_prev[i][j][k] 
                                - (dt / dx) * (flux_x_right[k] - flux_x_left[k])
                                - (dt / dy) * (flux_y_up[k] - flux_y_down[k]);
            }
            enforce_physical_state(u_half[i][j], g);
        }
    }

    // Установка граничных условий для u_half (по всем 4 сторонам)
    // Лево/Право (Axis 0)
    for(int j = 0; j < Ny; ++j) {
        for(int k = 0; k < fict_x; ++k) {
            u_half[k][j] = boundary(u_half[fict_x][j], left_bc_code, g, 0); 
            u_half[Nx - 1 - k][j] = boundary(u_half[Nx - fict_x - 1][j], right_bc_code, g, 0);
        }
    }
    // Верх/Низ (Axis 1)
    for(int i = 0; i < Nx; ++i) {
        for(int k = 0; k < fict_y; ++k) {
            u_half[i][k] = boundary(u_half[i][fict_y], down_bc_code, g, 1);
            u_half[i][Ny - 1 - k] = boundary(u_half[i][Ny - fict_y - 1], up_bc_code, g, 1);
        }
    }
}

void Rodionov_local(const std::vector<std::vector<std::vector<double>>>& u_prev,
                    std::vector<std::vector<std::vector<double>>>& u_half,
                    const std::vector<std::vector<std::vector<double>>>& delta_x,
                    const std::vector<std::vector<std::vector<double>>>& delta_y,
                    double dt, double dx, double dy, double g,
                    int i_start, int i_end, int j_start, int j_end) {
    for (int i = i_start; i < i_end; ++i) {
        for (int j = j_start; j < j_end; ++j) {
            std::vector<double> x_left_val(M), x_right_val(M);
            for (int var = 0; var < M; ++var) {
                x_left_val[var] = u_prev[i][j][var] - 0.5 * delta_x[i][j][var];
                x_right_val[var] = u_prev[i][j][var] + 0.5 * delta_x[i][j][var];
            }
            enforce_physical_state(x_left_val, g);
            enforce_physical_state(x_right_val, g);

            const std::vector<double> flux_x_left = diff_flux_ncons_x(cons_to_noncons(x_left_val, g), g);
            const std::vector<double> flux_x_right = diff_flux_ncons_x(cons_to_noncons(x_right_val, g), g);

            std::vector<double> y_down_val(M), y_up_val(M);
            for (int var = 0; var < M; ++var) {
                y_down_val[var] = u_prev[i][j][var] - 0.5 * delta_y[i][j][var];
                y_up_val[var] = u_prev[i][j][var] + 0.5 * delta_y[i][j][var];
            }
            enforce_physical_state(y_down_val, g);
            enforce_physical_state(y_up_val, g);

            const std::vector<double> flux_y_down = diff_flux_ncons_y(cons_to_noncons(y_down_val, g), g);
            const std::vector<double> flux_y_up = diff_flux_ncons_y(cons_to_noncons(y_up_val, g), g);

            for (int k = 0; k < M; ++k) {
                u_half[i][j][k] = u_prev[i][j][k]
                                - (dt / dx) * (flux_x_right[k] - flux_x_left[k])
                                - (dt / dy) * (flux_y_up[k] - flux_y_down[k]);
            }
            enforce_physical_state(u_half[i][j], g);
        }
    }
}












std::vector<double> rusanov_flux(const std::vector<double>& left_cons, 
                                 const std::vector<double>& right_cons, 
                                 double g) {
    // 1. Сразу считаем физические потоки F_L и F_R
    // Нам нужны примитивные переменные только для скорости звука и давления
    // Но чтобы не аллоцировать векторы, считаем локально
    
    // Left state extraction
    double rho_L = left_cons[r];
    double u_L   = left_cons[ru] / rho_L;
    double E_L   = left_cons[e];
    double p_L   = (g - 1.0) * (E_L - 0.5 * rho_L * u_L * u_L);
    if(p_L < 1e-8) p_L = 1e-8;
    double c_L   = std::sqrt(g * p_L / rho_L);

    // Right state extraction
    double rho_R = right_cons[r];
    double u_R   = right_cons[ru] / rho_R;
    double E_R   = right_cons[e];
    double p_R   = (g - 1.0) * (E_R - 0.5 * rho_R * u_R * u_R);
    if(p_R < 1e-8) p_R = 1e-8;
    double c_R   = std::sqrt(g * p_R / rho_R);

    // Physical fluxes F(U) = [rho*u, rho*u^2 + p, u*(E+p)]
    double fL_0 = left_cons[ru];
    double fL_1 = left_cons[ru] * u_L + p_L;
    double fL_2 = u_L * (E_L + p_L);

    double fR_0 = right_cons[ru];
    double fR_1 = right_cons[ru] * u_R + p_R;
    double fR_2 = u_R * (E_R + p_R);

    // 2. Спектральный радиус (максимальная скорость волны)
    // S_max = max(|u_L|+c_L, |u_R|+c_R)
    double S_L = std::abs(u_L) + c_L;
    double S_R = std::abs(u_R) + c_R;
    double S_max = std::max(S_L, S_R);

    // 3. Итоговый поток: F_Rus = 0.5 * (F_L + F_R) - 0.5 * S_max * (U_R - U_L)
    std::vector<double> flux(M);
    flux[0] = 0.5 * (fL_0 + fR_0 - S_max * (right_cons[r] - left_cons[r]));
    flux[1] = 0.5 * (fL_1 + fR_1 - S_max * (right_cons[ru] - left_cons[ru]));
    flux[2] = 0.5 * (fL_2 + fR_2 - S_max * (right_cons[e] - left_cons[e]));

    return flux;
}

// --- ОПТИМИЗИРОВАННЫЙ HLLC (Monolithic) ---
std::vector<double> hllc_flux_new(const std::vector<double>& left_cons, 
                              const std::vector<double>& right_cons, 
                              double g, int axis) {
    // Распаковка переменных (без создания лишних векторов)
    double rho_L = left_cons[r];
    double u_L   = left_cons[ru] / rho_L;
    double v_L   = left_cons[rv] / rho_L;
    double E_L   = left_cons[e];
    double p_L   = (g - 1.0) * (E_L - 0.5 * rho_L * (u_L * u_L + v_L * v_L));
    if(p_L < 1e-8) p_L = 1e-8;
    double a_L   = std::sqrt(g * p_L / rho_L);

    double rho_R = right_cons[r];
    double u_R   = right_cons[ru] / rho_R;
    double v_R   = right_cons[rv] / rho_R;
    double E_R   = right_cons[e];
    double p_R   = (g - 1.0) * (E_R - 0.5 * rho_R * (u_R * u_R + v_R * v_R));
    if(p_R < 1e-8) p_R = 1e-8;
    double a_R   = std::sqrt(g * p_R / rho_R);

    // --- Шаг 1: Оценка скоростей волн (PVRS / Hybrid) ---
    double rho_bar = 0.5 * (rho_L + rho_R);
    double a_bar   = 0.5 * (a_L + a_R);
    
    // PVRS оценки
    double p_pvrs, p_star;
    if (axis == 0){
        p_pvrs = 0.5 * (p_L + p_R) - 0.5 * (u_R - u_L) * rho_bar * a_bar;
        p_star = std::max(0.0, p_pvrs); 
    }
    else {
        p_pvrs = 0.5 * (p_L + p_R) - 0.5 * (v_R - v_L) * rho_bar * a_bar;
        p_star = std::max(0.0, p_pvrs); 
    }

    // Q-факторы для волновых скоростей
    double qL = (p_star <= p_L) ? 1.0 : std::sqrt(1.0 + (g + 1.0) / (2.0 * g) * (p_star / p_L - 1.0));
    double qR = (p_star <= p_R) ? 1.0 : std::sqrt(1.0 + (g + 1.0) / (2.0 * g) * (p_star / p_R - 1.0));
    double S_L, S_R, numer, denom;
    if (axis == 0){
        S_L = u_L - a_L * qL;
        S_R = u_R + a_R * qR;
        numer = p_R - p_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R);
        denom = rho_L * (S_L - u_L) - rho_R * (S_R - u_R);
    }

    else{
        S_L = v_L - a_L * qL;
        S_R = v_R + a_R * qR;
        numer = p_R - p_L + rho_L * v_L * (S_L - v_L) - rho_R * v_R * (S_R - v_R);
        denom = rho_L * (S_L - v_L) - rho_R * (S_R - v_R);
        }
    // Контактная скорость (S_star)
    
    double S_star = numer / denom;

    // --- Шаг 2: Выбор потока ---
    std::vector<double> flux(M);

    if (S_L >= 0.0) {
        // F_L
        if(axis == 0){
            flux[0] = left_cons[ru];
            flux[1] = left_cons[ru] * u_L + p_L;
            flux[2] = left_cons[ru] * v_L;
            flux[3] = u_L * (E_L + p_L);
        }
        else{
            flux[0] = left_cons[rv];
            flux[1] = left_cons[rv] * u_L;
            flux[2] = left_cons[rv] * v_L + p_L;
            flux[3] = v_L * (E_L + p_L);

        }
    } 
    else if (S_R <= 0.0) {
        // F_R
        if(axis == 0){
            flux[0] = right_cons[ru];
            flux[1] = right_cons[ru] * u_R + p_R;
            flux[2] = right_cons[ru] * v_R;
            flux[3] = u_R * (E_R + p_R);
        }

        else{
            flux[0] = right_cons[rv];
            flux[1] = right_cons[rv] * u_R;
            flux[2] = right_cons[rv] * v_R + p_R;
            flux[3] = v_R * (E_R + p_R); 
        }
    } 
    else {
        // HLLC поток (L* или R*)
        double S_K, u_K, v_K, p_K, rho_K, E_K; 
        double factor; // rho * (S_K - u_K) / (S_K - S_star)

        if (S_L <= 0.0 && S_star >= 0.0) {
            // Left Star State
            S_K = S_L; 
            u_K = u_L;
            v_K = v_L; 
            p_K = p_L; 
            rho_K = rho_L;
            E_K = E_L;
        } else {
            // Right Star State
            S_K = S_R; 
            u_K = u_R;
            v_K = v_R; 
            p_K = p_R; 
            rho_K = rho_R;
            E_K = E_R;
        }

        // Физический поток F_K
        double F_K_0, F_K_1, F_K_2, F_K_3;
        if (axis == 0){
            F_K_0 = rho_K * u_K;
            F_K_1 = rho_K * u_K * u_K + p_K;
            F_K_2 = rho_K * u_K * v_K;
            F_K_3 = u_K * (E_K + p_K);
        }

        else{
            F_K_0 = rho_K * v_K;
            F_K_1 = rho_K * u_K * v_K;
            F_K_2 = rho_K * v_K * v_K + p_K;
            F_K_3 = v_K * (E_K + p_K);
        }

        // Консервативные переменные U_star_K
        if (axis == 0){
            factor = rho_K * (S_K - u_K) / (S_K - S_star);
            double U_star_0 = factor;
            double U_star_1 = factor * S_star;
            double U_star_2 = factor * v_K;
            double U_star_3 = factor * (E_K / rho_K + (S_star - u_K) * (S_star + p_K / (rho_K * (S_K - u_K))));

        // HLLC Flux = F_K + S_K * (U_star_K - U_K)
            flux[0] = F_K_0 + S_K * (U_star_0 - rho_K);
            flux[1] = F_K_1 + S_K * (U_star_1 - rho_K * u_K);
            flux[2] = F_K_2 + S_K * (U_star_2 - rho_K * v_K);
            flux[3] = F_K_3 + S_K * (U_star_3 - E_K);
        }
        else {
            factor = rho_K * (S_K - v_K) / (S_K - S_star);
            double U_star_0 = factor;
            double U_star_1 = factor * u_K;
            double U_star_2 = factor * S_star;
            double U_star_3 = factor * (E_K / rho_K + (S_star - v_K) * (S_star + p_K / (rho_K * (S_K - v_K))));

        // HLLC Flux = F_K + S_K * (U_star_K - U_K)
            flux[0] = F_K_0 + S_K * (U_star_0 - rho_K);
            flux[1] = F_K_1 + S_K * (U_star_1 - rho_K * u_K);
            flux[2] = F_K_2 + S_K * (U_star_2 - rho_K * v_K);
            flux[3] = F_K_3 + S_K * (U_star_3 - E_K);

        }
    }

    return flux;
}







// --- ОПТИМИЗИРОВАННЫЙ HLL ---
std::vector<double> hll_flux_new(const std::vector<double>& left_cons, 
                             const std::vector<double>& right_cons, 
                             double g, int axis) {
    // 1. Распаковка левого состояния
    double rho_L = left_cons[r];
    double u_L   = left_cons[ru] / rho_L;
    double v_L = left_cons[rv] / rho_L;
    double E_L   = left_cons[e];
    double p_L   = (g - 1.0) * (E_L - 0.5 * rho_L * (u_L * u_L + v_L * v_L));
    if(p_L < 1e-8) p_L = 1e-8;
    double a_L   = std::sqrt(g * p_L / rho_L);

    // 2. Распаковка правого состояния
    double rho_R = right_cons[r];
    double u_R   = right_cons[ru] / rho_R;
    double v_R   = right_cons[rv] / rho_R;
    double E_R   = right_cons[e];
    double p_R   = (g - 1.0) * (E_R - 0.5 * rho_R * (u_R * u_R + v_R * v_R));
    if(p_R < 1e-8) p_R = 1e-8;
    double a_R   = std::sqrt(g * p_R / rho_R);

    // 3. Оценка волновых скоростей (Einfeldt or primitive based)
    // Используем простую оценку как в оригинале: min(uL-aL, uR-aR)

    double FL_0, FL_1, FL_2, FL_3, FR_0, FR_1, FR_2, FR_3;
    double S_R, S_L;
    if (axis == 0){
        S_L = std::min(u_L - a_L, u_R - a_R);
        S_R = std::max(u_L + a_L, u_R + a_R);
        FL_0 = left_cons[ru];
        FL_1 = left_cons[ru] * u_L + p_L;
        FL_2 = left_cons[ru] * v_L; 
        FL_3 = u_L * (E_L + p_L);

        FR_0 = right_cons[ru];
        FR_1 = right_cons[ru] * u_R + p_R;
        FR_2 = right_cons[ru] * v_R; 
        FR_3 = u_R * (E_R + p_R);

    }
    else {
        S_L = std::min(v_L - a_L, v_R - a_R);
        S_R = std::max(v_L + a_L, v_R + a_R);
        FL_0 = left_cons[rv];
        FL_1 = left_cons[rv] * u_L;
        FL_2 = left_cons[rv] * v_L + p_L;
        FL_3 = v_L * (E_L + p_L);

        FR_0 = right_cons[rv];
        FR_1 = right_cons[rv] * u_R;
        FR_2 = right_cons[rv] * v_R + p_R;
        FR_3 = v_R * (E_R + p_R);

    }
    

  

    // 5. Выбор потока
    std::vector<double> flux(M);

    if (S_L >= 0.0) {
        flux[0] = FL_0;
        flux[1] = FL_1;
        flux[2] = FL_2;
        flux[3] = FL_3;
    } 
    else if (S_R <= 0.0) {
        flux[0] = FR_0;
        flux[1] = FR_1;
        flux[2] = FR_2;
        flux[3] = FR_3;
    } 
    else {
        // HLL Formula: F_hll = (Sr*Fl - Sl*Fr + Sl*Sr*(Ur - Ul)) / (Sr - Sl)
        double inv_denom = 1.0 / (S_R - S_L);
        double term1 = S_L * S_R; 
        if(axis == 0){ 
        
            flux[0] = (S_R * FL_0 - S_L * FR_0 + term1 * (right_cons[r] - left_cons[r])) * inv_denom;
            flux[1] = (S_R * FL_1 - S_L * FR_1 + term1 * (right_cons[ru] - left_cons[ru])) * inv_denom;
            flux[2] = (S_R * FL_2 - S_L * FR_2 + term1 * (right_cons[rv] - left_cons[rv])) * inv_denom;
            flux[3] = (S_R * FL_3 - S_L * FR_3 + term1 * (right_cons[e] - left_cons[e])) * inv_denom;
        }
        else{
            flux[0] = (S_R * FL_0 - S_L * FR_0 + term1 * (right_cons[r] - left_cons[r])) * inv_denom;
            flux[1] = (S_R * FL_1 - S_L * FR_1 + term1 * (right_cons[ru] - left_cons[ru])) * inv_denom;
            flux[2] = (S_R * FL_2 - S_L * FR_2 + term1 * (right_cons[rv]  - left_cons[rv])) * inv_denom;
            flux[3] = (S_R * FL_3 - S_L * FR_3 + term1 * (right_cons[e] - left_cons[e])) * inv_denom;

        }
    }

    return flux;
}


struct OsherStates2D {
    double p_star, u_star;
    double c_1_3, c_2_3;
    // Состояния в примитивах: {rho, u, v, p}
    double U_1_3[4], U_2_3[4], U_s0[4], U_s1[4];
};

// Быстрая конвертация 2D cons -> prim
inline void cons_to_prim_2d(const std::vector<double>& v_cons, double g, double out[4]) {
    const double rho = v_cons[0];
    const double inv_rho = 1.0 / rho;
    const double u = v_cons[1] * inv_rho;
    const double v = v_cons[2] * inv_rho;
    double p = (g - 1.0) * (v_cons[3] - 0.5 * rho * (u * u + v * v));
    
    if (p < 1e-10) p = 1e-10; // P_FLOOR
    out[0] = rho; out[1] = u; out[2] = v; out[3] = p;
}

// Физический поток 2D (dir: 0 - X, 1 - Y)
inline void flux_from_prim_2d(const double q[4], double g, int dir, double out[4]) {
    const double rho = q[0], u = q[1], v = q[2], p = q[3];
    const double e_kin = 0.5 * rho * (u * u + v * v);
    const double energy = p / (g - 1.0) + e_kin;

    if (dir == 0) { // Поток по X
        out[0] = rho * u;
        out[1] = rho * u * u + p;
        out[2] = rho * u * v;
        out[3] = u * (energy + p);
    } else { // Поток по Y
        out[0] = rho * v;
        out[1] = rho * v * u;
        out[2] = rho * v * v + p;
        out[3] = v * (energy + p);
    }
}

OsherStates2D calc_osher_states_2d(const double qL[4], const double qR[4], double g, int axis) {
    OsherStates2D res;
    
    // Определяем нормальные и тангенциальные составляющие
    double unL = (axis == 0) ? qL[1] : qL[2];
    double utL = (axis == 0) ? qL[2] : qL[1];
    double unR = (axis == 0) ? qR[1] : qR[2];
    double utR = (axis == 0) ? qR[2] : qR[1];

    double cL = std::sqrt(g * qL[3] / qL[0]);
    double cR = std::sqrt(g * qR[3] / qR[0]);
    double z = (g - 1.0) / (2.0 * g);

    // 1. Звездное состояние (p* и u*)
    double numer = cL + cR - 0.5 * (g - 1.0) * (unR - unL);
    double denom = cL / std::pow(qL[3], z) + cR / std::pow(qR[3], z);
    if (numer < 0) numer = 0;
    
    res.p_star = std::pow(numer / denom, 1.0 / z);
    res.c_1_3 = cL * std::pow(res.p_star / qL[3], z);
    res.u_star = unL + 2.0 / (g - 1.0) * (cL - res.c_1_3);
    res.c_2_3 = cR * std::pow(res.p_star / qR[3], z);

    // 2. Промежуточные 1/3 и 2/3 (utL до контакта, utR после)
    auto fill_state = [&](double q_out[4], double rho, double un, double ut, double p) {
        q_out[0] = rho; q_out[3] = p;
        if (axis == 0) { q_out[1] = un; q_out[2] = ut; }
        else           { q_out[1] = ut; q_out[2] = un; }
    };

    fill_state(res.U_1_3, g * res.p_star / (res.c_1_3 * res.c_1_3), res.u_star, utL, res.p_star);
    fill_state(res.U_2_3, g * res.p_star / (res.c_2_3 * res.c_2_3), res.u_star, utR, res.p_star);

    // 3. Звуковые точки
    // Sonic L (un - c = 0)
    double JpL = unL + 2.0 * cL / (g - 1.0);
    double cs0 = JpL * (g - 1.0) / (g + 1.0);
    double rhos0 = qL[0] * std::pow(cs0 / cL, 2.0 / (g - 1.0));
    double ps0 = qL[3] * std::pow(rhos0 / qL[0], g);
    fill_state(res.U_s0, rhos0, cs0, utL, ps0);

    // Sonic R (un + c = 0)
    double JmR = unR - 2.0 * cR / (g - 1.0);
    double cs1 = -JmR * (g - 1.0) / (g + 1.0);
    double rhos1 = qR[0] * std::pow(cs1 / cR, 2.0 / (g - 1.0));
    double ps1 = qR[3] * std::pow(rhos1 / qR[0], g);
    fill_state(res.U_s1, rhos1, -cs1, utR, ps1);

    return res;
}

std::vector<double> osher_flux_2d(const std::vector<double>& left_cons, 
                                  const std::vector<double>& right_cons, 
                                  double g, int axis) {
    double qL[4], qR[4];
    cons_to_prim_2d(left_cons, g, qL);
    cons_to_prim_2d(right_cons, g, qR);

    OsherStates2D st = calc_osher_states_2d(qL, qR, g, axis);
    
    double F_net[4];
    flux_from_prim_2d(qL, g, axis, F_net); // Старт с F(U_L)

    double F_L[4], F_1_3[4], F_2_3[4], F_R[4], F_s0[4], F_s1[4];
    flux_from_prim_2d(qL, g, axis, F_L);
    flux_from_prim_2d(st.U_1_3, g, axis, F_1_3);
    flux_from_prim_2d(st.U_2_3, g, axis, F_2_3);
    flux_from_prim_2d(qR, g, axis, F_R);
    flux_from_prim_2d(st.U_s0, g, axis, F_s0);
    flux_from_prim_2d(st.U_s1, g, axis, F_s1);

    // Характеристические скорости
    double cL = std::sqrt(g * qL[3] / qL[0]);
    double cR = std::sqrt(g * qR[3] / qR[0]);
    double unL = (axis == 0) ? qL[1] : qL[2];
    double unR = (axis == 0) ? qR[1] : qR[2];

    double lamL = unL - cL;
    double lam13 = st.u_star - st.c_1_3;
    double lam_contact = st.u_star;
    double lam23 = st.u_star + st.c_2_3;
    double lamR = unR + cR;

    // Интегрирование (сумма приращений)
    auto add_diff = [&](const double F_high[4], const double F_low[4]) {
        for(int i=0; i<4; i++) F_net[i] += (F_high[i] - F_low[i]);
    };

    // 1. Волна un - c
    if (lamL < 0) {
        if (lam13 < 0) add_diff(F_1_3, F_L);
        else           add_diff(F_s0, F_L);
    }
    if (lamL > 0 && lam13 < 0) add_diff(F_1_3, F_s0);

    // 2. Контакт (un)
    if (lam_contact < 0) add_diff(F_2_3, F_1_3);

    // 3. Волна un + c
    if (lam23 < 0) {
        if (lamR < 0) add_diff(F_R, F_2_3);
        else          add_diff(F_s1, F_2_3);
    }
    if (lam23 > 0 && lamR < 0) add_diff(F_R, F_s1);

    return {F_net[0], F_net[1], F_net[2], F_net[3]};
}






// Вспомогательная функция для получения примитивных переменных
struct Prim {
    double rho, u, v, p, a, H;
};

Prim getPrim(const std::vector<double>& Z, double g) {
    double rho = Z[0];
    double u = Z[1] / rho;
    double v = Z[2] / rho;
    double e_kin = 0.5 * (u * u + v * v);
    double p = (g - 1.0) * (Z[3] - rho * e_kin);
    double a = std::sqrt(g * p / rho);
    double H = (Z[3] + p) / rho;
    return {rho, u, v, p, a, H};
}

// Физический поток в 2D (dir: 0 для X, 1 для Y)
std::vector<double> Flux(const std::vector<double>& Z, double g, int dir) {
    Prim p = getPrim(Z, g);
    if (dir == 0) { // Поток по X
        return {Z[1], Z[1] * p.u + p.p, Z[1] * p.v, (Z[3] + p.p) * p.u};
    } else { // Поток по Y
        return {Z[2], Z[2] * p.u, Z[2] * p.v + p.p, (Z[3] + p.p) * p.v};
    }
}

// --- RUSANOV (Local Lax-Friedrichs) ---
std::vector<double> rusanov_2d(const std::vector<double>& UL, const std::vector<double>& UR, double g, int dir) {
    Prim pL = getPrim(UL, g);
    Prim pR = getPrim(UR, g);
    
    std::vector<double> FL = Flux(UL, g, dir);
    std::vector<double> FR = Flux(UR, g, dir);
    
    double velL = (dir == 0) ? pL.u : pL.v;
    double velR = (dir == 0) ? pR.u : pR.v;
    
    double s_max = std::max(std::abs(velL) + pL.a, std::abs(velR) + pR.a);
    
    std::vector<double> F(4);
    for(int i = 0; i < 4; ++i) {
        F[i] = 0.5 * (FL[i] + FR[i]) - 0.5 * s_max * (UR[i] - UL[i]);
    }
    return F;
}


// Структура для примитивных переменных
struct State2D {
    double rho, u, v, p, a, H;
};

// Получение примитивных переменных и энтальпии
State2D getState(const std::vector<double>& Z, double g) {
    double rho = Z[0];
    double u = Z[1] / rho;
    double v = Z[2] / rho;
    double e_kin = 0.5 * (u * u + v * v);
    double p = (g - 1.0) * (Z[3] - rho * e_kin);
    double a = std::sqrt(g * p / rho);
    double H = (Z[3] + p) / rho;
    return {rho, u, v, p, a, H};
}

// Физический поток вдоль нормали (dir=0 -> X, dir=1 -> Y)
std::vector<double> PhysFlux(const std::vector<double>& Z, double g, int dir) {
    State2D s = getState(Z, g);
    if (dir == 0) return {Z[1], Z[1] * s.u + s.p, Z[1] * s.v, (Z[3] + s.p) * s.u};
    return {Z[2], Z[2] * s.u, Z[2] * s.v + s.p, (Z[3] + s.p) * s.v};
}


std::vector<double> roe_2d(const std::vector<double>& UL, const std::vector<double>& UR, double g, int dir) {
    State2D sL = getState(UL, g);
    State2D sR = getState(UR, g);

    // Осреднение Роэ
    double wl = std::sqrt(sL.rho);
    double wr = std::sqrt(sR.rho);
    double rho = wl * wr;
    double u = (wl * sL.u + wr * sR.u) / (wl + wr);
    double v = (wl * sL.v + wr * sR.v) / (wl + wr);
    double H = (wl * sL.H + wr * sR.H) / (wl + wr);
    double a = std::sqrt((g - 1.0) * (H - 0.5 * (u * u + v * v)));

    // Выбор нормальной (un) и тангенциальной (ut) скоростей
    double un = (dir == 0) ? u : v;
    double ut = (dir == 0) ? v : u;

    // Разности консервативных переменных
    std::vector<double> dU(4);
    for(int i=0; i<4; i++) dU[i] = UR[i] - UL[i];

    // Собственные значения
    double l1 = un - a, l2 = un, l3 = un, l4 = un + a;
    
    // Амплитуды волн (alpha)
    double dp = sR.p - sL.p;
    double drho = sR.rho - sL.rho;
    double dun = (dir == 0 ? sR.u - sL.u : sR.v - sL.v);
    double dut = (dir == 0 ? sR.v - sL.v : sR.u - sL.u);

    double a1 = (dp - rho * a * dun) / (2.0 * a * a);
    double a4 = (dp + rho * a * dun) / (2.0 * a * a);
    double a2 = drho - dp / (a * a);
    double a3 = rho * dut;

    // Векторы собственных состояний (K_i) и диссипация
    std::vector<double> diss(4, 0.0);
    auto add_diss = [&](double lam, double alpha, std::vector<double> K) {
        for(int i=0; i<4; i++) diss[i] += std::abs(lam) * alpha * K[i];
    };

    if (dir == 0) { // Нормаль по X
        add_diss(l1, a1, {1.0, un - a, ut, H - un * a});
        add_diss(l2, a2, {1.0, un, ut, 0.5 * (u * u + v * v)});
        add_diss(l3, a3, {0.0, 0.0, 1.0, ut});
        add_diss(l4, a4, {1.0, un + a, ut, H + un * a});
    } else { // Нормаль по Y
        add_diss(l1, a1, {1.0, ut, un - a, H - un * a});
        add_diss(l2, a2, {1.0, ut, un, 0.5 * (u * u + v * v)});
        add_diss(l3, a3, {0.0, 1.0, 0.0, ut});
        add_diss(l4, a4, {1.0, ut, un + a, H + un * a});
    }

    std::vector<double> FL = PhysFlux(UL, g, dir);
    std::vector<double> FR = PhysFlux(UR, g, dir);
    std::vector<double> Flux(4);
    for(int i=0; i<4; i++) Flux[i] = 0.5 * (FL[i] + FR[i] - diss[i]);

    return Flux;
}



// --- ОПТИМИЗИРОВАННЫЙ ROE 2D ---
std::vector<double> roe_flux_2d(const std::vector<double>& left_cons, 
                               const std::vector<double>& right_cons, 
                               double g, int axis) {
    // Индексы: 0: rho, 1: rho*u, 2: rho*v, 3: E
    const int RU = 1, RV = 2, ENER = 3;

    // 1. Примитивные переменные и Энтальпия
    double rhoL = left_cons[RHO];
    double uL = left_cons[RU] / rhoL;
    double vL = left_cons[RV] / rhoL;
    double pL = (g - 1.0) * (left_cons[ENER] - 0.5 * rhoL * (uL * uL + vL * vL));
    if(pL < 1e-8) pL = 1e-8;
    double HL = (left_cons[ENER] + pL) / rhoL;

    double rhoR = right_cons[RHO];
    double uR = right_cons[RU] / rhoR;
    double vR = right_cons[RV] / rhoR;
    double pR = (g - 1.0) * (right_cons[ENER] - 0.5 * rhoR * (uR * uR + vR * vR));
    if(pR < 1e-8) pR = 1e-8;
    double HR = (right_cons[ENER] + pR) / rhoR;

    // 2. Осреднение по Роэ
    double sqL = std::sqrt(rhoL);
    double sqR = std::sqrt(rhoR);
    double inv_denom = 1.0 / (sqL + sqR);

    double u_tilde = (sqL * uL + sqR * uR) * inv_denom;
    double v_tilde = (sqL * vL + sqR * vR) * inv_denom;
    double H_tilde = (sqL * HL + sqR * HR) * inv_denom;
    
    double vel2 = u_tilde * u_tilde + v_tilde * v_tilde;
    double a2_tilde = (g - 1.0) * (H_tilde - 0.5 * vel2);
    if (a2_tilde <= 0.0) a2_tilde = 1e-8;
    double a_tilde = std::sqrt(a2_tilde);

    // Выбор нормальной и тангенциальной скорости для расчетов волн
    double un_tilde = (axis == 0) ? u_tilde : v_tilde;
    double ut_tilde = (axis == 0) ? v_tilde : u_tilde;

    // 3. Собственные значения и Энтропийная поправка
    double l1 = un_tilde - a_tilde;
    double l2 = un_tilde; // Для энтропии
    double l3 = un_tilde; // Для тангенциальной скорости (shear wave)
    double l4 = un_tilde + a_tilde;

    double delta_e = 0.2 * (std::abs(un_tilde) + a_tilde);
    auto fix = [&](double l) {
        double abs_l = std::abs(l);
        return (abs_l < delta_e) ? (l*l + delta_e*delta_e)/(2.0*delta_e) : abs_l;
    };
    
    double abs_l1 = fix(l1);
    double abs_l2 = fix(l2);
    double abs_l3 = fix(l3);
    double abs_l4 = fix(l4);

    // 4. Скачки (Deltas) и Интенсивности волн (Alphas)
    double d_rho = rhoR - rhoL;
    double d_u   = uR - uL;
    double d_v   = vR - vL;
    double d_p   = pR - pL;
    
    double d_un = (axis == 0) ? d_u : d_v;
    double d_ut = (axis == 0) ? d_v : d_u;

    double rho_tilde = sqL * sqR;

    double alpha1 = (d_p - rho_tilde * a_tilde * d_un) / (2.0 * a2_tilde);
    double alpha2 = d_rho - d_p / a2_tilde;
    double alpha3 = rho_tilde * d_ut; // Интенсивность сдвиговой волны
    double alpha4 = (d_p + rho_tilde * a_tilde * d_un) / (2.0 * a2_tilde);

    // 5. Диссипация (Sum(|l| * alpha * K))
    // K - правые собственные векторы
    double diss[4] = {0, 0, 0, 0};

    auto add_to_diss = [&](double abs_l, double alpha, double k0, double k1, double k2, double k3) {
        diss[0] += abs_l * alpha * k0;
        diss[1] += abs_l * alpha * k1;
        diss[2] += abs_l * alpha * k2;
        diss[3] += abs_l * alpha * k3;
    };

    if (axis == 0) { // Поток по X
        add_to_diss(abs_l1, alpha1, 1.0, u_tilde - a_tilde, v_tilde, H_tilde - u_tilde * a_tilde);
        add_to_diss(abs_l2, alpha2, 1.0, u_tilde, v_tilde, 0.5 * vel2);
        add_to_diss(abs_l3, alpha3, 0.0, 0.0, 1.0, v_tilde); // Сдвиг по V
        add_to_diss(abs_l4, alpha4, 1.0, u_tilde + a_tilde, v_tilde, H_tilde + u_tilde * a_tilde);
    } else { // Поток по Y
        add_to_diss(abs_l1, alpha1, 1.0, u_tilde, v_tilde - a_tilde, H_tilde - v_tilde * a_tilde);
        add_to_diss(abs_l2, alpha2, 1.0, u_tilde, v_tilde, 0.5 * vel2);
        add_to_diss(abs_l3, alpha3, 0.0, 1.0, 0.0, u_tilde); // Сдвиг по U
        add_to_diss(abs_l4, alpha4, 1.0, u_tilde, v_tilde + a_tilde, H_tilde + v_tilde * a_tilde);
    }

    // 6. Физические потоки
    double FL[4], FR[4];
    auto get_phys_flux = [&](const std::vector<double>& cons, double p, double u, double v, double out[4]) {
        if (axis == 0) {
            out[0] = cons[RU];
            out[1] = cons[RU] * u + p;
            out[2] = cons[RU] * v;
            out[3] = u * (cons[ENER] + p);
        } else {
            out[0] = cons[RV];
            out[1] = cons[RV] * u;
            out[2] = cons[RV] * v + p;
            out[3] = v * (cons[ENER] + p);
        }
    };

    get_phys_flux(left_cons, pL, uL, vL, FL);
    get_phys_flux(right_cons, pR, uR, vR, FR);

    // Сборка итогового потока
    std::vector<double> flux(4);
    for(int i = 0; i < 4; ++i) {
        flux[i] = 0.5 * (FL[i] + FR[i] - diss[i]);
    }

    return flux;
}





// --- ОПТИМИЗИРОВАННЫЙ ROE ---
std::vector<double> roe_flux_new(const std::vector<double>& left_cons, 
                             const std::vector<double>& right_cons, 
                             double g) {
    // 1. Primitive variables & Enthalpy
    double rho_L = left_cons[r];
    double u_L   = left_cons[ru] / rho_L;
    double E_L   = left_cons[e];
    double p_L   = (g - 1.0) * (E_L - 0.5 * rho_L * u_L * u_L);
    if(p_L < 1e-8) p_L = 1e-8;
    double H_L   = (E_L + p_L) / rho_L;

    double rho_R = right_cons[r];
    double u_R   = right_cons[ru] / rho_R;
    double E_R   = right_cons[e];
    double p_R   = (g - 1.0) * (E_R - 0.5 * rho_R * u_R * u_R);
    if(p_R < 1e-8) p_R = 1e-8;
    double H_R   = (E_R + p_R) / rho_R;

    // 2. Roe Averages
    double sqL = std::sqrt(rho_L);
    double sqR = std::sqrt(rho_R);
    double inv_denom = 1.0 / (sqL + sqR);

    double u_tilde = (sqL * u_L + sqR * u_R) * inv_denom;
    double H_tilde = (sqL * H_L + sqR * H_R) * inv_denom;
    
    // Sound speed squared
    double a2_tilde = (g - 1.0) * (H_tilde - 0.5 * u_tilde * u_tilde);
    if (a2_tilde <= 0.0) a2_tilde = 1e-8;
    double a_tilde = std::sqrt(a2_tilde);

    // 3. Wave Speeds (Eigenvalues)
    double l1 = u_tilde - a_tilde;
    double l2 = u_tilde;
    double l3 = u_tilde + a_tilde;

    // Entropy Fix (Harten) - inline optimization
    double delta_e = 0.2 * (std::abs(u_tilde) + a_tilde);
    if(delta_e < 1e-6) delta_e = 1e-6;
    
    auto fix = [&](double l) {
        double abs_l = std::abs(l);
        return (abs_l < delta_e) ? (l*l + delta_e*delta_e)/(2.0*delta_e) : abs_l;
    };
    
    double abs_l1 = fix(l1);
    double abs_l2 = fix(l2);
    double abs_l3 = fix(l3);

    // 4. Wave Strengths (Alpha)
    double d_rho = rho_R - rho_L;
    double d_u   = u_R - u_L;
    double d_p   = p_R - p_L;
    
    double rho_tilde = sqL * sqR; // Geometric mean approximation for strength calc

    double alpha1 = (d_p - rho_tilde * a_tilde * d_u) / (2.0 * a2_tilde);
    double alpha2 = d_rho - d_p / a2_tilde;
    double alpha3 = (d_p + rho_tilde * a_tilde * d_u) / (2.0 * a2_tilde);

    // 5. Right Eigenvectors (K) components
    // K1: {1, u-a, H-ua}
    double k1_0 = 1.0; 
    double k1_1 = u_tilde - a_tilde; 
    double k1_2 = H_tilde - u_tilde * a_tilde;

    // K2: {1, u, 0.5u^2}
    double k2_0 = 1.0; 
    double k2_1 = u_tilde; 
    double k2_2 = 0.5 * u_tilde * u_tilde;

    // K3: {1, u+a, H+ua}
    double k3_0 = 1.0; 
    double k3_1 = u_tilde + a_tilde; 
    double k3_2 = H_tilde + u_tilde * a_tilde;

    // 6. Final Flux: 0.5 * (FL + FR - Sum(|l| * alpha * K))
    
    // Physical fluxes
    double FL_0 = left_cons[ru];
    double FL_1 = left_cons[ru] * u_L + p_L;
    double FL_2 = u_L * (E_L + p_L);

    double FR_0 = right_cons[ru];
    double FR_1 = right_cons[ru] * u_R + p_R;
    double FR_2 = u_R * (E_R + p_R);

    // Dissipation term
    double diss_0 = abs_l1 * alpha1 * k1_0 + abs_l2 * alpha2 * k2_0 + abs_l3 * alpha3 * k3_0;
    double diss_1 = abs_l1 * alpha1 * k1_1 + abs_l2 * alpha2 * k2_1 + abs_l3 * alpha3 * k3_1;
    double diss_2 = abs_l1 * alpha1 * k1_2 + abs_l2 * alpha2 * k2_2 + abs_l3 * alpha3 * k3_2;

    std::vector<double> flux(M);
    flux[0] = 0.5 * (FL_0 + FR_0 - diss_0);
    flux[1] = 0.5 * (FL_1 + FR_1 - diss_1);
    flux[2] = 0.5 * (FL_2 + FR_2 - diss_2);

    return flux;
}

void get_subdomain_bounds(int N, int parts, int rank, int& start, int& end) {
    const int base_count = N / parts;
    const int remainder = N % parts;

    if (rank < remainder) {
        start = rank * (base_count + 1);
        end = start + (base_count + 1);
    } else {
        start = rank * base_count + remainder;
        end = start + base_count;
    }
}

void exchange_halos_global_width(std::vector<std::vector<std::vector<double>>>& u,
                                 MPI_Comm cart_comm, int halo_x, int halo_y,
                                 int i_start, int i_end, int j_start, int j_end,
                                 int left_rank, int right_rank, int down_rank, int up_rank) {
    MPI_Status status;
    auto& scratch = halo_exchange_scratch();

    const int y_len = j_end - j_start;
    const int buf_size_x = halo_x * y_len * M;
    scratch.send_left.resize(buf_size_x);
    scratch.recv_left.resize(buf_size_x);
    scratch.send_right.resize(buf_size_x);
    scratch.recv_right.resize(buf_size_x);

    if (left_rank != MPI_PROC_NULL) {
        int idx_sl = 0;
        for (int i = 0; i < halo_x; ++i) {
            for (int j = j_start; j < j_end; ++j) {
                for (int k = 0; k < M; ++k) {
                    scratch.send_left[idx_sl++] = u[i_start + i][j][k];
                }
            }
        }
    }

    if (right_rank != MPI_PROC_NULL) {
        int idx_sr = 0;
        for (int i = 0; i < halo_x; ++i) {
            for (int j = j_start; j < j_end; ++j) {
                for (int k = 0; k < M; ++k) {
                    scratch.send_right[idx_sr++] = u[i_end - halo_x + i][j][k];
                }
            }
        }
    }

    MPI_Sendrecv(scratch.send_left.data(), buf_size_x, MPI_DOUBLE, left_rank, 0,
                 scratch.recv_right.data(), buf_size_x, MPI_DOUBLE, right_rank, 0,
                 cart_comm, &status);
    MPI_Sendrecv(scratch.send_right.data(), buf_size_x, MPI_DOUBLE, right_rank, 1,
                 scratch.recv_left.data(), buf_size_x, MPI_DOUBLE, left_rank, 1,
                 cart_comm, &status);

    if (left_rank != MPI_PROC_NULL) {
        int idx_rl = 0;
        for (int i = 0; i < halo_x; ++i) {
            for (int j = j_start; j < j_end; ++j) {
                for (int k = 0; k < M; ++k) {
                    u[i_start - halo_x + i][j][k] = scratch.recv_left[idx_rl++];
                }
            }
        }
    }

    if (right_rank != MPI_PROC_NULL) {
        int idx_rr = 0;
        for (int i = 0; i < halo_x; ++i) {
            for (int j = j_start; j < j_end; ++j) {
                for (int k = 0; k < M; ++k) {
                    u[i_end + i][j][k] = scratch.recv_right[idx_rr++];
                }
            }
        }
    }

    const int x_len = (i_end - i_start)
        + (left_rank != MPI_PROC_NULL ? halo_x : 0)
        + (right_rank != MPI_PROC_NULL ? halo_x : 0);
    const int x_start_exchange = (left_rank != MPI_PROC_NULL) ? i_start - halo_x : i_start;
    const int x_end_exchange = (right_rank != MPI_PROC_NULL) ? i_end + halo_x : i_end;

    const int buf_size_y = halo_y * x_len * M;
    scratch.send_down.resize(buf_size_y);
    scratch.recv_down.resize(buf_size_y);
    scratch.send_up.resize(buf_size_y);
    scratch.recv_up.resize(buf_size_y);

    if (down_rank != MPI_PROC_NULL) {
        int idx_sd = 0;
        for (int i = x_start_exchange; i < x_end_exchange; ++i) {
            for (int j = 0; j < halo_y; ++j) {
                for (int k = 0; k < M; ++k) {
                    scratch.send_down[idx_sd++] = u[i][j_start + j][k];
                }
            }
        }
    }

    if (up_rank != MPI_PROC_NULL) {
        int idx_su = 0;
        for (int i = x_start_exchange; i < x_end_exchange; ++i) {
            for (int j = 0; j < halo_y; ++j) {
                for (int k = 0; k < M; ++k) {
                    scratch.send_up[idx_su++] = u[i][j_end - halo_y + j][k];
                }
            }
        }
    }

    MPI_Sendrecv(scratch.send_down.data(), buf_size_y, MPI_DOUBLE, down_rank, 2,
                 scratch.recv_up.data(), buf_size_y, MPI_DOUBLE, up_rank, 2,
                 cart_comm, &status);
    MPI_Sendrecv(scratch.send_up.data(), buf_size_y, MPI_DOUBLE, up_rank, 3,
                 scratch.recv_down.data(), buf_size_y, MPI_DOUBLE, down_rank, 3,
                 cart_comm, &status);

    if (down_rank != MPI_PROC_NULL) {
        int idx_rd = 0;
        for (int i = x_start_exchange; i < x_end_exchange; ++i) {
            for (int j = 0; j < halo_y; ++j) {
                for (int k = 0; k < M; ++k) {
                    u[i][j_start - halo_y + j][k] = scratch.recv_down[idx_rd++];
                }
            }
        }
    }

    if (up_rank != MPI_PROC_NULL) {
        int idx_ru = 0;
        for (int i = x_start_exchange; i < x_end_exchange; ++i) {
            for (int j = 0; j < halo_y; ++j) {
                for (int k = 0; k < M; ++k) {
                    u[i][j_end + j][k] = scratch.recv_up[idx_ru++];
                }
            }
        }
    }
}

void exchange_halos_global(std::vector<std::vector<std::vector<double>>>& u,
                           MPI_Comm cart_comm, int fict_x, int fict_y,
                           int i_start, int i_end, int j_start, int j_end,
                           int left_rank, int right_rank, int down_rank, int up_rank) {
    exchange_halos_global_width(u, cart_comm, fict_x, fict_y,
                                i_start, i_end, j_start, j_end,
                                left_rank, right_rank, down_rank, up_rank);
}

void apply_physical_boundaries_local(std::vector<std::vector<std::vector<double>>>& u,
                                     int Nx, int Ny, int fict_x, int fict_y,
                                     int i_start, int i_end, int j_start, int j_end,
                                     int left_rank, int right_rank, int down_rank, int up_rank,
                                     int left_bc_code, int right_bc_code, int up_bc_code, int down_bc_code,
                                     double g) {
    if (left_rank == MPI_PROC_NULL) {
        for (int j = std::max(0, j_start - 1); j <= std::min(Ny - 1, j_end); ++j) {
            for (int k = 0; k < fict_x; ++k) {
                u[k][j] = boundary(u[fict_x][j], left_bc_code, g, 0);
            }
        }
    }
    if (right_rank == MPI_PROC_NULL) {
        for (int j = std::max(0, j_start - 1); j <= std::min(Ny - 1, j_end); ++j) {
            for (int k = 0; k < fict_x; ++k) {
                u[Nx - 1 - k][j] = boundary(u[Nx - fict_x - 1][j], right_bc_code, g, 0);
            }
        }
    }
    if (up_rank == MPI_PROC_NULL) {
        for (int i = std::max(0, i_start - 1); i <= std::min(Nx - 1, i_end); ++i) {
            for (int k = 0; k < fict_y; ++k) {
                u[i][Ny - 1 - k] = boundary(u[i][Ny - fict_y - 1], up_bc_code, g, 1);
            }
        }
    }
    if (down_rank == MPI_PROC_NULL) {
        for (int i = std::max(0, i_start - 1); i <= std::min(Nx - 1, i_end); ++i) {
            for (int k = 0; k < fict_y; ++k) {
                u[i][k] = boundary(u[i][fict_y], down_bc_code, g, 1);
            }
        }
    }
}

void gather_to_root(std::vector<std::vector<std::vector<double>>>& u,
                    int Nx, int Ny, int fict_x, int fict_y,
                    int i_start, int i_end, int j_start, int j_end,
                    int rank, int size, MPI_Comm cart_comm, int* dims) {
    const int local_nx = i_end - i_start;
    const int local_ny = j_end - j_start;
    const int buf_size = local_nx * local_ny * M;
    std::vector<double> send_buf(buf_size);

    int idx = 0;
    for (int i = i_start; i < i_end; ++i) {
        for (int j = j_start; j < j_end; ++j) {
            for (int k = 0; k < M; ++k) {
                send_buf[idx++] = u[i][j][k];
            }
        }
    }

    if (rank == 0) {
        for (int p = 1; p < size; ++p) {
            int p_coords[2];
            MPI_Cart_coords(cart_comm, p, 2, p_coords);

            int p_i_start = 0;
            int p_i_end = 0;
            int p_j_start = 0;
            int p_j_end = 0;
            get_subdomain_bounds(Nx, dims[0], p_coords[0], p_i_start, p_i_end);
            get_subdomain_bounds(Ny, dims[1], p_coords[1], p_j_start, p_j_end);
            p_i_start += fict_x;
            p_i_end += fict_x;
            p_j_start += fict_y;
            p_j_end += fict_y;

            const int p_buf_size = (p_i_end - p_i_start) * (p_j_end - p_j_start) * M;
            std::vector<double> recv_buf(p_buf_size);

            MPI_Recv(recv_buf.data(), p_buf_size, MPI_DOUBLE, p, 0, cart_comm, MPI_STATUS_IGNORE);

            int recv_idx = 0;
            for (int i = p_i_start; i < p_i_end; ++i) {
                for (int j = p_j_start; j < p_j_end; ++j) {
                    for (int k = 0; k < M; ++k) {
                        u[i][j][k] = recv_buf[recv_idx++];
                    }
                }
            }
        }
    } else {
        MPI_Send(send_buf.data(), buf_size, MPI_DOUBLE, 0, 0, cart_comm);
    }
}





