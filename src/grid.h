#ifndef CALC_AREA
#define CALC_AREA

#include <vector>
#include <iostream>
#include <fstream>
#include <tuple>

#include "point.h"
#include "solver.h"
#define RHO 0 // density
#define U 1 // speed x
#define V 2 // speed y
#define P 3 // pressure
// пїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ
#define r 0 // density
#define ru 1 // denity * speed x
#define rv 2 // denity * speed y
#define e 3 // energy

#define M 4


typedef std::tuple<int, int, int, int, double, double, double, double> grid_info;


enum class BoundaryType {
   None,
   Wall, 
   OpenFlow, 
   FixedValue,
};


struct Cell
{
public:
   BoundaryType type;

   double rho;
   double p;  
   Point u;

   Point center;
   std::vector<Point> vertices; 

   Cell(
      const BoundaryType _type, const double _rho, const double _p, 
      const Point _u, const Point _c, const std::vector<Point> _v) : 
      type(_type), rho(_rho), p(_p), u(_u), center(_c), vertices(_v) {};
   Cell();

   // double p_boundary;
   // Point u_boundary;
   // double rho_boundary;

   // void boundary(BoundaryType boundary, double p_old, Point u_old, double rho_old);
};


class Grid
{
friend class RimanSolver1D;

private:
   std::vector<Cell> grid;

   int fict_cells_x;
   int fict_cells_y;

   int cells_num_x;
   int cells_num_y;

   int total_cells_x;
   int total_cells_y;

   double left_boundary_x;
   double right_boundary_x;
   double left_boundary_y;
   double right_boundary_y;


public:
   Grid(
      int _fict_cells_x,
      int _fict_cells_y,
      int _cells_num_x,
      int _cells_num_y,
      double _left_boundary_x,
      double _right_boundary_x,
      double _left_boundary_y,
      double _right_boundary_y,
      const std::string &_left,
      const std::string &_right,
      const std::string &_up,
      const std::string &_down
   );

   Cell* get_cell(size_t i, size_t j);
   grid_info get_info();
   BoundaryType get_boundary_type(const std::string& type);

   void set_values(
      double rho_L, double rho_R, 
      double p_L, double p_R, 
      const Point u_L, const Point u_R
   );
   void WriteCSV(const std::string& filename);
};


#endif 