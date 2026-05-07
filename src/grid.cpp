#include "grid.h"

Cell::Cell()
{
    type = BoundaryType::None;
    rho = 0.0;
    p = 0.0;  
    u = Point();

    center = Point();
    vertices = std::vector<Point>(0, Point()); 
}

Grid::Grid(
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
    const std::string &_down)
{

    fict_cells_x = _fict_cells_x;
    fict_cells_y = _fict_cells_y;

    cells_num_x = _cells_num_x;
    cells_num_y = _cells_num_y;

    left_boundary_x  = _left_boundary_x;
    right_boundary_x = _right_boundary_x;
    left_boundary_y  = _left_boundary_y;
    right_boundary_y = _right_boundary_y;

    double hx = (right_boundary_x - left_boundary_x) / cells_num_x;
    int cells_num_fict_x = cells_num_x + 2 * fict_cells_x;
    total_cells_x = cells_num_fict_x;
    
    double hy = (right_boundary_y - left_boundary_y) / cells_num_y;
    int cells_num_fict_y = cells_num_y + 2 * fict_cells_y;
    total_cells_y = cells_num_fict_y;

    grid = std::vector<Cell>(cells_num_fict_x * cells_num_fict_y, Cell()); 
    
    for(int i = 0; i < cells_num_fict_x; i++){
        for(int j = 0; j < cells_num_fict_y; j++){
            Point left_up    = Point(left_boundary_x + i * hx, left_boundary_y + (j + 1) * hy);
            Point left_down  = Point(left_boundary_x + i * hx, left_boundary_y + j * hy);
            Point right_up   = Point(left_boundary_x + (i + 1) * hx, left_boundary_y + (j + 1) * hy);
            Point right_down = Point(left_boundary_x + (i + 1) * hx, left_boundary_y + j * hy);
            Cell* cell;
            cell = (this -> get_cell(i, j));
            cell -> vertices.push_back(left_up);
            cell -> vertices.push_back(left_down);
            cell -> vertices.push_back(right_up);
            cell -> vertices.push_back(right_down);            
            cell -> center = Point((left_up.x + right_up.x) / 2, (left_up.y + left_down.y) / 2);
            if (i < fict_cells_x) {
                cell -> type = this -> get_boundary_type(_left);
            } else if (j < fict_cells_y) {
                cell -> type = this -> get_boundary_type(_down);
            } else if (i >= (cells_num_fict_x - fict_cells_x)) {
                cell -> type = this -> get_boundary_type(_right);
            } else if (j >= (cells_num_fict_y - fict_cells_y)) {
                cell -> type = this -> get_boundary_type(_up);
            } 
        }
    }
}

Cell* Grid::get_cell(size_t i, size_t j)
{
    return &grid[i * total_cells_y + j];
}

void Grid::set_values(
    double rho_L, double rho_R, 
    double p_L, double p_R, 
    const Point u_L, const Point u_R
)
{
    for (int i = fict_cells_x; i < cells_num_x - fict_cells_x; ++i) {
        for (int j = fict_cells_y; j < cells_num_y - fict_cells_y; ++j) {
            Cell* cell;
            cell = (this -> get_cell(i, j));
            if (cell -> center.x < 0.5 * (left_boundary_x + right_boundary_x)) {
                cell -> p = p_L;
                cell -> u = u_L;
                cell -> rho = rho_L;
            } else {
                cell -> p = p_R;
                cell -> u = u_R;
                cell -> rho = rho_R;
            }
        }
    }
}

BoundaryType Grid::get_boundary_type(const std::string& type)
{
    if (type == "Wall") return BoundaryType::Wall;
    else if (type == "OpenFlow") return BoundaryType::OpenFlow;
    else if (type == "FixedValue") return BoundaryType::FixedValue;
    return BoundaryType::None;
}

void Grid::WriteCSV(const std::string& filename)
{
    std::ofstream fout(filename);
    fout << "x,y,p,vx,vy,rho\n";
    for (int i = fict_cells_x; i < cells_num_x - fict_cells_x; ++i) {
        for (int j = fict_cells_y; j < cells_num_y - fict_cells_y; ++j) {
            Cell* cell;
            cell = (this -> get_cell(i, j));
            fout << cell -> center.x << ","  << cell -> center.y << "," << cell -> p << "," << (cell -> u.x) << "," << cell -> u.y  << "," << cell -> rho << "\n";
        }
    }
    fout.close();
}

grid_info Grid::get_info()
{
    return std::make_tuple(fict_cells_x, fict_cells_y, cells_num_x, cells_num_y, left_boundary_x, right_boundary_x, left_boundary_y, right_boundary_y); 
}