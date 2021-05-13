#include "Boundary.hpp"
#include <cmath>
#include <iostream>
#include <cassert>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {
    for (const auto & cell : _cells) {
        std::vector<border_position> cell_borders = cell->borders();
        int i = cell->i();
        int j = cell->j();

        if (cell_borders.size() == 1) {

            for (const auto & elem : cell_borders) {
                if (static_cast<int>(elem) == 0) {      // border is in TOP position (fluid cell above ghost cell)
                    field.u(i,j) = -field.u(i,j+1);     
                    field.v(i,j) = 0;
                    field.p(i,j) = field.p(i,j+1);
                    field.f(i,j) = field.u(i,j);
                    field.g(i,j) = field.v(i,j);
                }
                else if (static_cast<int>(elem) == 1) { // border is in BOTTOM position (fluid cell under ghost cell)
                    field.u(i,j)   = - field.u(i,j-1);
                    field.v(i,j-1) = 0;
                    field.p(i,j)   = field.p(i,j-1);
                    field.f(i,j)   = field.u(i,j);
                    field.g(i,j-1)   = field.v(i,j-1);
                }
                else if (static_cast<int>(elem) == 2) { // border is in LEFT position (fluid cell is left from ghost cell)
                    field.u(i-1,j) = 0;
                    field.v(i,j)   = - field.v(i-1,j);
                    field.p(i,j)   = field.p(i-1,j);
                    field.f(i-1,j)   = field.u(i-1,j);
                    field.g(i,j)   = field.v(i,j);
                }
                else if (static_cast<int>(elem) == 3) { // border is in RIGHT position (fluid cell is right from ghost cell)
                    field.u(i,j) = 0;
                    field.v(i,j) = - field.v(i+1,j);
                    field.p(i,j) = field.p(i+1,j);
                    field.f(i,j) = field.u(i,j);
                    field.g(i,j) = field.v(i,j);
                }
            }
        }

        else if (cell_borders.size() == 2) {
            
            if ( (cell_borders[0] == border_position::TOP && cell_borders[1] == border_position::RIGHT)
                || (cell_borders[0] == border_position::RIGHT && cell_borders[1] == border_position::TOP) ) {   // borders in TOP and RIGHT position

                field.u(i,j) = 0;
                field.v(i,j) = 0;
                field.u(i-1,j) = - field.u(i-1,j+1);
                field.v(i,j-1) = - field.v(i+1,j-1);
                field.p(i,j) = 0.5*(field.p(i+1,j) + field.p(i,j+1));
                field.f(i,j) = field.u(i,j);
                field.g(i,j) = field.v(i,j);
            }
            else if ( (cell_borders[0] == border_position::TOP && cell_borders[1] == border_position::LEFT) 
                || (cell_borders[0] == border_position::LEFT && cell_borders[1] == border_position::TOP) ) {    // borders in TOP and LEFT position
                
                field.u(i,j) = - field.u(i,j+1);
                field.v(i,j) = 0;
                field.u(i-1,j) = 0;
                field.v(i,j-1) = - field.v(i-1,j-1);
                field.p(i,j) = 0.5*(field.p(i-1,j) + field.p(i,j+1));
                field.f(i,j) = field.u(i,j);
                field.g(i,j) = field.v(i,j);
            }
            else if ( (cell_borders[0] == border_position::BOTTOM && cell_borders[1] == border_position::RIGHT) 
                || (cell_borders[0] == border_position::RIGHT && cell_borders[1] == border_position::BOTTOM) ) { // borders in BOTTOM and RIGHT position
                
                field.u(i,j) = 0;
                field.v(i,j) = - field.v(i+1,j);
                field.u(i-1,j) = - field.u(i-1,j-1);
                field.v(i,j-1) = 0;
                field.p(i,j) = 0.5*(field.p(i+1,j) + field.p(i,j-1));
                field.f(i,j) = field.u(i,j);
                field.g(i,j) = field.v(i,j);
                }
            else if ( (cell_borders[0] == border_position::BOTTOM && cell_borders[1] == border_position::LEFT) 
                || (cell_borders[0] == border_position::LEFT && cell_borders[1] == border_position::BOTTOM) ) { // borders in BOTTOM and LEFT position

                field.u(i,j) = - field.u(i,j-1);
                field.v(i,j) = - field.v(i-1,j);
                field.u(i-1,j) = 0;
                field.v(i,j-1) = 0;
                field.p(i,j) = 0.5*(field.p(i-1,j) + field.p(i,j-1));
                field.f(i,j) = field.u(i,j);
                field.g(i,j) = field.v(i,j);
                }
        }
    }
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {
     for (const auto & cell : _cells) {
        std::vector<border_position> cell_borders = cell->borders();
        int i = cell->i();
        int j = cell->j();

        if (cell_borders.size() == 1) {

            for (const auto & elem : cell_borders) {
                if (static_cast<int>(elem) == 0) {      // border is in TOP position (fluid cell above ghost cell)
                    field.u(i,j) = 2 * _wall_velocity[LidDrivenCavity::moving_wall_id] - field.u(i,j+1);     
                    field.v(i,j) = 0;
                    field.p(i,j) = field.p(i,j+1);
                    field.f(i,j) = field.u(i,j);
                    field.g(i,j) = field.v(i,j);
                }
                else if (static_cast<int>(elem) == 1) { // border is in BOTTOM position (fluid cell under ghost cell)
                    field.u(i,j) = 2 * _wall_velocity[LidDrivenCavity::moving_wall_id] - field.u(i,j-1);
                    field.v(i,j-1) = 0;
                    field.p(i,j) = field.p(i,j-1);
                    field.f(i,j) = field.u(i,j);
                    field.g(i,j-1) = field.v(i,j-1);
                }
                else if (static_cast<int>(elem) == 2) { // border is in LEFT position (fluid cell left from ghost cell)
                    field.u(i-1,j) = 0;
                    field.v(i,j) = 2 * _wall_velocity[LidDrivenCavity::moving_wall_id] - field.v(i-1,j);
                    field.p(i,j) = field.p(i-1,j);
                    field.f(i-1,j) = field.u(i-1,j);
                    field.g(i,j) = field.v(i,j);
                }
                else if (static_cast<int>(elem) == 3) { // border is in RIGHT position (fluid cell right from ghost cell)
                    field.u(i,j) = 0;
                    field.v(i,j) = 2 * _wall_velocity[LidDrivenCavity::moving_wall_id] - field.v(i+1,j);
                    field.p(i,j) = field.p(i+1,j);
                    field.f(i,j) = field.u(i,j);
                    field.g(i,j) = field.v(i,j);
                }
            }
        }

        else if (cell_borders.size() == 2) {

        }
    }

}
