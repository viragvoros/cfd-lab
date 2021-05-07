#include "Boundary.hpp"
#include <cmath>
#include <iostream>
#include <cassert>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field, int imax, int jmax) {
    for (const auto & cell : _cells) {
        std::vector<border_position> cell_borders = cell->borders();
        int i = cell->i();
        int j = cell->j();

        if (not cell_borders.empty()) {

            for (const auto & elem : cell_borders) {
                if (elem == border_position::TOP) {      // border is in TOP position (fluid cell above ghost cell)
                    assert(j == 0);
                    field.u(i,j) = -field.u(i,j+1);     
                    field.v(i,j) = 0;
                    field.p(i,j) = field.p(i,j+1);
                    field.f(i,j) = field.u(i,j);
                    field.g(i,j) = field.v(i,j);
                }
                else if (elem == border_position::BOTTOM) { // border is in BOTTOM position (fluid cell under ghost cell)
                    assert(j == jmax + 1);
                    field.u(i,j)   = - field.u(i,j-1);
                    field.v(i,j-1) = 0;
                    field.p(i,j)   = field.p(i,j-1);
                    field.f(i,j)   = field.u(i,j);
                    field.g(i,j-1)   = field.v(i,j-1);
                }
                else if (elem == border_position::LEFT) { // border is in LEFT position (fluid cell is left from ghost cell)
                    assert( i == imax + 1);
                    field.u(i-1,j) = 0;
                    field.v(i,j)   = - field.v(i-1,j);
                    field.p(i,j)   = field.p(i-1,j);
                    field.f(i-1,j)   = field.u(i-1,j);
                    field.g(i,j)   = field.v(i,j);
                }
                else if (elem == border_position::RIGHT) { // border is in RIGHT position (fluid cell is right from ghost cell)
                    assert( i == 0);
                    field.u(i,j) = 0;
                    field.v(i,j) = - field.v(i+1,j);
                    field.p(i,j) = field.p(i+1,j);
                    field.f(i,j) = field.u(i,j);
                    field.g(i,j) = field.v(i,j);
                }
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

void MovingWallBoundary::apply(Fields &field, int imax, int jmax) {
     for (const auto & cell : _cells) {
        std::vector<border_position> cell_borders = cell->borders();
        int i = cell->i();
        int j = cell->j();

        if (not cell_borders.empty()) {

            for (const auto & elem : cell_borders) {
                if (elem == border_position::TOP) {      // border is in TOP position (fluid cell above ghost cell)
                    assert( j == 0);
                    field.u(i,j) = 2 * _wall_velocity[LidDrivenCavity::moving_wall_id] - field.u(i,j+1);     
                    field.v(i,j) = 0;
                    field.p(i,j) = field.p(i,j+1);
                    field.f(i,j) = field.u(i,j);
                    field.g(i,j) = field.v(i,j);
                }
                else if (elem == border_position::BOTTOM) { // border is in BOTTOM position (fluid cell under ghost cell)
                    assert( j == jmax + 1);
                    field.u(i,j) = 2 * _wall_velocity[LidDrivenCavity::moving_wall_id] - field.u(i,j-1);
                    field.v(i,j-1) = 0;
                    field.p(i,j) = field.p(i,j-1);
                    field.f(i,j) = field.u(i,j);
                    field.g(i,j-1) = field.v(i,j-1);
                }
                else if (elem == border_position::LEFT) { // border is in LEFT position (fluid cell left from ghost cell)
                    assert( i == imax + 1);
                    field.u(i-1,j) = 0;
                    field.v(i,j) = 2 * _wall_velocity[LidDrivenCavity::moving_wall_id] - field.v(i-1,j);
                    field.p(i,j) = field.p(i-1,j);
                    field.f(i-1,j) = field.u(i-1,j);
                    field.g(i,j) = field.v(i,j);
                }
                else if (elem == border_position::RIGHT) { // border is in RIGHT position (fluid cell right from ghost cell)
                    assert( i == 0);
                    field.u(i,j) = 0;
                    field.v(i,j) = 2 * _wall_velocity[LidDrivenCavity::moving_wall_id] - field.v(i+1,j);
                    field.p(i,j) = field.p(i+1,j);
                    field.f(i,j) = field.u(i,j);
                    field.g(i,j) = field.v(i,j);
                }
            }
        }
    }

}
