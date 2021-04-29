#include "Boundary.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {
    for (auto & cell : _cells) {
        std::vector<border_position> cell_borders = cell->borders();
        int i = cell->i();
        int j = cell->j();

        if (not cell_borders.empty()) {

            for (const auto & elem : cell_borders) {
                if (static_cast<int>(elem) == 0) {      // border is in TOP position
                    field.u(i,j) = -field.u(i,j+1);     
                    field.v(i,j) = 0;
                    field.p(i,j) = field.p(i,j+1);
                    field.f(i,j) = field.u(i,j);
                    field.g(i,j) = field.v(i,j);
                }
                else if (static_cast<int>(elem) == 1) { // border is in BOTTOM position
                    field.u(i,j)   = - field.u(i,j-1);
                    field.v(i,j-1) = 0;
                    field.p(i,j)   = field.p(i,j-1);
                    field.f(i,j)   = field.u(i,j);
                    field.g(i,j-1)   = field.v(i,j-1);
                }
                else if (static_cast<int>(elem) == 2) { // border is in LEFT position
                    field.u(i-1,j) = 0;
                    field.v(i,j)   = - field.v(i-1,j);
                    field.p(i,j)   = field.p(i-1,j);
                    field.f(i-1,j)   = field.u(i-1,j);
                    field.g(i,j)   = field.v(i,j);
                }
                else if (static_cast<int>(elem) == 3) { // border is in RIGHT position
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

void MovingWallBoundary::apply(Fields &field) {

     for (auto & cell : _cells) {
        std::vector<border_position> cell_borders = cell->borders();
        int i = cell->i();
        int j = cell->j();

        if (not cell_borders.empty()) {

            for (auto & elem : cell_borders) {
                if (static_cast<int>(elem) == 0) {      // border is in TOP position
                    field.u(i,j) = 2 * 1 - field.u(i,j+1);     
                    field.v(i,j) = 0;
                    field.p(i,j) = field.p(i,j+1);
                    field.f(i,j) = field.u(i,j);
                    field.g(i,j) = field.v(i,j);
                }
                else if (static_cast<int>(elem) == 1) { // border is in BOTTOM position
                    field.u(i,j) = 2 * 1 - field.u(i,j-1);
                    field.v(i,j-1) = 0;
                    field.p(i,j) = field.p(i,j-1);
                    field.f(i,j) = field.u(i,j);
                    field.g(i,j-1) = field.v(i,j-1);
                }
                else if (static_cast<int>(elem) == 2) { // border is in LEFT position
                    field.u(i-1,j) = 0;
                    field.v(i,j) = 2 * 0 - field.v(i-1,j);
                    field.p(i,j) = field.p(i-1,j);
                    field.f(i-1,j) = field.u(i-1,j);
                    field.g(i,j) = field.v(i,j);
                }
                else if (static_cast<int>(elem) == 3) { // border is in RIGHT position
                    field.u(i,j) = 0;
                    field.v(i,j) = 2 * 0 - field.v(i+1,j);
                    field.p(i,j) = field.p(i+1,j);
                    field.f(i,j) = field.u(i,j);
                    field.g(i,j) = field.v(i,j);
                }
            }
        }
    }

}
