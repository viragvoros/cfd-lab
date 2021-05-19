#include "Boundary.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, double wall_temperature) : _cells(cells) {
    _wall_temperature.insert(std::pair(boundary_ids::fixed_wall_cell_3_id, wall_temperature));
}

void FixedWallBoundary::apply(Fields &field) {
    for (const auto &cell : _cells) {
        std::vector<border_position> cell_borders = cell->borders();
        int i = cell->i();
        int j = cell->j();

        if (cell_borders.empty()){
            field.u(i, j) = 0;
                field.v(i, j) = 0;
                field.p(i, j) = 0;
                field.f(i, j) = 0;
                field.g(i, j) = 0;
        }

        if (cell_borders.size() == 1) { // One  neighbour

            if (cell_borders[0] == border_position::TOP) { // border is in TOP position (fluid cell above ghost cell)
                field.u(i, j) = -field.u(i, j + 1);
                field.v(i, j) = 0;
                field.p(i, j) = field.p(i, j + 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            } else if (cell_borders[0] ==
                       border_position::BOTTOM) { // border is in BOTTOM position (fluid cell under ghost cell)
                field.u(i, j) = -field.u(i, j - 1);
                field.v(i, j - 1) = 0;
                field.p(i, j) = field.p(i, j - 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j - 1) = field.v(i, j - 1);
            } else if (cell_borders[0] ==
                       border_position::LEFT) { // border is in LEFT position (fluid cell is left from ghost cell)
                field.u(i - 1, j) = 0;
                field.v(i, j) = -field.v(i - 1, j);
                field.p(i, j) = field.p(i - 1, j);
                field.f(i - 1, j) = field.u(i - 1, j);
                field.g(i, j) = field.v(i, j);
            } else if (cell_borders[0] ==
                       border_position::RIGHT) { // border is in RIGHT position (fluid cell is right from ghost cell)
                field.u(i, j) = 0;
                field.v(i, j) = -field.v(i + 1, j);
                field.p(i, j) = field.p(i + 1, j);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            }
        }

        else if (cell_borders.size() == 2) {

            

            if ((cell_borders[0] == border_position::TOP && cell_borders[1] == border_position::RIGHT) ||
                (cell_borders[0] == border_position::RIGHT &&
                 cell_borders[1] == border_position::TOP)) { // borders in TOP and RIGHT position
                 //std::cout << "Top Right" << std::endl;

                field.u(i, j) = 0;
                field.v(i, j) = 0;
                field.u(i - 1, j) = -field.u(i - 1, j + 1);
                field.v(i, j - 1) = -field.v(i + 1, j - 1);
                field.p(i, j) = 0.5 * (field.p(i + 1, j) + field.p(i, j + 1));
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            } else if ((cell_borders[0] == border_position::TOP && cell_borders[1] == border_position::LEFT) ||
                       (cell_borders[0] == border_position::LEFT &&
                        cell_borders[1] == border_position::TOP)) { // borders in TOP and LEFT position
                        //std::cout << "Top Left" << std::endl;

                field.u(i, j) = -field.u(i, j + 1);
                field.v(i, j) = 0;
                field.u(i - 1, j) = 0;
                field.v(i, j - 1) = -field.v(i - 1, j - 1);
                field.p(i, j) = 0.5 * (field.p(i - 1, j) + field.p(i, j + 1));
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            } else if ((cell_borders[0] == border_position::BOTTOM && cell_borders[1] == border_position::RIGHT) ||
                       (cell_borders[0] == border_position::RIGHT &&
                        cell_borders[1] == border_position::BOTTOM)) { // borders in BOTTOM and RIGHT position
                        //std::cout << "Bottom Right" << std::endl;

                field.u(i, j) = 0;
                field.v(i, j) = -field.v(i + 1, j);
                field.u(i - 1, j) = -field.u(i - 1, j - 1);
                field.v(i, j - 1) = 0;
                field.p(i, j) = 0.5 * (field.p(i + 1, j) + field.p(i, j - 1));
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            } else if ((cell_borders[0] == border_position::BOTTOM && cell_borders[1] == border_position::LEFT) ||
                       (cell_borders[0] == border_position::LEFT &&
                        cell_borders[1] == border_position::BOTTOM)) { // borders in BOTTOM and LEFT position
                        //std::cout << "Bottom Left" << std::endl;

                field.u(i, j) = -field.u(i, j - 1);
                field.v(i, j) = -field.v(i - 1, j);
                field.u(i - 1, j) = 0;
                field.v(i, j - 1) = 0;
                field.p(i, j) = 0.5 * (field.p(i - 1, j) + field.p(i, j - 1));
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            }
        }
    }


    // std::cout << "---------------------- u field ------------------------------" << std::endl;
    //  for (int jx = 0; jx < 22; jx++ ){
    //              for (int ix = 0; ix < 102; ix++) {
    //                  std::cout << field.u(ix, jx) << " " ;
    //                 }
    //             std::cout << "\n";
    //          }

    // std::cout << "---------------------- v field ------------------------------" << std::endl;
    //  for (int jx = 0; jx < 22; jx++ ){
    //              for (int ix = 0; ix < 102; ix++) {
    //                  std::cout << field.v(ix, jx) << " " ;
    //                 }
    //             std::cout << "\n";
    //          }
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {
    for (const auto &cell : _cells) {
        std::vector<border_position> cell_borders = cell->borders();
        int i = cell->i();
        int j = cell->j();
        int id = cell->wall_id();

        if (cell_borders.size() == 1) {

            // TODO: replace moving_wall_id with cell->id()
            if (cell_borders[0] == border_position::TOP) { // border is in TOP position (fluid cell above ghost cell)
                field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j + 1);
                field.v(i, j) = 0;
                field.p(i, j) = field.p(i, j + 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            } else if (cell_borders[0] ==
                       border_position::BOTTOM) { // border is in BOTTOM position (fluid cell under ghost cell)
                field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j - 1);
                field.v(i, j - 1) = 0;
                field.p(i, j) = field.p(i, j - 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j - 1) = field.v(i, j - 1);
            } else if (cell_borders[0] ==
                       border_position::LEFT) { // border is in LEFT position (fluid cell left from ghost cell)
                field.u(i - 1, j) = 0;
                field.v(i, j) = 2 * _wall_velocity[id] - field.v(i - 1, j);
                field.p(i, j) = field.p(i - 1, j);
                field.f(i - 1, j) = field.u(i - 1, j);
                field.g(i, j) = field.v(i, j);
            } else if (cell_borders[0] ==
                       border_position::RIGHT) { // border is in RIGHT position (fluid cell right from ghost cell)
                field.u(i, j) = 0;
                field.v(i, j) = 2 * _wall_velocity[id] - field.v(i + 1, j);
                field.p(i, j) = field.p(i + 1, j);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            }
        }

        /* NOT NEEDED YET, TODO LATER
        else if (cell_borders.size() == 2) {
            // do same as above
        }
        */
    }
}

InFlowBoundary::InFlowBoundary(std::vector<Cell *> cells, double inflow_velocity) : _cells(cells) {
    _inflow_velocity.insert(std::pair(boundary_ids::inflow_cell_id, inflow_velocity));
}

InFlowBoundary::InFlowBoundary(std::vector<Cell *> cells, double inflow_velocity, double inflow_temperature)
    : _cells(cells) {
    _inflow_velocity.insert(std::pair(boundary_ids::inflow_cell_id, inflow_velocity));
    _inflow_temperature.insert(std::pair(boundary_ids::inflow_cell_id, inflow_temperature));
}

InFlowBoundary::InFlowBoundary(std::vector<Cell *> cells, std::map<int, double> inflow_velocity,
                               std::map<int, double> inflow_temperature)
    : _cells(cells), _inflow_velocity(inflow_velocity), _inflow_temperature(inflow_temperature) {}

void InFlowBoundary::apply(Fields &field) {
    for (const auto cell : _cells) {
        std::vector<border_position> cell_borders = cell->borders();
        int i = cell->i();
        int j = cell->j();
        int id = cell->wall_id();

        if (cell_borders.size() == 1) {

            if (cell_borders[0] == border_position::TOP) { // border is in TOP position
                field.u(i, j) = -field.u(i, j + 1);
                field.v(i, j) = 1;
                field.p(i, j) = field.p(i, j + 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            } else if (cell_borders[0] == border_position::BOTTOM) { // border is in BOTTOM position
                field.u(i, j) = -field.u(i, j - 1);
                field.v(i, j - 1) = -1;
                field.p(i, j) = field.p(i, j - 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            } else if (cell_borders[0] == border_position::LEFT) { // border is in LEFT position
                field.u(i - 1, j) = -1;
                field.v(i, j) = -field.v(i - 1, j);
                field.p(i, j) = field.p(i - 1, j);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            } else if (cell_borders[0] == border_position::RIGHT) { // border is in RIGHT position
                field.u(i, j) = 1;
                field.v(i, j) = -field.v(i + 1, j);
                field.p(i, j) = field.p(i + 1, j);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            }
        }

        /* NOT NEEDED YET, TODO LATER

        else if (cell_borders.size() == 2) {

            if ( (cell_borders[0] == border_position::TOP && cell_borders[1] == border_position::RIGHT)
                || (cell_borders[0] == border_position::RIGHT && cell_borders[1] == border_position::TOP) ) {   //
        borders in TOP and RIGHT position


            }
        }
        */
    }
}

OutFlowBoundary::OutFlowBoundary(std::vector<Cell *> cells, double outflow_velocity) : _cells(cells) {
    _outflow_velocity.insert(std::pair(boundary_ids::outflow_cell_id, outflow_velocity));
}

OutFlowBoundary::OutFlowBoundary(std::vector<Cell *> cells, double outflow_velocity, double outflow_temperature)
    : _cells(cells) {
    _outflow_velocity.insert(std::pair(boundary_ids::outflow_cell_id, outflow_velocity));
    _outflow_temperature.insert(std::pair(boundary_ids::outflow_cell_id, outflow_temperature));
}

OutFlowBoundary::OutFlowBoundary(std::vector<Cell *> cells, std::map<int, double> outflow_velocity,
                                 std::map<int, double> outflow_temperature)
    : _cells(cells), _outflow_velocity(outflow_velocity), _outflow_temperature(outflow_temperature) {}

void OutFlowBoundary::apply(Fields &field) {
    for (const auto cell : _cells) {
        std::vector<border_position> cell_borders = cell->borders();
        int i = cell->i();
        int j = cell->j();
        int id = cell->wall_id();

        if (cell_borders.size() == 1) {

            if (cell_borders[0] == border_position::TOP) { // border is in TOP position
                field.u(i, j) = field.u(i, j + 1);
                field.v(i, j) = field.v(i,j+1);
                field.p(i, j) = 0;  // TODO remove hard coded value
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            } else if (cell_borders[0] == border_position::BOTTOM) { // border is in BOTTOM position
                field.u(i, j) = field.u(i, j - 1);
                field.v(i, j) = field.v(i,j-1);
                field.p(i, j) = 0; // TODO remove hard coded value
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            } else if (cell_borders[0] == border_position::LEFT) { // border is in LEFT position
                //std::cout << "Applying Outflow" << std::endl;
                field.u(i, j) = field.u(i - 1, j);
                field.v(i, j) = field.v(i - 1, j);
                field.p(i, j) = 0.0;
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            } else if (cell_borders[0] == border_position::RIGHT) { // border is in RIGHT position
                field.u(i, j) = field.u(i+1,j);
                field.v(i, j) = field.v(i + 1, j);
                field.p(i, j) = 0;
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
            }
        }

            //  for (int jx = 0; jx < 22; jx++ ){
            //      for (int ix = 0; ix < 102; ix++) {
            //          std::cout << field.u(ix, jx) << " " ;
            //         }
            //     std::cout << "\n";
            //  }



        /* NOT NEEDED YET, TODO LATER

        else if (cell_borders.size() == 2) {

            if ( (cell_borders[0] == border_position::TOP && cell_borders[1] == border_position::RIGHT)
                || (cell_borders[0] == border_position::RIGHT && cell_borders[1] == border_position::TOP) ) {   //
        borders in TOP and RIGHT position


            }
        }
        */
    }
}