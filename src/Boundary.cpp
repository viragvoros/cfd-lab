#include "Boundary.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature, double kappa,
                                     std::map<int, double> wall_heatflux)
    : _cells(cells), _wall_temperature(wall_temperature), _kappa(kappa), _wall_heatflux(wall_heatflux) {}

// Constructor only used for LidDrivenCavity, where no geometry is available
FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, double wall_temperature) : _cells(cells) {
    _wall_temperature.insert(std::pair(boundary_ids::fixed_wall_cell_3_id, wall_temperature));
}

void FixedWallBoundary::apply(Fields &field) {
    for (const auto &cell : _cells) {
        std::vector<border_position> cell_borders = cell->borders();
        int i = cell->i();
        int j = cell->j();
        int id = cell->cell_id();

        if (cell_borders.empty()) { // Set all values of cells without a fluid neighbour to 0
            field.u(i, j) = 0;
            field.v(i, j) = 0;
            field.p(i, j) = 0;
            field.f(i, j) = 0;
            field.g(i, j) = 0;
        }

        if (cell_borders.size() == 1) { // One fluid neighbour

            if (cell_borders[0] == border_position::TOP) { // Border is in TOP position (fluid cell above ghost cell)
                field.u(i, j) = -field.u(i, j + 1);
                field.v(i, j) = 0;
                field.p(i, j) = field.p(i, j + 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);

                field.ca(i, j) = field.ca(i, j + 1);
                field.cb(i, j) = field.cb(i, j + 1);
                field.cc(i, j) = field.cc(i, j + 1);

                if (!_wall_temperature.empty() && _wall_temperature[id] != -1) {
                    field.t(i, j) = 2 * _wall_temperature[id] - field.t(i, j + 1);
                } else if (!_wall_temperature.empty() && _wall_temperature[id] == -1) {
                    field.t(i, j) = field.t(i, j + 1) + _wall_heatflux[id] / _kappa;
                }

            } else if (cell_borders[0] ==
                       border_position::BOTTOM) { // Border is in BOTTOM position (fluid cell under ghost cell)
                field.u(i, j) = -field.u(i, j - 1);
                field.v(i, j) = 0;
                field.v(i, j - 1) = 0;
                field.p(i, j) = field.p(i, j - 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j - 1) = field.v(i, j - 1);

                field.ca(i, j) = field.ca(i, j - 1);
                field.cb(i, j) = field.cb(i, j - 1);
                field.cc(i, j) = field.cc(i, j - 1);

                if (!_wall_temperature.empty() && _wall_temperature[id] != -1) {
                    field.t(i, j) = 2 * _wall_temperature[id] - field.t(i, j - 1);
                } else if (!_wall_temperature.empty() && _wall_temperature[id] == -1) {
                    field.t(i, j) = field.t(i, j - 1) + _wall_heatflux[id] / _kappa;
                }
            } else if (cell_borders[0] ==
                       border_position::LEFT) { // Border is in LEFT position (fluid cell is left from ghost cell)
                field.u(i - 1, j) = 0;
                field.u(i, j) = 0;
                field.v(i, j) = -field.v(i - 1, j);
                field.p(i, j) = field.p(i - 1, j);
                field.f(i - 1, j) = field.u(i - 1, j);
                field.g(i, j) = field.v(i, j);

                field.ca(i, j) = field.ca(i - 1, j);
                field.cb(i, j) = field.cb(i - 1, j);
                field.cc(i, j) = field.cc(i - 1, j);

                if (!_wall_temperature.empty() && _wall_temperature[id] != -1) {
                    field.t(i, j) = 2 * _wall_temperature[id] - field.t(i - 1, j);
                } else if (!_wall_temperature.empty() && _wall_temperature[id] == -1) {
                    field.t(i, j) = field.t(i - 1, j) + _wall_heatflux[id] / _kappa;
                }

            } else if (cell_borders[0] ==
                       border_position::RIGHT) { // Border is in RIGHT position (fluid cell is right from ghost cell)
                field.u(i, j) = 0;
                field.v(i, j) = -field.v(i + 1, j);
                field.p(i, j) = field.p(i + 1, j);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);

                field.ca(i, j) = field.ca(i + 1, j);
                field.cb(i, j) = field.cb(i + 1, j);
                field.cc(i, j) = field.cc(i + 1, j);

                if (!_wall_temperature.empty() && _wall_temperature[id] != -1) {
                    field.t(i, j) = 2 * _wall_temperature[id] - field.t(i + 1, j);
                } else if (!_wall_temperature.empty() && _wall_temperature[id] == -1) {
                    field.t(i, j) = field.t(i + 1, j) + _wall_heatflux[id] / _kappa;
                }
            }
        }

        else if (cell_borders.size() == 2) {

            if ((cell_borders[0] == border_position::TOP && cell_borders[1] == border_position::RIGHT) ||
                (cell_borders[0] == border_position::RIGHT &&
                 cell_borders[1] == border_position::TOP)) { // Borders in TOP and RIGHT position

                field.u(i, j) = 0;
                field.v(i, j) = 0;
                field.u(i - 1, j) = -field.u(i - 1, j + 1);
                field.v(i, j - 1) = -field.v(i + 1, j - 1);
                field.p(i, j) = 0.5 * (field.p(i + 1, j) + field.p(i, j + 1));
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);

                field.ca(i, j) = 0.5 * (field.ca(i + 1, j) + field.ca(i, j + 1));
                field.cb(i, j) = 0.5 * (field.cb(i + 1, j) + field.cb(i, j + 1));
                field.cc(i, j) = 0.5 * (field.cc(i + 1, j) + field.cc(i, j + 1));

                if (!_wall_temperature.empty() && _wall_temperature[id] != -1) {
                    field.t(i, j) = 2 * _wall_temperature[id] - 0.5 * (field.t(i, j + 1) + field.t(i + 1, j));
                } else if (!_wall_temperature.empty() && _wall_temperature[id] == -1) {
                    field.t(i, j) = 0.5 * (field.t(i + 1, j) + field.t(i, j + 1)) + _wall_heatflux[id] / _kappa;
                }

            } else if ((cell_borders[0] == border_position::TOP && cell_borders[1] == border_position::LEFT) ||
                       (cell_borders[0] == border_position::LEFT &&
                        cell_borders[1] == border_position::TOP)) { // Borders in TOP and LEFT position

                field.u(i, j) = -field.u(i, j + 1);
                field.v(i, j) = 0;
                field.u(i - 1, j) = 0;
                field.v(i, j - 1) = -field.v(i - 1, j - 1);
                field.p(i, j) = 0.5 * (field.p(i - 1, j) + field.p(i, j + 1));
                field.f(i - 1, j) = field.u(i - 1, j);
                field.g(i, j) = field.v(i, j);

                field.ca(i, j) = 0.5 * (field.ca(i - 1, j) + field.ca(i, j + 1));
                field.cb(i, j) = 0.5 * (field.cb(i - 1, j) + field.cb(i, j + 1));
                field.cc(i, j) = 0.5 * (field.cc(i - 1, j) + field.cc(i, j + 1));

                if (!_wall_temperature.empty() && _wall_temperature[id] != -1) {
                    field.t(i, j) = 2 * _wall_temperature[id] - 0.5 * (field.t(i, j + 1) + field.t(i - 1, j));
                } else if (!_wall_temperature.empty() && _wall_temperature[id] == -1) {
                    field.t(i, j) = 0.5 * (field.t(i - 1, j) + field.t(i, j + 1)) + _wall_heatflux[id] / _kappa;
                }

            } else if ((cell_borders[0] == border_position::BOTTOM && cell_borders[1] == border_position::RIGHT) ||
                       (cell_borders[0] == border_position::RIGHT &&
                        cell_borders[1] == border_position::BOTTOM)) { // Borders in BOTTOM and RIGHT position

                field.u(i, j) = 0;
                field.v(i, j) = -field.v(i + 1, j);
                field.u(i - 1, j) = -field.u(i - 1, j - 1);
                field.v(i, j - 1) = 0;
                field.p(i, j) = 0.5 * (field.p(i + 1, j) + field.p(i, j - 1));
                field.f(i, j) = field.u(i, j);
                field.g(i, j - 1) = field.v(i, j - 1);

                field.ca(i, j) = 0.5 * (field.ca(i + 1, j) + field.ca(i, j - 1));
                field.cb(i, j) = 0.5 * (field.cb(i + 1, j) + field.cb(i, j - 1));
                field.cc(i, j) = 0.5 * (field.cc(i + 1, j) + field.cc(i, j - 1));

                if (!_wall_temperature.empty() && _wall_temperature[id] != -1) {
                    field.t(i, j) = 2 * _wall_temperature[id] - 0.5 * (field.t(i, j - 1) + field.t(i + 1, j));
                } else if (!_wall_temperature.empty() && _wall_temperature[id] == -1) {
                    field.t(i, j) = 0.5 * (field.t(i + 1, j) + field.t(i, j - 1)) + _wall_heatflux[id] / _kappa;
                }

            } else if ((cell_borders[0] == border_position::BOTTOM && cell_borders[1] == border_position::LEFT) ||
                       (cell_borders[0] == border_position::LEFT &&
                        cell_borders[1] == border_position::BOTTOM)) { // Borders in BOTTOM and LEFT position

                field.u(i, j) = -field.u(i, j - 1);
                field.v(i, j) = -field.v(i - 1, j);
                field.u(i - 1, j) = 0;
                field.v(i, j - 1) = 0;
                field.p(i, j) = 0.5 * (field.p(i - 1, j) + field.p(i, j - 1));
                field.f(i - 1, j) = field.u(i - 1, j);
                field.g(i, j - 1) = field.v(i, j - 1);

                field.ca(i, j) = 0.5 * (field.ca(i + 1, j) + field.ca(i, j - 1));
                field.cb(i, j) = 0.5 * (field.cb(i + 1, j) + field.cb(i, j - 1));
                field.cc(i, j) = 0.5 * (field.cc(i + 1, j) + field.cc(i, j - 1));

                if (!_wall_temperature.empty() && _wall_temperature[id] != -1) {
                    field.t(i, j) = 2 * _wall_temperature[id] - 0.5 * (field.t(i, j - 1) + field.t(i - 1, j));
                } else if (!_wall_temperature.empty() && _wall_temperature[id] == -1) {
                    field.t(i, j) = 0.5 * (field.t(i - 1, j) + field.t(i, j - 1)) + _wall_heatflux[id] / _kappa;
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
    for (const auto &cell : _cells) {
        std::vector<border_position> cell_borders = cell->borders();
        int i = cell->i();
        int j = cell->j();
        int id = cell->cell_id();

        if (cell_borders.size() == 1) {

            if (cell_borders[0] == border_position::TOP) { // Border is in TOP position (fluid cell above ghost cell)
                field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j + 1);
                field.v(i, j) = 0;
                field.p(i, j) = field.p(i, j + 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);

                field.ca(i, j) = field.ca(i, j + 1);
                field.cb(i, j) = field.cb(i, j + 1);
                field.cc(i, j) = field.cc(i, j + 1);

                if (!_wall_temperature.empty()) {
                    field.t(i, j) = 2 * _wall_temperature[id] - field.t(i, j + 1);
                }
            } else if (cell_borders[0] ==
                       border_position::BOTTOM) { // Border is in BOTTOM position (fluid cell under ghost cell)
                field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j - 1);
                field.v(i, j - 1) = 0;
                field.p(i, j) = field.p(i, j - 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j - 1) = field.v(i, j - 1);

                field.ca(i, j) = field.ca(i, j - 1);
                field.cb(i, j) = field.cb(i, j - 1);
                field.cc(i, j) = field.cc(i, j - 1);

                if (!_wall_temperature.empty()) {
                    field.t(i, j) = 2 * _wall_temperature[id] - field.t(i, j - 1);
                }
            } else if (cell_borders[0] ==
                       border_position::LEFT) { // Border is in LEFT position (fluid cell left from ghost cell)
                field.u(i - 1, j) = 0;
                field.v(i, j) = 2 * _wall_velocity[id] - field.v(i - 1, j);
                field.p(i, j) = field.p(i - 1, j);
                field.f(i - 1, j) = field.u(i - 1, j);
                field.g(i, j) = field.v(i, j);

                field.ca(i, j) = field.ca(i - 1, j);
                field.cb(i, j) = field.cb(i - 1, j);
                field.cc(i, j) = field.cc(i - 1, j);

                if (!_wall_temperature.empty()) {
                    field.t(i, j) = 2 * _wall_temperature[id] - field.t(i - 1, j);
                }
            } else if (cell_borders[0] ==
                       border_position::RIGHT) { // Border is in RIGHT position (fluid cell right from ghost cell)
                field.u(i, j) = 0;
                field.v(i, j) = 2 * _wall_velocity[id] - field.v(i + 1, j);
                field.p(i, j) = field.p(i + 1, j);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);

                field.ca(i, j) = field.ca(i + 1, j);
                field.cb(i, j) = field.cb(i + 1, j);
                field.cc(i, j) = field.cc(i + 1, j);

                if (!_wall_temperature.empty()) {
                    field.t(i, j) = 2 * _wall_temperature[id] - field.t(i + 1, j);
                }
            }
        }
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

InFlowBoundary::InFlowBoundary(std::vector<Cell *> cells, double inflow_velocity, double inflow_concentration_a,
                               double inflow_concentration_b)
    : _cells(cells), _inflow_concentration_a(inflow_concentration_a), _inflow_concentration_b(inflow_concentration_b) {
    _inflow_velocity.insert(std::pair(boundary_ids::inflow_cell_id, inflow_velocity));
    _inflow_velocity.insert(std::pair(boundary_ids::inlet_a_cell_id, inflow_velocity));
    _inflow_velocity.insert(std::pair(boundary_ids::inlet_b_cell_id, inflow_velocity));
}

InFlowBoundary::InFlowBoundary(std::vector<Cell *> cells, double inflow_velocity, double inflow_temperature,
                               double inflow_concentration_a, double inflow_concentration_b)
    : _cells(cells), _inflow_concentration_a(inflow_concentration_a), _inflow_concentration_b(inflow_concentration_b) {
    _inflow_velocity.insert(std::pair(boundary_ids::inflow_cell_id, inflow_velocity));
    _inflow_velocity.insert(std::pair(boundary_ids::inlet_a_cell_id, inflow_velocity));
    _inflow_velocity.insert(std::pair(boundary_ids::inlet_b_cell_id, inflow_velocity));

    _inflow_temperature.insert(std::pair(boundary_ids::inflow_cell_id, inflow_temperature));
    _inflow_temperature.insert(std::pair(boundary_ids::inlet_a_cell_id, inflow_temperature));
    _inflow_temperature.insert(std::pair(boundary_ids::inlet_b_cell_id, inflow_temperature));
}

InFlowBoundary::InFlowBoundary(std::vector<Cell *> cells, std::map<int, double> inflow_velocity,
                               std::map<int, double> inflow_temperature)
    : _cells(cells), _inflow_velocity(inflow_velocity), _inflow_temperature(inflow_temperature) {}

void InFlowBoundary::apply(Fields &field) {
    for (const auto cell : _cells) {
        std::vector<border_position> cell_borders = cell->borders();
        int i = cell->i();
        int j = cell->j();
        int id = cell->cell_id();

        if (cell_borders.size() == 1) {

            if (cell_borders[0] == border_position::TOP) { // Border is in TOP position
                field.u(i, j) = -field.u(i, j + 1);
                field.v(i, j) = _inflow_velocity[id];
                field.p(i, j) = field.p(i, j + 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);

                if (id == 4) {
                    field.ca(i, j) = _inflow_concentration_a;
                    field.cb(i, j) = 0;
                    field.cc(i, j) = 0;
                } else if (id == 5) {
                    field.ca(i, j) = 0;
                    field.cb(i, j) = _inflow_concentration_b;
                    field.cc(i, j) = 0;
                }

                if (!_inflow_temperature.empty()) {
                    field.t(i, j) = 2 * _inflow_temperature[id] - field.t(i, j + 1);
                }
            } else if (cell_borders[0] == border_position::BOTTOM) { // Border is in BOTTOM position
                field.u(i, j) = -field.u(i, j - 1);
                field.v(i, j - 1) = -_inflow_velocity[id];
                field.p(i, j) = field.p(i, j - 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);

                if (id == 4) {
                    field.ca(i, j) = _inflow_concentration_a;
                    field.cb(i, j) = 0;
                    field.cc(i, j) = 0;
                } else if (id == 5) {
                    field.ca(i, j) = 0;
                    field.cb(i, j) = _inflow_concentration_b;
                    field.cc(i, j) = 0;
                }

                if (!_inflow_temperature.empty()) {
                    field.t(i, j) = 2 * _inflow_temperature[id] - field.t(i, j - 1);
                }
            } else if (cell_borders[0] == border_position::LEFT) { // Border is in LEFT position
                field.u(i - 1, j) = -_inflow_velocity[id];
                field.v(i, j) = -field.v(i - 1, j);
                field.p(i, j) = field.p(i - 1, j);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);

                if (id == 4) {
                    field.ca(i, j) = _inflow_concentration_a;
                    field.cb(i, j) = 0;
                    field.cc(i, j) = 0;
                } else if (id == 5) {
                    field.ca(i, j) = 0;
                    field.cb(i, j) = _inflow_concentration_b;
                    field.cc(i, j) = 0;
                }

                if (!_inflow_temperature.empty()) {
                    field.t(i, j) = 2 * _inflow_temperature[id] - field.t(i - 1, j);
                }
            } else if (cell_borders[0] == border_position::RIGHT) { // Border is in RIGHT position
                field.u(i, j) = _inflow_velocity[id];
                field.v(i, j) = -field.v(i + 1, j);
                field.p(i, j) = field.p(i + 1, j);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);

                if (id == 4) {
                    field.ca(i, j) = _inflow_concentration_a;
                    field.cb(i, j) = 0;
                    field.cc(i, j) = 0;
                } else if (id == 5) {
                    field.ca(i, j) = 0;
                    field.cb(i, j) = _inflow_concentration_b;
                    field.cc(i, j) = 0;
                }

                if (!_inflow_temperature.empty()) {
                    field.t(i, j) = 2 * _inflow_temperature[id] - field.t(i + 1, j);
                }
            }
        }
    }
}

OutFlowBoundary::OutFlowBoundary(std::vector<Cell *> cells, double outflow_pressure) : _cells(cells) {
    _outflow_pressure.insert(std::pair(boundary_ids::outflow_cell_id, outflow_pressure));
}

OutFlowBoundary::OutFlowBoundary(std::vector<Cell *> cells, std::map<int, double> outflow_pressure)
    : _cells(cells), _outflow_pressure(outflow_pressure) {}

void OutFlowBoundary::apply(Fields &field) {
    for (const auto cell : _cells) {
        std::vector<border_position> cell_borders = cell->borders();
        int i = cell->i();
        int j = cell->j();
        int id = cell->cell_id();

        if (cell_borders.size() == 1) {

            if (cell_borders[0] == border_position::TOP) { // Border is in TOP position
                field.u(i, j) = field.u(i, j + 1);
                field.v(i, j) = field.v(i, j + 1);
                field.p(i, j) = _outflow_pressure[id];
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
                field.t(i, j) = field.t(i, j + 1);

                field.ca(i, j) = field.ca(i, j + 1);
                field.cb(i, j) = field.cb(i, j + 1);
                field.cc(i, j) = field.cc(i, j + 1);

            } else if (cell_borders[0] == border_position::BOTTOM) { // Border is in BOTTOM position
                field.u(i, j) = field.u(i, j - 1);
                field.v(i, j) = field.v(i, j - 1);
                field.p(i, j) = _outflow_pressure[id];
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
                field.t(i, j) = field.t(i, j - 1);

                field.ca(i, j) = field.ca(i, j - 1);
                field.cb(i, j) = field.cb(i, j - 1);
                field.cc(i, j) = field.cc(i, j - 1);

            } else if (cell_borders[0] == border_position::LEFT) { // Border is in LEFT position
                field.u(i, j) = field.u(i - 1, j);
                field.v(i, j) = field.v(i - 1, j);
                field.p(i, j) = _outflow_pressure[id];
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
                field.t(i, j) = field.t(i - 1, j);

                field.ca(i, j) = field.ca(i - 1, j);
                field.cb(i, j) = field.cb(i - 1, j);
                field.cc(i, j) = field.cc(i - 1, j);

            } else if (cell_borders[0] == border_position::RIGHT) { // Border is in RIGHT position
                field.u(i, j) = field.u(i + 1, j);
                field.v(i, j) = field.v(i + 1, j);
                field.p(i, j) = _outflow_pressure[id];
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
                field.t(i, j) = field.t(i + 1, j);

                field.ca(i, j) = field.ca(i + 1, j);
                field.cb(i, j) = field.cb(i + 1, j);
                field.cc(i, j) = field.cc(i + 1, j);
            }
        }
    }
}