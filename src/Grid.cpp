#include "Grid.hpp"
#include "Enums.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

Grid::Grid(std::string geom_name, Domain &domain) {
    _domain = domain;
    _cells = Matrix<Cell>(_domain.size_x + 2, _domain.size_y + 2);

    if (geom_name.compare("NONE")) {
        std::vector<std::vector<int>> geom_data(_domain.size_x + 2, std::vector<int>(_domain.size_y + 2, 0));
        parse_geometry_file(geom_name, geom_data);
        assign_cell_types(geom_data);
        geometry_data = geom_data;
    } else {
        build_lid_driven_cavity();
    }
}

void Grid::build_lid_driven_cavity() {
    std::vector<std::vector<int>> geometry_data(_domain.size_x + 2, std::vector<int>(_domain.size_y + 2, 0));

    int array_lid[_domain.domain_size_x + 2][_domain.domain_size_y + 2];

    // Build the Lid Geometry in the array
    for (int i = 0; i < _domain.domain_size_x + 2; ++i) {
        for (int j = 0; j < _domain.domain_size_y + 2; ++j) {
            // Bottom, left and right walls: no-slip
            if (i == 0 || j == 0 || i == _domain.domain_size_x + 1) {
                array_lid[i][j] = LidDrivenCavity::fixed_wall_id;
            }
            // Top wall: moving wall
            else if (j == _domain.domain_size_y + 1) {
                array_lid[i][j] = LidDrivenCavity::moving_wall_id;
            }
            // Inner cells: fluid
            else {
                array_lid[i][j] = 0;
            }
        }
    }

    // Copy the local domain information into geometry_data
    for (int col = _domain.jmin; col < _domain.jmax; ++col) {
        for (int row = _domain.imin; row < _domain.imax; ++row) {
            geometry_data[row - _domain.imin][col - _domain.jmin] = array_lid[row][col];
        }
    }

    assign_cell_types(geometry_data);
}

void Grid::assign_cell_types(std::vector<std::vector<int>> &geometry_data) {
    int i = 0;
    int j = 0;

    for (int j_geom = 0; j_geom <= _domain.size_y + 1; ++j_geom) {
        {
            i = 0;
        }

        for (int i_geom = 0; i_geom <= _domain.size_x + 1; ++i_geom) {
            if ((geometry_data.at(i_geom).at(j_geom) == 0 || geometry_data.at(i_geom).at(j_geom) == 7) &&
                (i == 0 or j == 0 or i == _domain.size_x + 1 or j == _domain.size_y + 1)) {
                _cells(i, j) = Cell(i, j, cell_type::FLUID_BUFFER);
                _fluidbuffer_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 8) {
                _cells(i, j) = Cell(i, j, cell_type::MOVING_WALL, geometry_data.at(i_geom).at(j_geom));
                _moving_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 7) {
                _cells(i, j) = Cell(i, j, cell_type::FLUID, geometry_data.at(i_geom).at(j_geom));
                _fluid_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 3) {
                _cells(i, j) = Cell(i, j, cell_type::FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                _fixed_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 9) {
                _cells(i, j) = Cell(i, j, cell_type::FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                _fixed_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 4) {
                _cells(i, j) = Cell(i, j, cell_type::INFLOW, geometry_data.at(i_geom).at(j_geom));
                _inflow_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 5) {
                _cells(i, j) = Cell(i, j, cell_type::INFLOW, geometry_data.at(i_geom).at(j_geom));
                _inflow_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 1) {
                _cells(i, j) = Cell(i, j, cell_type::INFLOW, geometry_data.at(i_geom).at(j_geom));
                _inflow_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 2) {
                _cells(i, j) = Cell(i, j, cell_type::OUTFLOW, geometry_data.at(i_geom).at(j_geom));
                _outflow_cells.push_back(&_cells(i, j));
            } else {
                _cells(i, j) = Cell(i, j, cell_type::FLUID, geometry_data.at(i_geom).at(j_geom));
                _fluid_cells.push_back(&_cells(i, j));
            }
            ++i;
        }
        ++j;
    }

    // Corner cell neighbour assigment
    // Bottom-Left Corner
    i = 0;
    j = 0;
    _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
    _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
    if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID or
        _cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID_BUFFER) {
        _cells(i, j).add_border(border_position::TOP);
    }
    if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID or
        _cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID_BUFFER) {
        _cells(i, j).add_border(border_position::RIGHT);
    }
    // Top-Left Corner
    i = 0;
    j = _domain.size_y + 1;
    _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
    _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
    if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID or
        _cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID_BUFFER) {
        _cells(i, j).add_border(border_position::BOTTOM);
    }
    if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID or
        _cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID_BUFFER) {
        _cells(i, j).add_border(border_position::RIGHT);
    }

    // Top-Right Corner
    i = _domain.size_x + 1;
    j = Grid::_domain.size_y + 1;
    _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
    _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
    if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID or
        _cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID_BUFFER) {
        _cells(i, j).add_border(border_position::BOTTOM);
    }
    if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID or
        _cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID_BUFFER) {
        _cells(i, j).add_border(border_position::LEFT);
    }

    // Bottom-Right Corner
    i = Grid::_domain.size_x + 1;
    j = 0;
    _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
    _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
    if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID or
        _cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID_BUFFER) {
        _cells(i, j).add_border(border_position::TOP);
    }
    if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID or
        _cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID_BUFFER) {
        _cells(i, j).add_border(border_position::LEFT);
    }
    // Bottom cells
    j = 0;
    for (int i = 1; i < _domain.size_x + 1; ++i) {
        _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
        _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
        _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID or
            _cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID_BUFFER) {
            _cells(i, j).add_border(border_position::RIGHT);
        }
        if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID or
            _cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID_BUFFER) {
            _cells(i, j).add_border(border_position::LEFT);
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID or
            _cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID_BUFFER) {
            _cells(i, j).add_border(border_position::TOP);
        }
    }

    // Top Cells
    j = Grid::_domain.size_y + 1;

    for (int i = 1; i < _domain.size_x + 1; ++i) {
        _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
        _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
        _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID or
            _cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID_BUFFER) {
            _cells(i, j).add_border(border_position::RIGHT);
        }
        if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID or
            _cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID_BUFFER) {
            _cells(i, j).add_border(border_position::LEFT);
        }
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID or
            _cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID_BUFFER) {
            _cells(i, j).add_border(border_position::BOTTOM);
        }
    }

    // Left Cells
    i = 0;
    for (int j = 1; j < _domain.size_y + 1; ++j) {
        _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
        _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
        _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID or
            _cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID_BUFFER) {
            _cells(i, j).add_border(border_position::RIGHT);
        }
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID or
            _cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID_BUFFER) {
            _cells(i, j).add_border(border_position::BOTTOM);
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID or
            _cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID_BUFFER) {
            _cells(i, j).add_border(border_position::TOP);
        }
    }
    // Right Cells
    i = Grid::_domain.size_x + 1;
    for (int j = 1; j < _domain.size_y + 1; ++j) {
        _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
        _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
        _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
        if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID or
            _cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID_BUFFER) {
            _cells(i, j).add_border(border_position::LEFT);
        }
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID or
            _cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID_BUFFER) {
            _cells(i, j).add_border(border_position::BOTTOM);
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID or
            _cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID_BUFFER) {
            _cells(i, j).add_border(border_position::TOP);
        }
    }

    // Inner cells
    for (int i = 1; i < _domain.size_x + 1; ++i) {
        for (int j = 1; j < _domain.size_y + 1; ++j) {
            _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
            _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
            _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
            _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);

            if (_cells(i, j).type() != cell_type::FLUID) {
                if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID or
                    _cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID_BUFFER) {
                    _cells(i, j).add_border(border_position::LEFT);
                }
                if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID or
                    _cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID_BUFFER) {
                    _cells(i, j).add_border(border_position::RIGHT);
                }
                if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID or
                    _cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID_BUFFER) {
                    _cells(i, j).add_border(border_position::BOTTOM);
                }
                if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID or
                    _cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID_BUFFER) {
                    _cells(i, j).add_border(border_position::TOP);
                }
            }
        }
    }
}

void Grid::parse_geometry_file(std::string filedoc, std::vector<std::vector<int>> &geometry_data) {

    int numcols, numrows, depth;

    std::ifstream infile(filedoc);
    std::stringstream ss;
    std::string inputLine = "";

    // First line : version
    getline(infile, inputLine);
    if (inputLine.compare("P2") != 0) {
        std::cerr << "First line of the PGM file should be P2" << std::endl;
    }

    // Second line : comment
    getline(infile, inputLine);

    // Continue with a stringstream
    ss << infile.rdbuf();
    // Third line : size
    ss >> numrows >> numcols;
    // Fourth line : depth
    ss >> depth;

    int original_array[numrows][numcols];

    // Following lines : data
    for (int col = numcols - 1; col > -1; --col) {
        for (int row = 0; row < numrows; ++row) {
            ss >> original_array[row][col];
        }
    }

    // Refine the grid
    int numrows_ref = _domain.ref_factor * (numrows - 2) + 2;
    int numcols_ref = _domain.ref_factor * (numcols - 2) + 2;
    int array[numrows_ref][numcols_ref];

    for (int j_small = 1; j_small < numcols - 1; j_small++) { // inner cells
        for (int i_small = 1; i_small < numrows - 1; i_small++) {
            for (int i = (_domain.ref_factor * i_small - _domain.ref_factor + 1); i <= (_domain.ref_factor * i_small);
                 i++) {
                for (int j = (_domain.ref_factor * j_small - _domain.ref_factor + 1);
                     j <= (_domain.ref_factor * j_small); j++) {
                    array[i][j] = original_array[i_small][j_small];
                }
            }
        }
    }

    for (int j_small = 1; j_small < numcols - 1; j_small++) {
        for (int j = (_domain.ref_factor * j_small - _domain.ref_factor + 1); j <= (_domain.ref_factor * j_small);
             j++) { // upper and lower row
            array[0][j] = original_array[0][j_small];
            array[numrows_ref - 1][j] = original_array[numrows - 1][j_small];
        }
    }

    for (int i_small = 0; i_small < numrows - 1; i_small++) {
        for (int i = (_domain.ref_factor * i_small - _domain.ref_factor); i <= (_domain.ref_factor * i_small + 1);
             i++) {
            array[i][0] = original_array[i_small][0];
            array[i][numcols_ref - 1] = original_array[i_small][numcols - 1];
        }
    }

    // Copy the local domain information into geometry_data
    for (int col = _domain.jmin; col < _domain.jmax; ++col) {
        for (int row = _domain.imin; row < _domain.imax; ++row) {
            geometry_data[row - _domain.imin][col - _domain.jmin] = array[row][col];
        }
    }

    infile.close();
}

int Grid::imax() const { return _domain.size_x; }
int Grid::jmax() const { return _domain.size_y; }

int Grid::imaxb() const { return _domain.size_x + 2; }
int Grid::jmaxb() const { return _domain.size_y + 2; }

Cell Grid::cell(int i, int j) const { return _cells(i, j); }

double Grid::dx() const { return _domain.dx; }

double Grid::dy() const { return _domain.dy; }

const Domain &Grid::domain() const { return _domain; }

const std::vector<std::vector<int>> &Grid::get_geometry_data() const { return geometry_data; }

const std::vector<Cell *> &Grid::fluid_cells() const { return _fluid_cells; }

const std::vector<Cell *> &Grid::fixed_wall_cells() const { return _fixed_wall_cells; }

const std::vector<Cell *> &Grid::fluidbuffer_cells() const { return _fluidbuffer_cells; }

const std::vector<Cell *> &Grid::outflow_cells() const { return _outflow_cells; }

const std::vector<Cell *> &Grid::inflow_cells() const { return _inflow_cells; }

const std::vector<Cell *> &Grid::free_slip() const { return _free_slip_cells; }

const std::vector<Cell *> &Grid::moving_wall_cells() const { return _moving_wall_cells; }
