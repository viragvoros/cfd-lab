#include "PressureSolver.hpp"

#include <cmath>
#include <iostream>
#include <mpi.h>

SOR::SOR(double omega) : _omega(omega) {}

double SOR::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {

    double dx = grid.dx();
    double dy = grid.dy();

    // Additional application of boundary conditions as mentioned in tutorial.
    for (auto &boundary : boundaries) {
        boundary->apply(field);
    }

    double coeff = _omega / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy))); // = _omega * h^2 / 4.0, if dx == dy == h
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        field.p(i, j) = (1.0 - _omega) * field.p(i, j) +
                        coeff * (Discretization::sor_helper(field.p_matrix(), i, j) - field.rs(i, j));
    }

    double res = 0.0;
    double rloc = 0.0;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += (val * val);
    }

    int num_fluid_cell = grid.fluid_cells().size();

    MPI_Allreduce(MPI_IN_PLACE, &num_fluid_cell, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &rloc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    {
        res = rloc / num_fluid_cell;
        res = std::sqrt(res);
    }
    return res;
}
