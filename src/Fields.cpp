#include "Fields.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

Fields::Fields(double nu, double dt, double tau, double alpha, double beta, std::vector<Cell *> cells, int imax,
               int jmax, double UI, double VI, double PI, double TI, std::string energy_eq)
    : _nu(nu), _dt(dt), _tau(tau), _alpha(alpha), _beta(beta), _cells(cells), _energy_eq(energy_eq) {
    _U = Matrix<double>(imax + 2, jmax + 2);
    _V = Matrix<double>(imax + 2, jmax + 2);
    _P = Matrix<double>(imax + 2, jmax + 2);
    _T = Matrix<double>(imax + 2, jmax + 2);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);

    for (const auto &cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->type() == cell_type::FLUID) {
            u(i, j) = UI;
            v(i, j) = VI;
            p(i, j) = PI;
            t(i, j) = TI;
        }
    }
    // for (int jx = 0; jx < jmax + 2; jx++ ){

    //     for (int ix = 0; ix < imax + 2; ix++) {
    //         std::cout << _P(ix, jx) << " " ;
    //     }
    //     std::cout << "\n";
    //  }
}

void Fields::calculate_fluxes(Grid &grid) {
    // std::cout << _energy_eq << std::endl;
    if (_energy_eq.compare("NONE") == 0) {
        std::cout << "NO_TEMP_FLUX" << std::endl;
        for (int j = 1; j <= grid.jmax(); j++) {
            for (int i = 1; i <= (grid.imax() - 1); i++) {
                f(i, j) = u(i, j) + _dt * (_nu * Discretization::diffusion(_U, i, j) -
                                           Discretization::convection_u(_U, _V, i, j) + _gx);
            }
        }
        for (int j = 1; j <= (grid.jmax() - 1); j++) {
            for (int i = 1; i <= grid.imax(); i++) {
                g(i, j) = v(i, j) + _dt * (_nu * Discretization::diffusion(_V, i, j) -
                                           Discretization::convection_v(_U, _V, i, j) + _gy);
            }
        }
    }

    else {
        std::cout << "WITH_TEMP_FLUX" << std::endl;
        for (int j = 1; j <= grid.jmax(); j++) {
            for (int i = 1; i <= (grid.imax() - 1); i++) {
                f(i, j) = u(i, j) +
                          _dt * (_nu * Discretization::diffusion(_U, i, j) -
                                 Discretization::convection_u(_U, _V, i, j) + _gx) -
                          _beta * _dt * 0.5 * (t(i, j) + t(i + 1, j)) * _gx;
            }
        }
        for (int j = 1; j <= (grid.jmax() - 1); j++) {
            for (int i = 1; i <= grid.imax(); i++) {
                g(i, j) = v(i, j) +
                          _dt * (_nu * Discretization::diffusion(_V, i, j) -
                                 Discretization::convection_v(_U, _V, i, j) + _gy) -
                          _beta * _dt * 0.5 * (t(i, j) + t(i, j + 1)) * _gy;
            }
        }
    }
}

// Calculate righ-hand-side of PPE
void Fields::calculate_rs(Grid &grid) {
    for (int j = 1; j <= grid.jmax(); j++) {
        for (int i = 1; i <= grid.imax(); i++) {
            rs(i, j) = 1 / _dt * ((f(i, j) - f(i - 1, j)) / grid.dx() + (g(i, j) - g(i, j - 1)) / grid.dy());
        }
    }
}

void Fields::calculate_velocities(Grid &grid) {
    for (int j = 1; j <= grid.jmax(); j++) {
        for (int i = 1; i <= (grid.imax() - 1); i++) {
            u(i, j) = f(i, j) - _dt / grid.dx() * (p(i + 1, j) - p(i, j));
        }
    }
    for (int j = 1; j <= (grid.jmax() - 1); j++) {
        for (int i = 1; i <= grid.imax(); i++) {
            v(i, j) = g(i, j) - _dt / grid.dy() * (p(i, j + 1) - p(i, j));
        }
    }
}

void Fields::calculate_temperature(Grid &grid) {
    if (_energy_eq.compare("NONE") != 0) {
        std::cout << "CALC_TEMP" << std::endl;
        for (int j = 1; j <= grid.jmax(); j++) {
            for (int i = 1; i <= grid.imax(); i++) {
                t(i, j) = t(i, j) + _dt * (_alpha * Discretization::diffusion(_T, i, j) -
                                           Discretization::convection_t(_T, _U, _V, i, j));
            }
        }
    }
}

// Function was implemented to get maximum value of matrix
double Fields::find_max(const Matrix<double> &M, const int &imaxb, const int &jmaxb) {
    double maximum = 0;
    for (int j = 0; j <= jmaxb; ++j) {
        for (int i = 0; i <= imaxb; ++i) {
            maximum = std::max(M(i, j), maximum);
        }
    }
    return maximum;
}

// Getting optimum dt based on CFL condition
double Fields::calculate_dt(Grid &grid) {
    double val_1 = (1 / (2 * _nu)) * 1 / (1 / (grid.dx() * grid.dx()) + 1 / (grid.dy() * grid.dy()));
    double val_2 = grid.dx() / std::abs(find_max(_U, grid.imax(), grid.jmax()));
    double val_3 = grid.dy() / std::abs(find_max(_V, grid.imax(), grid.jmax()));
    double max_dt;

    if (_energy_eq.compare("NONE") == 0) {
        std::cout << "NO_TEMP_DT" << std::endl;
        max_dt = _tau * std::min({val_1, val_2, val_3});
    } else {
        std::cout << "WITH_TEMP_DT" << std::endl;
        double val_4 = (1 / (2 * _alpha)) * 1 / (1 / (grid.dx() * grid.dx()) + 1 / (grid.dy() * grid.dy()));
        max_dt = _tau * std::min({val_1, val_2, val_3, val_4});
    }

    _dt = max_dt;
    return max_dt;
}

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::t(int i, int j) { return _T(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }