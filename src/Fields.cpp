#include "Fields.hpp"

#include <algorithm>
#include <iostream>

Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI)
    : _nu(nu), _dt(dt), _tau(tau) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI);
    _V = Matrix<double>(imax + 2, jmax + 2, VI);
    _P = Matrix<double>(imax + 2, jmax + 2, PI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
}

void Fields::calculate_fluxes(Grid &grid) {
    for(int j = 1; j < grid.jmax(); j++){
        for(int i = 1; i < grid.imax(); i++){
            f(i, j) = u(i, j) + _dt * (_nu * Discretization::diffusion(_U, i, j) - Discretization::convection_u(_U, _V, i, j) + _gx);
            g(i, j) = v(i, j) + _dt * (_nu * Discretization::diffusion(_V, i, j) - Discretization::convection_v(_U, _V, i, j) + _gy);
        }
    }
}

void Fields::calculate_rs(Grid &grid) {
    for(int j = 0; j < grid.jmaxb(); j++){
        for(int i = 0; i < grid.imaxb(); i++){
            rs(i, j) = 1 / _dt * ((f(i, j) - f(i - 1, j)) / grid.dx() + (g(i, j) - g(i, j - 1)) / grid.dy());
        }
    }
}

void Fields::calculate_velocities(Grid &grid) {
    for(int j = 1; j < (grid.jmax()-1); j++){
        for(int i = 1; i < (grid.imax()-1); i++){
            u(i, j) = f(i, j) + _dt / grid.dx() * (p(i + 1, j) - p(i, j));
            v(i, j) = g(i, j) + _dt / grid.dy() * (p(i, j + 1) - p(i, j));
        }
    } 
}

double Fields::calculate_dt(Grid &grid) {
    double val_1 = 1 / (2 * _nu) * 1 / (1 / (grid.dx()* grid.dx()) + 1 / (grid.dy()* grid.dy()));
    double val_2 = grid.dx() / std::max(_U);
    double val_3 = grid.dy() / std::max(_V);
    double max_dt = _tau * std::min(val_1, val_2, val_3);
    return max_dt;
}

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }