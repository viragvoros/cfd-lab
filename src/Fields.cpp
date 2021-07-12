#include "Fields.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

Fields::Fields(double nu, double dt, double tau, double alpha, double beta, double diffusivity, double rate_const,
               double order_a, double order_b, double pre_exp_factor, double act_energy, double react_temp_increase,
               std::vector<Cell *> fluid_cells, int imax, int jmax, double UI, double VI, double PI, double TI,
               double CAI, double CBI, double CCI, std::string energy_eq, double GX, double GY)
    : _nu(nu), _dt(dt), _tau(tau), _alpha(alpha), _beta(beta), _diffusivity(diffusivity), _rate_const(rate_const),
      _order_a(order_a), _order_b(order_b), _pre_exp_factor(pre_exp_factor), _act_energy(act_energy),
      _react_temp_increase(react_temp_increase), _fluid_cells(fluid_cells), _gx(GX), _gy(GY) {
    _U = Matrix<double>(imax + 2, jmax + 2);
    _V = Matrix<double>(imax + 2, jmax + 2);
    _P = Matrix<double>(imax + 2, jmax + 2);
    _T = Matrix<double>(imax + 2, jmax + 2);
    _TEMP = Matrix<double>(imax + 2, jmax + 2);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);

    _CA = Matrix<double>(imax + 2, jmax + 2);
    _CB = Matrix<double>(imax + 2, jmax + 2);
    _CC = Matrix<double>(imax + 2, jmax + 2);
    _TEMPCA = Matrix<double>(imax + 2, jmax + 2);
    _TEMPCB = Matrix<double>(imax + 2, jmax + 2);
    _TEMPCC = Matrix<double>(imax + 2, jmax + 2);

    _energy_eq = energy_eq;

    for (const auto &cell : _fluid_cells) {
        int i = cell->i();
        int j = cell->j();

        u(i, j) = UI;
        v(i, j) = VI;
        p(i, j) = PI;
        t(i, j) = TI;
        ca(i, j) = CAI;
        cb(i, j) = CBI;
        cc(i, j) = CCI;

        if (cell->cell_id() == 7) {
            _conversion_cells.push_back(cell);
        }
    }
}

void Fields::calculate_fluxes(Grid &grid) {
    if (_energy_eq.compare("NONE") == 0) {
        for (const auto &cell : _fluid_cells) {
            int i = cell->i();
            int j = cell->j();

            f(i, j) = u(i, j) + _dt * (_nu * Discretization::diffusion(_U, i, j) -
                                       Discretization::convection_u(_U, _V, i, j) + _gx);
            g(i, j) = v(i, j) + _dt * (_nu * Discretization::diffusion(_V, i, j) -
                                       Discretization::convection_v(_U, _V, i, j) + _gy);
        }
    } else {
        for (const auto &cell : _fluid_cells) {
            int i = cell->i();
            int j = cell->j();

            f(i, j) = u(i, j) +
                      _dt * (_nu * Discretization::diffusion(_U, i, j) -
                             Discretization::convection_u(_U, _V, i, j)) - // removed _gx from NS
                      _beta * _dt * 0.5 * (t(i, j) + t(i + 1, j)) * _gx;
            g(i, j) = v(i, j) +
                      _dt * (_nu * Discretization::diffusion(_V, i, j) -
                             Discretization::convection_v(_U, _V, i, j)) - // removed _gy from NS
                      _beta * _dt * 0.5 * (t(i, j) + t(i, j + 1)) * _gy;
        }
    }
}

// Calculate right-hand-side of PPE
void Fields::calculate_rs(Grid &grid) {
    for (const auto &cell : _fluid_cells) {
        int i = cell->i();
        int j = cell->j();

        rs(i, j) = 1 / _dt * ((f(i, j) - f(i - 1, j)) / grid.dx() + (g(i, j) - g(i, j - 1)) / grid.dy());
    }
}

void Fields::calculate_velocities(Grid &grid) {
    for (const auto &cell : _fluid_cells) {
        int i = cell->i();
        int j = cell->j();

        u(i, j) = f(i, j) - _dt / grid.dx() * (p(i + 1, j) - p(i, j));
        _u_avg += std::abs(u(i, j)); // Collect field values for later relative update calculation

        v(i, j) = g(i, j) - _dt / grid.dy() * (p(i, j + 1) - p(i, j));
        _v_avg += std::abs(v(i, j)); // Collect field values for later relative update calculation
    }
    _u_avg = std::sqrt(_u_avg) /
             (grid.jmax() * (grid.imax() - 1)); // Calculate average to prepare for relative update calculation
    _v_avg = std::sqrt(_v_avg) /
             (grid.jmax() * (grid.imax() - 1)); // Calculate average to prepare for relative update calculation
}

// Function was implemented to copy matrices
void Fields::copy_matrix(Grid &grid, const Matrix<double> &FROM, Matrix<double> &TO) {
    for (int i = 0; i <= (grid.imax() + 1); i++) {
        for (int j = 0; j <= (grid.jmax() + 1); j++) {
            TO(i, j) = FROM(i, j);
        }
    }
}

void Fields::calculate_temperature(Grid &grid) {
    if (_energy_eq.compare("NONE") != 0) {
        copy_matrix(grid, _T, _TEMP);
        for (const auto &cell : _fluid_cells) {
            int i = cell->i();
            int j = cell->j();

            t(i, j) = temp(i, j) + _dt * (_alpha * Discretization::diffusion(_TEMP, i, j) -
                                          Discretization::convection_t(_TEMP, _U, _V, i, j));

            if (t(i, j) < 0) {
                t(i, j) = 0;
            }

            _t_avg += std::abs(t(i, j)); // Collect field values for later relative update calculation
        }
        _t_avg = std::sqrt(_t_avg) /
                 (grid.jmax() * (grid.imax() - 1)); // Calculate average to prepare for relative update calculation
    }
}

void Fields::calculate_concentrations(Grid &grid) {
    copy_matrix(grid, _CA, _TEMPCA);
    copy_matrix(grid, _CB, _TEMPCB);
    copy_matrix(grid, _CC, _TEMPCC);
    for (const auto &cell : _fluid_cells) {
        int i = cell->i();
        int j = cell->j();

        ca(i, j) = tempca(i, j) + _dt * (_diffusivity * Discretization::diffusion(_TEMPCA, i, j) -
                                         Discretization::convection_t(_TEMPCA, _U, _V, i, j));
        cb(i, j) = tempcb(i, j) + _dt * (_diffusivity * Discretization::diffusion(_TEMPCB, i, j) -
                                         Discretization::convection_t(_TEMPCB, _U, _V, i, j));
        cc(i, j) = tempcc(i, j) + _dt * (_diffusivity * Discretization::diffusion(_TEMPCC, i, j) -
                                         Discretization::convection_t(_TEMPCC, _U, _V, i, j));

        if (ca(i, j) < 0) {
            ca(i, j) = 0;
        }
        if (cb(i, j) < 0) {
            cb(i, j) = 0;
        }
        if (cc(i, j) < 0) {
            cc(i, j) = 0;
        }
    }
}

// Function to calculate concentrations after reaction
void Fields::react(Grid &grid) {
    if (_energy_eq.compare("NONE") != 0) {
        for (const auto cell : _conversion_cells) {
            int i = cell->i();
            int j = cell->j();
            double tempca = ca(i, j);
            double tempcb = cb(i, j);

            double exponent = -(_act_energy / (8.314 * t(i, j)));
            double rate_const_t = _pre_exp_factor * std::exp(exponent);

            cc(i, j) += rate_const_t * std::pow(tempca, _order_a) * std::pow(tempcb, _order_b);
            ca(i, j) -= rate_const_t * std::pow(tempca, _order_a) * std::pow(tempcb, _order_b);
            cb(i, j) -= rate_const_t * std::pow(tempca, _order_a) * std::pow(tempcb, _order_b);

            if (ca(i, j) < 0) {
                ca(i, j) = 0;
            }
            if (cb(i, j) < 0) {
                cb(i, j) = 0;
            }
            if (cc(i, j) < 0) {
                cc(i, j) = 0;
            }

            t(i, j) += _react_temp_increase;
        }
    } else {
        for (const auto cell : _conversion_cells) {
            int i = cell->i();
            int j = cell->j();
            double tempca = ca(i, j);
            double tempcb = cb(i, j);
            cc(i, j) += _rate_const * std::pow(tempca, _order_a) * std::pow(tempcb, _order_b);
            ca(i, j) -= _rate_const * std::pow(tempca, _order_a) * std::pow(tempcb, _order_b);
            cb(i, j) -= _rate_const * std::pow(tempca, _order_a) * std::pow(tempcb, _order_b);

            if (ca(i, j) < 0) {
                ca(i, j) = 0;
            }
            if (cb(i, j) < 0) {
                cb(i, j) = 0;
            }
            if (cc(i, j) < 0) {
                cc(i, j) = 0;
            }
        }
    }
}

// Function was implemented to get maximum value of matrix
double Fields::find_max(const Matrix<double> &M) {
    double maximum = 0;
    for (const auto &cell : _fluid_cells) {
        int i = cell->i();
        int j = cell->j();
        maximum = std::max(M(i, j), maximum);
    }
    return maximum;
}

// Getting optimum dt based on CFL condition
double Fields::calculate_dt(Grid &grid, const double &max_u, const double &max_v) {
    double val_1 = (1 / (2 * _nu)) * 1 / (1 / (grid.dx() * grid.dx()) + 1 / (grid.dy() * grid.dy()));
    double val_2 = grid.dx() / std::abs(max_u);
    double val_3 = grid.dy() / std::abs(max_v);
    double val_4 = (1 / (2 * _diffusivity)) * 1 / (1 / (grid.dx() * grid.dx()) + 1 / (grid.dy() * grid.dy()));

    double max_dt;

    if (_energy_eq.compare("NONE") == 0) {
        max_dt = _tau * std::min({val_1, val_2, val_3, val_4});
    } else {
        double val_5 = (1 / (2 * _alpha)) * 1 / (1 / (grid.dx() * grid.dx()) + 1 / (grid.dy() * grid.dy()));
        max_dt = _tau * std::min({val_1, val_2, val_3, val_4, val_5});
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
double &Fields::temp(int i, int j) { return _TEMP(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }
double &Fields::ca(int i, int j) { return _CA(i, j); }
double &Fields::cb(int i, int j) { return _CB(i, j); }
double &Fields::cc(int i, int j) { return _CC(i, j); }
double &Fields::tempca(int i, int j) { return _TEMPCA(i, j); }
double &Fields::tempcb(int i, int j) { return _TEMPCB(i, j); }
double &Fields::tempcc(int i, int j) { return _TEMPCC(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }
Matrix<double> &Fields::u_matrix() { return _U; }
Matrix<double> &Fields::v_matrix() { return _V; }
Matrix<double> &Fields::f_matrix() { return _F; }
Matrix<double> &Fields::g_matrix() { return _G; }
Matrix<double> &Fields::t_matrix() { return _T; }
Matrix<double> &Fields::ca_matrix() { return _CA; }
Matrix<double> &Fields::cb_matrix() { return _CB; }
Matrix<double> &Fields::cc_matrix() { return _CC; }

double Fields::dt() const { return _dt; }
double &Fields::u_avg() { return _u_avg; }
double &Fields::v_avg() { return _v_avg; }
double &Fields::p_avg() { return _p_avg; }
double &Fields::t_avg() { return _t_avg; }
void Fields::set_p_avg(double p_avg) { _p_avg = p_avg; }
