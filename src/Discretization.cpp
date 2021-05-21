#include "Discretization.hpp"
#include <iostream>

#include <cmath>

double Discretization::_dx = 0.0;
double Discretization::_dy = 0.0;
double Discretization::_gamma = 0.0;

Discretization::Discretization(double dx, double dy, double gamma) {
    _dx = dx;
    _dy = dy;
    _gamma = gamma;
}

double Discretization::convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {

    double dUU_dx = 1 / _dx * (std::pow(interpolate(U, i, j, 1, 0), 2.0) - std::pow(interpolate(U, i, j, -1, 0), 2.0)) +
                    _gamma / _dx *
                        (std::abs(interpolate(U, i, j, 1, 0)) * (U(i, j) - U(i + 1, j)) * 0.5 -
                         std::abs(interpolate(U, i, j, -1, 0)) * (U(i - 1, j) - U(i, j)) * 0.5);

    double dUV_dy = 1 / _dy *
                        ((interpolate(V, i, j, 1, 0) * interpolate(U, i, j, 0, 1)) -
                         (interpolate(V, i, j - 1, 1, 0) * interpolate(U, i, j, 0, -1))) +
                    _gamma / _dy *
                        ((std::abs(interpolate(V, i, j, 1, 0)) * (U(i, j) - U(i, j + 1)) * 0.5) -
                         (std::abs(interpolate(V, i, j - 1, 1, 0)) * (U(i, j - 1) - U(i, j)) * 0.5));

    double result = dUU_dx + dUV_dy;
    return result;
}

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {

    double dVV_dy = 1 / _dy * (std::pow(interpolate(V, i, j, 0, 1), 2.0) - std::pow(interpolate(V, i, j, 0, -1), 2.0)) +
                    _gamma / _dy *
                        (std::abs(interpolate(V, i, j, 0, 1)) * (V(i, j) - V(i, j + 1)) * 0.5 -
                         std::abs(interpolate(V, i, j, 0, -1)) * (V(i, j - 1) - V(i, j)) * 0.5);

    double dUV_dx = 1 / _dx *
                        ((interpolate(U, i, j, 0, 1) * interpolate(V, i, j, 1, 0)) -
                         (interpolate(U, i - 1, j, 0, 1) * interpolate(V, i, j, -1, 0))) +
                    _gamma / _dx *
                        ((std::abs(interpolate(U, i, j, 0, 1)) * (V(i, j) - V(i + 1, j)) * 0.5) -
                         (std::abs(interpolate(U, i - 1, j, 0, 1)) * (V(i - 1, j) - V(i, j)) * 0.5));

    double result = dVV_dy + dUV_dx;
    return result;
}

double Discretization::convection_t(const Matrix<double> &T, const Matrix<double> &U, const Matrix<double> &V, int i,
                                    int j) {

    double dUT_dx =
        1 / _dx * (U(i, j) * interpolate(T, i, j, 1, 0) - U(i - 1, j) * interpolate(T, i, j, -1, 0)) +
        _gamma / _dx *
            (std::abs(U(i, j)) * (T(i, j) - T(i + 1, j)) * 0.5 - std::abs(U(i - 1, j)) * (T(i - 1, j) - T(i, j)) * 0.5);

    double dVT_dy =
        1 / _dy * (V(i, j) * interpolate(T, i, j, 0, 1) - V(i, j - 1) * interpolate(T, i, j, 0, -1)) +
        _gamma / _dy *
            (std::abs(V(i, j)) * (T(i, j) - T(i, j + 1)) * 0.5 - std::abs(V(i, j - 1)) * (T(i, j - 1) - T(i, j)) * 0.5);

    double result = dUT_dx + dVT_dy;
    return result;
}

double Discretization::diffusion(const Matrix<double> &A, int i, int j) {
    double result = ((A(i + 1, j) - 2.0 * A(i, j) + A(i - 1, j)) / (_dx * _dx) +
                     (A(i, j + 1) - 2.0 * A(i, j) + A(i, j - 1)) / (_dy * _dy));
    return result;
}

double Discretization::laplacian(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (_dx * _dx) +
                    (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) + (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

// Take the "average" of two values
double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {
    double interpolated = (A(i, j) + A(i + i_offset, j + j_offset)) / 2;
    return interpolated;
}