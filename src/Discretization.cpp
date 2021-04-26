#include "Discretization.hpp"

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
    double dUU = 1 / _dx * (pow((((U(i, j) + U(i + 1, j)) * 0.5)), 2.0) - pow((((U(i - 1, j) + U(i, j)) * 0.5)), 2.0)) +
                    _gamma / _dx * (abs((U(i, j) + U(i + 1, j)) * 0.5) * ((U(i, j) - U(i + 1, j)) * 0.5) - 
                    abs((U(i - 1, j) + U(i, j)) * 0.5) * ((U(i - 1, j) - U(i, j)) * 0.5));

    double dUV = 1 / _dy * ((V(i, j) + V(i + 1, j) * 0.5) * ((U(i, j) + U(i, j + 1)) * 0.5) - 
                    (V(i, j - 1) + V(i + 1, j - 1) * 0.5) * ((U(i, j - 1) + U(i, j)) * 0.5)) + 
                    _gamma / _dy * (abs(V(i, j) + V(i + 1, j) * 0.5) * ((U(i, j) - U(i, j + 1)) * 0.5) - 
                    abs(V(i, j - 1) + V(i + 1, j - 1) * 0.5) * ((U(i, j - 1) - U(i, j)) * 0.5));
    
    double result = dUU + dUV;
    return result;
}

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    double dVV = 1 / _dy * (pow((((V(i, j) + V(i, j + 1)) * 0.5)), 2.0) - pow((((V(i, j - 1) + V(i, j)) * 0.5)), 2.0)) +
                    _gamma / _dy * (abs((V(i, j) + V(i, j + 1)) * 0.5) * ((V(i, j) - V(i, j + 1)) * 0.5) - 
                    abs((V(i, j - 1) + V(i, j)) * 0.5) * ((V(i, j - 1) - V(i, j)) * 0.5));

    double dUV = 1 / _dx * ((U(i, j) + U(i, j + 1) * 0.5) * ((V(i, j) + V(i + 1, j)) * 0.5) - 
                    (U(i - 1, j) + U(i - 1, j + 1) * 0.5) * ((V(i - 1, j) + V(i, j)) * 0.5)) + 
                    _gamma / _dx * (abs(U(i, j) + U(i, j + 1) * 0.5) * ((V(i, j) - V(i + 1, j)) * 0.5) - 
                    abs(U(i - 1, j) + U(i - 1, j + 1) * 0.5) * ((V(i - 1, j) - V(i, j)) * 0.5));
    
    double result = dVV + dUV;
    return result;
}

double Discretization::diffusion(const Matrix<double> &A, int i, int j) {
    double result = (A(i + 1, j) - 2.0 * A(i, j) + A(i - 1, j) / (_dx * _dx) + 
                    A(i, j + 1) - 2.0 * A(i, j) + A(i, j - 1) / (_dy * _dy));
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

double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {}