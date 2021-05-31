#pragma once

#include "Datastructures.hpp"
#include "Discretization.hpp"
#include "Grid.hpp"

/**
 * @brief Class of container and modifier for the physical fields
 *
 */
class Fields {
  public:
    Fields() = default;

    /**
     * @brief Constructor for the fields
     *
     * @param[in] kinematic viscosity
     * @param[in] initial timestep size
     * @param[in] adaptive timestep coefficient
     * @param[in] number of cells in x direction
     * @param[in] number of cells in y direction
     * @param[in] initial x-velocity
     * @param[in] initial y-velocity
     * @param[in] initial pressure
     *
     */
    Fields(double _nu, double _dt, double _tau, double _alpha, double _beta, std::vector<Cell *> cells, int imax,
           int jmax, double UI, double VI, double PI, double TI, std::string _energy_eq, double GX, double GY);

    /**
     * @brief Calculates the convective and diffusive fluxes in x and y
     * direction based on explicit discretization of the momentum equations
     *
     * @param[in] grid in which the fluxes are calculated
     *
     */
    void calculate_fluxes(Grid &grid);

    /**
     * @brief Right hand side calculations using the fluxes for the pressure
     * Poisson equation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_rs(Grid &grid);

    /**
     * @brief Velocity calculation using pressure values
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_velocities(Grid &grid);

    /**
     * @brief Copying matrices
     *
     * @param[in] grid in which the calculations are done
     * @param[in] data to be copied and copied to
     *
     */
    void copy_matrix(Grid &grid, const Matrix<double> &FROM, Matrix<double> &TO);

    /**
     * @brief Temperature calculation based on explicit discretization of the heat equation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_temperature(Grid &grid);

    /**
     * @brief Adaptive step size calculation using x-velocity condition,
     * y-velocity condition and CFL condition
     *
     * @param[in] grid in which the calculations are done
     *
     */
    double calculate_dt(Grid &grid);

    /**
     * @brief Finding maximum element of matrices
     *
     * @param[in] grid in which the calculations are done
     * @param[in] data the maximum is to be found in
     *
     */
    double find_max(const Matrix<double> &M, const Grid &grid);

    /// x-velocity index based access and modify
    double &u(int i, int j);

    /// y-velocity index based access and modify
    double &v(int i, int j);

    /// pressure index based access and modify
    double &p(int i, int j);

    /// RHS index based access and modify
    double &rs(int i, int j);

    /// x-momentum flux index based access and modify
    double &f(int i, int j);

    /// y-momentum flux index based access and modify
    double &g(int i, int j);

    /// temperature index based access and modify
    double &t(int i, int j);

    /// temporary temperature index based access and modify
    double &temp(int i, int j);

    /// current mean of the u velocity field
    double &u_avg();

    /// current mean of the v velocity field
    double &v_avg();

    /// current mean of the pressure field
    double &p_avg();

    /// current mean of the pressure field
    void set_p_avg(double p_avg);

    /// current mean of the temperature field
    double &t_avg();

    /// get timestep size
    double dt() const;

    /// calculate mean value of field matrix for relative update calculation
    double calculate_mean(Grid &grid);

    /// pressure matrix access and modify
    Matrix<double> &p_matrix();

  private:
    /// x-velocity matrix
    Matrix<double> _U;
    /// x-velocity matrix
    Matrix<double> _V;
    /// pressure matrix
    Matrix<double> _P;
    /// x-momentum flux matrix
    Matrix<double> _F;
    /// y-momentum flux matrix
    Matrix<double> _G;
    /// right hand side matrix
    Matrix<double> _RS;
    /// temperature matrix
    Matrix<double> _T;
    /// temporary temperature update
    Matrix<double> _TEMP;

    /// kinematic viscosity
    double _nu;
    /// gravitional accelearation in x direction
    double _gx;
    /// gravitional accelearation in y direction
    double _gy;
    /// timestep size
    double _dt;
    /// adaptive timestep coefficient
    double _tau;
    /// thermal diffusivity
    double _alpha;
    /// thermal expansion coefficient
    double _beta;
    /// fluid cells
    std::vector<Cell *> _cells;
    /// Heat energy on
    std::string _energy_eq{"NONE"};
    /// Average u velocity for relative update calculation
    double _u_avg;

    /// Average v velocity for relative update calculation
    double _v_avg;

    /// Average pressure for relative update calculation
    double _p_avg;

    /// Average temperature for relative update calculation
    double _t_avg;
};
