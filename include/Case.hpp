#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Boundary.hpp"
#include "Communication.hpp"
#include "Discretization.hpp"
#include "Domain.hpp"
#include "Fields.hpp"
#include "Grid.hpp"
#include "PressureSolver.hpp"

/**
 * @brief Class to hold and orchestrate the simulation flow.
 *
 */
class Case {
  public:
    /**
     * @brief Parallel constructor for the Case.
     *
     * Reads input file, creates Fields, Grid, Boundary, Solver and sets
     * Discretization parameters Creates output directory
     *
     * @param[in] Input file name
     */
    Case(std::string file_name, int argn, char **args);

    /**
     * @brief Main function to simulate the flow until the end time.
     *
     * Calculates the fluxes
     * Calculates the right hand side
     * Solves pressure
     * Calculates velocities
     * Outputs the solution files
     */
    void simulate();

  private:
    /// Plain case name without paths
    std::string _case_name;
    /// Output directiory name
    std::string _dict_name;
    /// Geometry file name
    std::string _geom_name{"NONE"};
    /// Relative input file path
    std::string _prefix;
    /// Heat energy on
    std::string energy_eq{"NONE"};
    /// Wall (id 3) temperature
    double wall_temp_3;
    /// Wall (id 9) temperature
    double wall_temp_9;

    /// Wall (id 3) velocity
    double wall_vel_3;

    /// Number of processes in x direction
    int iproc{1};
    /// Number of processes in y direction
    int jproc{1};

    /// Thermal conductivity
    double kappa;

    /// Simulation time
    double _t_end;
    /// Solution file outputting frequency
    double _output_freq;
    /// Refinement factor
    int ref_factor{1};

    Communication _communication;
    Fields _field;
    Grid _grid;
    Discretization _discretization;
    Domain domain;
    std::unique_ptr<PressureSolver> _pressure_solver;
    std::vector<std::unique_ptr<Boundary>> _boundaries;

    /// Solver convergence tolerance
    double _tolerance;

    /// Maximum number of iterations for the solver
    int _max_iter;

    /**
     * @brief Creating file names from given input data file
     *
     * Extracts path of the case file and creates code-readable file names
     * for outputting directory and geometry file.
     *
     * @param[in] input data file
     */
    void set_file_names(std::string file_name);

    /**
     * @brief Solution file outputter
     *
     * Outputs the solution files in .vtk format. Ghost cells are excluded.
     * Pressure is cell variable while velocity is point variable while being
     * interpolated to the cell faces
     *
     * @param[in] Timestep of the solution
     */
    void output_vtk(int file_number, int my_rank);

    void build_domain(Domain &domain, int imax, int jmax);
};
