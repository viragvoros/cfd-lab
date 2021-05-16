#pragma once

#include <vector>

#include "Cell.hpp"
#include "Fields.hpp"

/**
 * @brief Abstact of boundary conditions.
 *
 * This class patches the physical values to the given field.
 */
class Boundary {
  public:
    /**
     * @brief Main method to patch the boundary conditons to given field and
     * grid
     *
     * @param[in] Field to be applied
     */

    // We changed signature to implement assert, see implementation in Boundary.cpp
    virtual void apply(Fields &field) = 0;
    virtual ~Boundary() = default;
};

/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure
 */
class FixedWallBoundary : public Boundary {
  public:
    FixedWallBoundary(std::vector<Cell *> cells);
    FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature);
    FixedWallBoundary(std::vector<Cell *> cells, double wall_temperature);
    virtual ~FixedWallBoundary() = default;
    virtual void apply(Fields &field);

  private:
    std::vector<Cell *> _cells;
    std::map<int, double> _wall_temperature;
};

/**
 * @brief Moving wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities for the given velocity parallel to the fluid,
 * Neumann for pressure
 */
class MovingWallBoundary : public Boundary {
  public:
    MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity);
    MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                       std::map<int, double> wall_temperature);
    virtual ~MovingWallBoundary() = default;
    virtual void apply(Fields &field);

  private:
    std::vector<Cell *> _cells;
    std::map<int, double> _wall_velocity;
    std::map<int, double> _wall_temperature;
};

/**
 * @brief Inflow cell boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities for the given velocity parallel to the fluid,
 * Neumann for pressure
 */
class InFlowBoundary : public Boundary {
  public:
    InFlowBoundary(std::vector<Cell *> cells, double inflow_velocity);
    InFlowBoundary(std::vector<Cell *> cells, double inflow_velocity, double inflow_temperature);
    InFlowBoundary(std::vector<Cell *> cells, std::map<int, double> inflow_velocity,
                       std::map<int, double> inflow_temperature);
    virtual ~InFlowBoundary() = default;
    virtual void apply(Fields &field);

  private:
    std::vector<Cell *> _cells;
    std::map<int, double> _inflow_velocity;
    std::map<int, double> _inflow_temperature;
};

/**
 * @brief Outflow cell boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities for the given velocity parallel to the fluid,
 * Neumann for pressure
 */
class OutFlowBoundary : public Boundary {
  public:
    OutFlowBoundary(std::vector<Cell *> cells, double outflow_velocity);
    OutFlowBoundary(std::vector<Cell *> cells, double outflow_velocity, double outflow_temperature);
    OutFlowBoundary(std::vector<Cell *> cells, std::map<int, double> outflow_velocity,
                       std::map<int, double> outflow_temperature);
    virtual ~OutFlowBoundary() = default;
    virtual void apply(Fields &field);

  private:
    std::vector<Cell *> _cells;
    std::map<int, double> _outflow_velocity;
    std::map<int, double> _outflow_temperature;
};
