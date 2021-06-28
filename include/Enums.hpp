#pragma once

// If no geometry file is provided in the input file, lid driven cavity case
// will run by default. In the Grid.cpp, geometry will be created following
// PGM convention, which is:
// 0: fluid, 3: fixed wall, 8: moving wall
namespace LidDrivenCavity {
const int moving_wall_id = 8;
const int fixed_wall_id = 3;
const double wall_velocity = 1.0;
} // namespace LidDrivenCavity

enum class border_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT,
};

namespace border {
const int TOP = 0;
const int BOTTOM = 1;
const int LEFT = 2;
const int RIGHT = 3;
} // namespace border

enum class cell_type {

    FLUID,        // 0
    INFLOW,       // 1
    INLET_A,      // 4, considered an inflow cell
    INLET_B,      // 5, considered an inflow cell
    OUTFLOW,      // 2
    FIXED_WALL,   // 3
    FLUID_BUFFER, // 6
    CONVERSION,   // 7, considered a fluid cell
    MOVING_WALL,  // 8
    DEFAULT       // 9
};

namespace boundary_ids {
const int fluid_cell_id = 0;
const int inflow_cell_id = 1;
const int outflow_cell_id = 2;
const int fixed_wall_cell_3_id = 3;
const int inlet_a_cell_id = 4;
const int inlet_b_cell_id = 5;
const int fluidbuffer_cell_id = 6;
const int conversion_cell_id = 7;
const int moving_wall_cell_id = 8;
const int fixed_wall_cell_9_id = 9;
} // namespace boundary_ids
