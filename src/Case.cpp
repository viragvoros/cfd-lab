#include "Case.hpp"
#include "Enums.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

namespace filesystem = std::filesystem;

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>

Case::Case(std::string file_name, int argn, char **args) {
    // Read input parameters
    const int MAX_LINE_LENGTH = 1024;
    std::ifstream file(file_name);
    double nu;      /* viscosity   */
    double UI;      /* velocity x-direction */
    double VI;      /* velocity y-direction */
    double PI;      /* pressure */
    double GX;      /* gravitation x-direction */
    double GY;      /* gravitation y-direction */
    double xlength; /* length of the domain x-dir.*/
    double ylength; /* length of the domain y-dir.*/
    double dt;      /* time step */
    int imax;       /* number of cells x-direction*/
    int jmax;       /* number of cells y-direction*/
    double gamma;   /* uppwind differencing factor*/
    double omg;     /* relaxation factor */
    double tau;     /* safety factor for time step*/
    int itermax;    /* max. number of iterations for pressure per time step */
    double eps;     /* accuracy bound for pressure*/
    int iproc;
    int jproc;
    double UIN;       /* inlet velocity x-direction */
    double VIN;       /* inlet velocity y-direction */
    int num_of_walls; /* number of walls */
    double TI;        /* initial temperature */
    double TIN;       /* inlet temperature */
    double beta;      /* thermal expansion coefficient */
    double alpha;     /* thermal diffusivity */

    if (file.is_open()) {

        std::string var;
        while (!file.eof() && file.good()) {
            file >> var;
            if (var[0] == '#') { /* ignore comment line*/
                file.ignore(MAX_LINE_LENGTH, '\n');
            } else {
                if (var == "geo_file") file >> _geom_name;
                if (var == "xlength") file >> xlength;
                if (var == "ylength") file >> ylength;
                if (var == "nu") file >> nu;
                if (var == "t_end") file >> _t_end;
                if (var == "dt") file >> dt;
                if (var == "omg") file >> omg;
                if (var == "eps") file >> eps;
                if (var == "tau") file >> tau;
                if (var == "gamma") file >> gamma;
                if (var == "dt_value") file >> _output_freq;
                if (var == "UI") file >> UI;
                if (var == "VI") file >> VI;
                if (var == "GX") file >> GX;
                if (var == "GY") file >> GY;
                if (var == "PI") file >> PI;
                if (var == "itermax") file >> itermax;
                if (var == "imax") file >> imax;
                if (var == "jmax") file >> jmax;
                if (var == "iproc") file >> iproc;
                if (var == "jproc") file >> jproc;
                if (var == "UIN") file >> UIN;
                if (var == "VIN") file >> VIN;
                if (var == "num_of_walls") file >> num_of_walls;
                if (var == "energy_eq") file >> energy_eq;
                if (var == "TI") file >> TI;
                if (var == "TIN") file >> TIN;
                if (var == "beta") file >> beta;
                if (var == "alpha") file >> alpha;
                if (var == "wall_temp_3") file >> wall_temp_3;
                if (var == "wall_temp_4") file >> wall_temp_4;
                if (var == "wall_temp_5") file >> wall_temp_5;
                if (var == "wall_vel_3") file >> wall_vel_3;
                if (var == "wall_vel_4") file >> wall_vel_4;
                if (var == "wall_vel_5") file >> wall_vel_5;
            }
        }
    }
    file.close();

    std::map<int, double> wall_vel;
    std::map<int, double> wall_temp;
    if (_geom_name.compare("NONE") == 0) {
        wall_vel.insert(std::pair<int, double>(LidDrivenCavity::moving_wall_id, LidDrivenCavity::wall_velocity));
    } else {
        // IDEA: construct maps of ids and velocities/temperatures for constructors of boundaries
        wall_vel[boundary_ids::fixed_wall_cell_3_id] = wall_vel_3;
        wall_vel[boundary_ids::fixed_wall_cell_4_id] = wall_vel_4;
        wall_vel[boundary_ids::fixed_wall_cell_5_id] = wall_vel_5;
        wall_temp[boundary_ids::fixed_wall_cell_3_id] = wall_temp_3;
        wall_temp[boundary_ids::fixed_wall_cell_4_id] = wall_temp_4;
        wall_temp[boundary_ids::fixed_wall_cell_5_id] = wall_temp_5;
    }

    // Set file names for geometry file and output directory
    set_file_names(file_name);

    // Build up the domain
    Domain domain;
    domain.dx = xlength / (double)imax;
    domain.dy = ylength / (double)jmax;
    domain.domain_size_x = imax;
    domain.domain_size_y = jmax;

    build_domain(domain, imax, jmax);

    _grid = Grid(_geom_name, domain);
    _field = Fields(nu, dt, tau, alpha, beta, _grid.fluid_cells(), _grid.domain().size_x, _grid.domain().size_y, UI, VI,
                    PI, TI, energy_eq);

    _discretization = Discretization(domain.dx, domain.dy, gamma);
    _pressure_solver = std::make_unique<SOR>(omg);
    _max_iter = itermax;
    _tolerance = eps;

    // TODO fix outlet Construct boundaries
    if (not _grid.moving_wall_cells().empty()) {
        _boundaries.push_back(
            std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), LidDrivenCavity::wall_velocity));
        std::cout << _boundaries.size() << std::endl;
    }
    if (not _grid.fixed_wall_cells_3().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells_3(),
                                                                  wall_temp[boundary_ids::fixed_wall_cell_3_id]));
    }
    if (not _grid.fixed_wall_cells_4().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells_4(),
                                                                  wall_temp[boundary_ids::fixed_wall_cell_4_id]));
    }
    if (not _grid.fixed_wall_cells_5().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells_5(),
                                                                  wall_temp[boundary_ids::fixed_wall_cell_5_id]));
    }
    /*
    if (not _grid.fixed_wall_cells_6().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells_6()));
    }
    */
    if (not _grid.inflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<InFlowBoundary>(_grid.inflow_cells(), UIN, TIN));
    }
    if (not _grid.outflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<OutFlowBoundary>(_grid.outflow_cells(), PI));
    }

    // std::cout << _boundaries.size() << std::endl; // is this supposed to be only 3? How is this organised?
}

void Case::set_file_names(std::string file_name) {
    std::string temp_dir;
    bool case_name_flag = true;
    bool prefix_flag = false;

    for (int i = file_name.size() - 1; i > -1; --i) {
        if (file_name[i] == '/') {
            case_name_flag = false;
            prefix_flag = true;
        }
        if (case_name_flag) {
            _case_name.push_back(file_name[i]);
        }
        if (prefix_flag) {
            _prefix.push_back(file_name[i]);
        }
    }

    for (int i = file_name.size() - _case_name.size() - 1; i > -1; --i) {
        temp_dir.push_back(file_name[i]);
    }

    std::reverse(_case_name.begin(), _case_name.end());
    std::reverse(_prefix.begin(), _prefix.end());
    std::reverse(temp_dir.begin(), temp_dir.end());

    _case_name.erase(_case_name.size() - 4);
    _dict_name = temp_dir;
    _dict_name.append(_case_name);
    _dict_name.append("_Output");

    if (_geom_name.compare("NONE") != 0) {
        _geom_name = _prefix + _geom_name;
    }

    // Create output directory
    filesystem::path folder(_dict_name);
    try {
        filesystem::create_directory(folder);
    } catch (const std::exception &e) {
        std::cerr << "Output directory could not be created." << std::endl;
        std::cerr << "Make sure that you have write permissions to the "
                     "corresponding location"
                  << std::endl;
    }
}

/**
 * This function is the main simulation loop. In the simulation loop, following steps are required
 * - Calculate and apply boundary conditions for all the boundaries in _boundaries container
 *   using apply() member function of Boundary class
 * - Calculate fluxes (F and G) using calculate_fluxes() member function of Fields class.
 *   Flux consists of diffusion and convection part, which are located in Discretization class
 * - Calculate right-hand-side of PPE using calculate_rs() member function of Fields class
 * - Iterate the pressure poisson equation until the residual becomes smaller than the desired tolerance
 *   or the maximum number of the iterations are performed using solve() member function of PressureSolver class
 * - Calculate the velocities u and v using calculate_velocities() member function of Fields class
 * - Calculat the maximal timestep size for the next iteration using calculate_dt() member function of Fields class
 * - Write vtk files using output_vtk() function
 *
 * Please note that some classes such as PressureSolver, Boundary are abstract classes which means they only provide the
 * interface. No member functions should be defined in abstract classes. You need to define functions in inherited
 * classes such as MovingWallBoundary class.
 *
 * For information about the classes and functions, you can check the header files.
 */
void Case::simulate() {

    std::cout << "Entering simulate" << std::endl;
    double t = 0.0;
    double dt = _field.dt();
    int timestep = 0;
    // Changed to integer to use as input parameter for vtk file generation
    int output_counter = 0;

    while (t <= _t_end) {
        // Calculate optimum timestep
        dt = _field.calculate_dt(_grid);

        std::cout << dt << std::endl;
        std::cout << energy_eq << std::endl;
        // Application of boundary conditions
        for (auto &boundary : _boundaries) {
            boundary->apply(_field);
            // std::cout << "Entering apply boundaries" << std::endl;
        }

        _field.calculate_temperature(_grid);
        _field.calculate_fluxes(_grid);
        _field.calculate_rs(_grid);

        int nb_iter = 0;
        while (nb_iter <= _max_iter) {
            double res = _pressure_solver->solve(_field, _grid, _boundaries);
            if (res <= _tolerance) {
                // std::cout << res << std::endl;
                break;
            }
            nb_iter++;
        }

        if (nb_iter == _max_iter + 1) {
            std::cout << "WARNING: SOR SOLVER DID NOT CONVERGE IN TIMESTEP " << output_counter + 1 << "\n"
                      << "OBTAINED RESULTS MIGHT BE ERRONOUS. \n";
        }

        _field.calculate_velocities(_grid);

        // --------------DEBUG: Printig fields to console ----------------

        //     std::cout << "---------------------- u field ------------------------------" << std::endl;
        //  for (int jx = 0; jx < 22; jx++ ){
        //              for (int ix = 0; ix < 73; ix++) {
        //                  std::cout << _field.u(ix, jx) << " " ;
        //                 }
        //             std::cout << "\n";
        //          }

        // std::cout << "---------------------- v field ------------------------------" << std::endl;
        //  for (int jx = 0; jx < 22; jx++ ){
        //              for (int ix = 0; ix < 73; ix++) {
        //                  std::cout << _field.v(ix, jx) << " " ;
        //                 }
        //             std::cout << "\n";
        //          }

        output_counter++;
        t = t + dt;

        if (output_counter == 20 || output_counter % 100 == 0) {
            for (auto &boundary : _boundaries) {
                boundary->apply(_field);
                // std::cout << "Entering apply boundaries" << std::endl;
            }
            output_vtk(output_counter);
        }

        //     // Additional application of boundary conditions as mentioned in tutorial.
        // for (auto &boundary : _boundaries) {
        //     boundary->apply(_field);
        // }
    }

    output_vtk(output_counter);
}

void Case::output_vtk(int file_number) {
    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    double dx = _grid.dx();
    double dy = _grid.dy();

    double x = _grid.domain().imin * dx;
    double y = _grid.domain().jmin * dy;

    { y += dy; }
    { x += dx; }

    double z = 0;
    for (int col = 0; col < _grid.domain().size_y + 1; col++) {
        x = _grid.domain().imin * dx;
        { x += dx; }
        for (int row = 0; row < _grid.domain().size_x + 1; row++) {

            // if(_grid.get_geometry_data().at(x).at(y) == 0){ // maybe TODO .vtk does not work with this
            points->InsertNextPoint(x, y, z);
            //}

            x += dx;
        }
        y += dy;
    }

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(_grid.domain().size_x + 1, _grid.domain().size_y + 1, 1);
    structuredGrid->SetPoints(points);

    // Pressure Array
    vtkDoubleArray *Pressure = vtkDoubleArray::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);

    // Temperature Array
    vtkDoubleArray *Temperature = vtkDoubleArray::New();
    Temperature->SetName("temperature");
    Temperature->SetNumberOfComponents(1);

    // Velocity Array
    vtkDoubleArray *Velocity = vtkDoubleArray::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Print pressure and temperature from bottom to top
    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            if (_geom_name.compare("NONE") == 0) {
                double pressure = _field.p(i, j);
                Pressure->InsertNextTuple(&pressure);

                double temperature = _field.t(i, j);
                Temperature->InsertNextTuple(&temperature);
            } else if (_grid.get_geometry_data().at(i).at(j) == 0 || _grid.get_geometry_data().at(i).at(j) == 1 ||
                       _grid.get_geometry_data().at(i).at(j) == 2) {
                double pressure = _field.p(i, j);
                Pressure->InsertNextTuple(&pressure);

                double temperature = _field.t(i, j);
                Temperature->InsertNextTuple(&temperature);
            } else {
                double pressure = 0;
                Pressure->InsertNextTuple(&pressure);

                if (_grid.get_geometry_data().at(i).at(j) == 3) {
                    double temperature = wall_temp_3;
                    Temperature->InsertNextTuple(&temperature);
                } else if (_grid.get_geometry_data().at(i).at(j) == 4) {
                    double temperature = wall_temp_4;
                    Temperature->InsertNextTuple(&temperature);
                } else if (_grid.get_geometry_data().at(i).at(j) == 5) {
                    double temperature = wall_temp_5;
                    Temperature->InsertNextTuple(&temperature);
                } else {
                    // TODO ID 7/8/9 how to be represented?
                    double temperature = 0;
                    Temperature->InsertNextTuple(&temperature);
                }
            }
        }
    }

    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0

    // Print Velocity from bottom to top
    for (int j = 0; j < _grid.domain().size_y + 1; j++) {
        for (int i = 0; i < _grid.domain().size_x + 1; i++) {
            if (_geom_name.compare("NONE") == 0) {
                vel[0] = (_field.u(i, j) + _field.u(i, j + 1)) * 0.5;
                vel[1] = (_field.v(i, j) + _field.v(i + 1, j)) * 0.5;
                Velocity->InsertNextTuple(vel);
            } else if (_grid.get_geometry_data().at(i).at(j) == 0 || _grid.get_geometry_data().at(i).at(j) == 1 ||
                       _grid.get_geometry_data().at(i).at(j) == 2) {
                vel[0] = (_field.u(i, j) + _field.u(i, j + 1)) * 0.5;
                vel[1] = (_field.v(i, j) + _field.v(i + 1, j)) * 0.5;
                Velocity->InsertNextTuple(vel);
            } else {
                vel[0] = 0;
                vel[1] = 0;
                Velocity->InsertNextTuple(vel);
            }
        }
    }

    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);

    // Add Temperature to Structured Grid
    structuredGrid->GetCellData()->AddArray(Temperature);

    // Add Velocity to Structured Grid
    structuredGrid->GetPointData()->AddArray(Velocity);

    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname = _dict_name + '/' + _case_name + "_" + std::to_string(file_number) + ".vtk";

    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}

void Case::build_domain(Domain &domain, int imax_domain, int jmax_domain) {
    domain.imin = 0;
    domain.jmin = 0;
    domain.imax = imax_domain + 2;
    domain.jmax = jmax_domain + 2;
    domain.size_x = imax_domain;
    domain.size_y = jmax_domain;
}
