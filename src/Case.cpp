#include "Case.hpp"
#include "Enums.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <mpi.h>
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
    double nu;                  /* viscosity   */
    double UI;                  /* velocity x-direction */
    double VI;                  /* velocity y-direction */
    double PI;                  /* pressure */
    double CAI{0};              /* initial concentration of component A */
    double CBI{0};              /* initial concentration of component B */
    double CCI{0};              /* initial concentration of component C */
    double CAIN{0};             /* inlet concentration of component A */
    double CBIN{0};             /* inlet concentration of component B */
    double GX;                  /* gravitation x-direction */
    double GY;                  /* gravitation y-direction */
    double xlength;             /* length of the domain x-dir.*/
    double ylength;             /* length of the domain y-dir.*/
    double dt;                  /* time step */
    int imax;                   /* number of cells x-direction*/
    int jmax;                   /* number of cells y-direction*/
    double gamma;               /* uppwind differencing factor*/
    double omg;                 /* relaxation factor */
    double tau;                 /* safety factor for time step*/
    int itermax;                /* max. number of iterations for pressure per time step */
    double eps;                 /* accuracy bound for pressure*/
    double UIN;                 /* inlet velocity x-direction */
    double VIN;                 /* inlet velocity y-direction */
    int num_of_walls;           /* number of walls */
    double TI;                  /* initial temperature */
    double TIN;                 /* inlet temperature */
    double beta;                /* thermal expansion coefficient */
    double alpha;               /* thermal diffusivity */
    double diffusivity;         /* diffusivity */
    double rate_const;          /* reaction rate constant */
    double order_a;             /* order of reaction with regard to A */
    double order_b;             /* order of reaction with regard to B */
    double kappa{1};            /* Thermal conductivity */
    double wall_heatflux_3{0};  /* heat flux specified at boundary */
    double wall_heatflux_9{0};  /* heat flux specified at boundary */
    double pre_exp_factor;      /* pre-exponential factor */
    double act_energy;          /* activation energy */
    double react_temp_increase; /* increase in temperature from reaction heat */

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
                if (var == "CAI") file >> CAI;
                if (var == "CBI") file >> CBI;
                if (var == "CCI") file >> CCI;
                if (var == "CAIN") file >> CAIN;
                if (var == "CBIN") file >> CBIN;
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
                if (var == "diffusivity") file >> diffusivity;
                if (var == "rate_const") file >> rate_const;
                if (var == "order_a") file >> order_a;
                if (var == "order_b") file >> order_b;
                if (var == "wall_temp_3") file >> wall_temp_3;
                if (var == "wall_temp_9") file >> wall_temp_9;
                if (var == "wall_vel_3") file >> wall_vel_3;
                if (var == "kappa") file >> kappa;
                if (var == "wall_heatflux_3") file >> wall_heatflux_3;
                if (var == "wall_heatflux_9") file >> wall_heatflux_9;
                if (var == "pre_exp_factor") file >> pre_exp_factor;
                if (var == "act_energy") file >> act_energy;
                if (var == "react_temp_increase") file >> react_temp_increase;
                if (var == "ref_factor") file >> ref_factor;
            }
        }
    }
    file.close();

    // Redefine imax and jmax for refinement
    imax = imax * ref_factor;
    jmax = jmax * ref_factor;

    int error_size;
    MPI_Comm_size(MPI_COMM_WORLD, &error_size);
    int error_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &error_rank);

    if (error_rank == 0 && error_size != jproc * iproc) {
        std::cout << "-------------------------- ERROR: INVALID INPUT --------------------------" << std::endl;
        std::cout << "Current value of jproc * jproc: " << jproc * iproc
                  << ". Your given input for number of processes is " << error_size
                  << ". Change it so it matches jproc * jproc." << std::endl;
        std::cout << "Example usage: mpirun -np (iproc*jproc) /path/to/fluidchen /path/to/input_data.dat" << std::endl;
        exit(1);
    }

    std::map<int, double> wall_vel;
    std::map<int, double> wall_temp;
    std::map<int, double> wall_heatflux;
    if (_geom_name.compare("NONE") == 0) {
        wall_vel.insert(std::pair<int, double>(LidDrivenCavity::moving_wall_id, LidDrivenCavity::wall_velocity));
    } else {
        // Construct maps of ids and velocities/temperatures for constructors of boundaries
        wall_vel[boundary_ids::fixed_wall_cell_3_id] = wall_vel_3;
        wall_temp[boundary_ids::fixed_wall_cell_3_id] = wall_temp_3;
        wall_temp[boundary_ids::fixed_wall_cell_9_id] = wall_temp_9;
        wall_heatflux[boundary_ids::fixed_wall_cell_3_id] = wall_heatflux_3;
        wall_heatflux[boundary_ids::fixed_wall_cell_9_id] = wall_heatflux_9;
    }

    // Set file names for geometry file and output directory
    set_file_names(file_name);

    // Build up the domain
    domain.dx = xlength / (double)imax; // physical length of one cell
    domain.dy = ylength / (double)jmax;
    domain.domain_size_x = imax; // global domain size x (amount of cells in x direction) without the ghost cell frame
    domain.domain_size_y = jmax; // global domain size y (amount of cells in y direction) without the ghost cell frame
    domain.ref_factor = ref_factor;

    build_domain(domain, imax, jmax);

    _grid = Grid(_geom_name, domain);
    _field = Fields(nu, dt, tau, alpha, beta, diffusivity, rate_const, order_a, order_b, pre_exp_factor, act_energy,
                    react_temp_increase, _grid.fluid_cells(), _grid.domain().size_x, _grid.domain().size_y, UI, VI, PI,
                    TI, CAI, CBI, CBI, energy_eq, GX, GY);

    _discretization = Discretization(domain.dx, domain.dy, gamma);
    _pressure_solver = std::make_unique<SOR>(omg);
    _max_iter = itermax;
    _tolerance = eps;

    // Construct boundaries
    if (not _grid.moving_wall_cells().empty()) {
        _boundaries.push_back(
            std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), LidDrivenCavity::wall_velocity));
    }
    if (not _grid.fixed_wall_cells().empty()) {
        _boundaries.push_back(
            std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells(), wall_temp, kappa, wall_heatflux));
    }
    if (not _grid.inflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<InFlowBoundary>(_grid.inflow_cells(), UIN, TIN, CAIN, CBIN));
    }
    if (not _grid.outflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<OutFlowBoundary>(_grid.outflow_cells(), PI));
    }
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
    _dict_name.append("_Output_");
    std::string string_iproc = std::to_string(iproc);
    std::string string_jproc = std::to_string(jproc);
    _dict_name.append(string_iproc);
    _dict_name.append("x");
    _dict_name.append(string_jproc);

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

    double t = 0.0;
    double dt = _field.dt();

    // Used as input parameter for vtk file generation
    int output_counter = 0;

    // Getting rank for vtk output name
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0) {
        std::cout << "----------------Start of Simulation--------------\n"
                  << "Printed simulation statistics are valid for the master rank only\n"
                  << std::endl;
    }

    // Calculate local maximum for velocities
    double local_max_v;
    double local_max_u;

    // Declare parameters for relative update calculation
    double previous_mean_u;
    double previous_mean_v;
    double mean_p;
    double previous_mean_t;
    double u_rel_update;
    double v_rel_update;
    double t_rel_update;

    // Counter for SOR fails
    int SOR_fail_counter{0};

    while (t <= _t_end) {
        // Store old mean values for later relative update calculation
        previous_mean_u = _field.u_avg();
        previous_mean_v = _field.v_avg();
        previous_mean_t = _field.t_avg();

        // Application of boundary conditions
        for (auto &boundary : _boundaries) {
            boundary->apply(_field);
        }

        _field.calculate_temperature(_grid);
        _field.calculate_concentrations(_grid);
        _field.react(_grid);

        _communication.communicate(_field.t_matrix(), domain, iproc, jproc);
        _communication.communicate(_field.ca_matrix(), domain, iproc, jproc);
        _communication.communicate(_field.cb_matrix(), domain, iproc, jproc);
        _communication.communicate(_field.cc_matrix(), domain, iproc, jproc);

        _field.calculate_fluxes(_grid);

        // Application of boundary conditions
        for (auto &boundary : _boundaries) {
            boundary->apply(_field);
        }

        // Communicate fluxes
        _communication.communicate(_field.f_matrix(), domain, iproc, jproc);
        _communication.communicate(_field.g_matrix(), domain, iproc, jproc);

        _field.calculate_rs(_grid);

        int nb_iter = 0;
        while (nb_iter <= _max_iter) {
            double res = _pressure_solver->solve(_field, _grid, _boundaries);
            _communication.communicate(_field.p_matrix(), domain, iproc, jproc);
            if (res <= _tolerance) {
                mean_p = res; // pressure relative update
                break;
            }
            nb_iter++;
        }

        if (nb_iter == _max_iter + 1) {
            if (my_rank == 0) {
                std::cout << "WARNING: SOR SOLVER DID NOT CONVERGE IN TIMESTEP " << output_counter + 1 << "\t"
                          << "OBTAINED RESULTS MIGHT BE ERRONOUS. \n";
            }
            SOR_fail_counter++;
        }

        _field.calculate_velocities(_grid);

        _communication.communicate(_field.u_matrix(), domain, iproc, jproc);
        _communication.communicate(_field.v_matrix(), domain, iproc, jproc);

        local_max_u = _field.find_max(_field.u_matrix());
        local_max_v = _field.find_max(_field.v_matrix());

        MPI_Allreduce(MPI_IN_PLACE, &local_max_u, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &local_max_v, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // Calculate optimum timestep
        dt = _field.calculate_dt(_grid, local_max_u, local_max_v);

        output_counter++;
        t = t + dt;

        if (output_counter % 10 == 0) {
            if (my_rank == 0) {

                u_rel_update = std::abs(1 - previous_mean_u / _field.u_avg());
                v_rel_update = std::abs(1 - previous_mean_v / _field.v_avg());
                t_rel_update = std::abs(1 - previous_mean_t / _field.t_avg());

                std::cout << std::fixed;
                std::cout << std::setprecision(3);
                std::cout << "Time: " << t << "\t";
                std::cout << std::setprecision(5) << "dt: " << dt << "\t"
                          << "SOR-Iter: " << nb_iter << "\t";
                std::cout << std::setprecision(3) << "U-Rel-Update: " << std::scientific << u_rel_update << "\t\t"
                          << "V-Rel-Update: " << v_rel_update << "\t\t"
                          << "P-Rel-Update: " << mean_p << "\t\t"
                          << "T-Rel-Update: " << t_rel_update << std::endl;
            }
            output_vtk(output_counter, my_rank);
        }
    }

    output_vtk(output_counter, my_rank);
    if (my_rank == 0) {
        std::cout << "\nEnd of Simulation \n\nSimulation Report: "
                  << "End time: " << t << "\t || SOR Solver failed " << SOR_fail_counter
                  << " times. Check previous terminal information to find the corresponding timesteps." << std::endl;
    }
}

void Case::output_vtk(int file_number, int my_rank) {
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
            points->InsertNextPoint(x, y, z);
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

    // Concentration A Array
    vtkDoubleArray *CA = vtkDoubleArray::New();
    CA->SetName("concentration_a");
    CA->SetNumberOfComponents(1);

    // Concentration B Array
    vtkDoubleArray *CB = vtkDoubleArray::New();
    CB->SetName("concentration_b");
    CB->SetNumberOfComponents(1);

    // Concentration C Array
    vtkDoubleArray *CC = vtkDoubleArray::New();
    CC->SetName("concentration_c");
    CC->SetNumberOfComponents(1);

    // Print pressure and temperature from bottom to top
    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            if (_geom_name.compare("NONE") == 0) {
                double pressure = _field.p(i, j);
                Pressure->InsertNextTuple(&pressure);

                double concentration_a = _field.ca(i, j);
                double concentration_b = _field.cb(i, j);
                double concentration_c = _field.cc(i, j);
                CA->InsertNextTuple(&concentration_a);
                CB->InsertNextTuple(&concentration_b);
                CC->InsertNextTuple(&concentration_c);

            } else if (_grid.get_geometry_data().at(i).at(j) == 0 || _grid.get_geometry_data().at(i).at(j) == 2 ||
                       _grid.get_geometry_data().at(i).at(j) == 7) {
                double pressure = _field.p(i, j);
                Pressure->InsertNextTuple(&pressure);

                double concentration_a = _field.ca(i, j);
                double concentration_b = _field.cb(i, j);
                double concentration_c = _field.cc(i, j);
                CA->InsertNextTuple(&concentration_a);
                CB->InsertNextTuple(&concentration_b);
                CC->InsertNextTuple(&concentration_c);
            } else {
                double pressure = 0;
                Pressure->InsertNextTuple(&pressure);

                double concentration_a = 0;
                double concentration_b = 0;
                double concentration_c = 0;
                CA->InsertNextTuple(&concentration_a);
                CB->InsertNextTuple(&concentration_b);
                CC->InsertNextTuple(&concentration_c);
            }

            if (energy_eq.compare("NONE") != 0) {
                if (_geom_name.compare("NONE") == 0) {
                    double temperature = _field.t(i, j);
                    Temperature->InsertNextTuple(&temperature);
                } else if (_grid.get_geometry_data().at(i).at(j) == 0 || _grid.get_geometry_data().at(i).at(j) == 2 ||
                           _grid.get_geometry_data().at(i).at(j) == 7) {
                    double temperature = _field.t(i, j);
                    Temperature->InsertNextTuple(&temperature);
                } else {
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
            if (_geom_name.compare("NONE") == 0) { // Lid driven cavity
                vel[0] = (_field.u(i, j) + _field.u(i, j + 1)) * 0.5;
                vel[1] = (_field.v(i, j) + _field.v(i + 1, j)) * 0.5;
                Velocity->InsertNextTuple(vel);

            } else if (_grid.get_geometry_data().at(i).at(j) == 0 || _grid.get_geometry_data().at(i).at(j) == 2 ||
                       _grid.get_geometry_data().at(i).at(j) == 7 || j == 0) {
                vel[0] = (_field.u(i, j) + _field.u(i, j + 1)) * 0.5;
                vel[1] = (_field.v(i, j) + _field.v(i + 1, j)) * 0.5;
                Velocity->InsertNextTuple(vel);

            } else { // Obstacles become zero
                vel[0] = 0;
                vel[1] = 0;
                // Velocity->InsertNextTuple(vel);
                structuredGrid->BlankPoint(Velocity->InsertNextTuple(vel));
            }
        }
    }

    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);

    // Add Temperature to Structured Grid
    if (energy_eq.compare("NONE") != 0) {
        structuredGrid->GetCellData()->AddArray(Temperature);
    }

    // Add Concentration to Structured Grid
    structuredGrid->GetCellData()->AddArray(CA);
    // Add Concentration to Structured Grid
    structuredGrid->GetCellData()->AddArray(CB);
    // Add Concentration to Structured Grid
    structuredGrid->GetCellData()->AddArray(CC);

    // Add Velocity to Structured Grid
    structuredGrid->GetPointData()->AddArray(Velocity);

    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname =
        _dict_name + '/' + _case_name + "_" + std::to_string(my_rank) + "." + std::to_string(file_number) + ".vtk";

    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}

void Case::build_domain(Domain &domain, int imax, int jmax) {
    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        domain.size_x = floor(imax / iproc);
        domain.size_y = floor(jmax / jproc);
        domain.imin = 0;
        domain.jmin = 0;
        domain.imax = domain.imin + domain.size_x + 2;
        domain.jmax = domain.jmin + domain.size_y + 2;

        int imin_loc = 1;
        int jmin_loc = 1;
        int imax_loc;
        int jmax_loc = jmin_loc + domain.size_y - 1;

        int send_size_x;
        int send_size_y;
        int send_imin_loc;
        int send_jmin_loc;
        int send_imax_loc;
        int send_jmax_loc;

        // Send data of domain
        for (int k = 1; k < jproc * iproc; k++) {
            send_size_x = domain.size_x;
            send_size_y = domain.size_y;

            imin_loc += domain.size_x;
            imax_loc = imin_loc + domain.size_x - 1;

            if ((k % iproc) == 0) {
                imin_loc = 1;
                imax_loc = imin_loc + domain.size_x - 1; // Reassigning because of imin_loc
                jmin_loc += domain.size_y;
                jmax_loc = jmin_loc + domain.size_y - 1;
            }

            // On TOP edge in case of uneven j division
            if ((k + iproc) >= (iproc * jproc) && (jmax % jproc) != 0) {
                send_size_y = domain.size_y + jmax % jproc;
                jmax_loc = jmin_loc + send_size_y - 1;
            }
            // On RIGHT edge in case of uneven i division
            if (((k + 1) % iproc) == 0 && (imax % iproc) != 0) {
                send_size_x = domain.size_x + imax % iproc;
                imax_loc = imin_loc + send_size_x - 1;
            }

            send_imin_loc = imin_loc - 1;
            send_jmin_loc = jmin_loc - 1;
            send_imax_loc = imax_loc + 2;
            send_jmax_loc = jmax_loc + 2;

            MPI_Send(&send_size_x, 1, MPI_INT, k, 1, MPI_COMM_WORLD);
            MPI_Send(&send_size_y, 1, MPI_INT, k, 2, MPI_COMM_WORLD);
            MPI_Send(&send_imin_loc, 1, MPI_INT, k, 3, MPI_COMM_WORLD);
            MPI_Send(&send_jmin_loc, 1, MPI_INT, k, 4, MPI_COMM_WORLD);
            MPI_Send(&send_imax_loc, 1, MPI_INT, k, 5, MPI_COMM_WORLD);
            MPI_Send(&send_jmax_loc, 1, MPI_INT, k, 6, MPI_COMM_WORLD);
        }

    } else {
        MPI_Recv(&domain.size_x, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&domain.size_y, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&domain.imin, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&domain.jmin, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&domain.imax, 1, MPI_INT, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&domain.jmax, 1, MPI_INT, 0, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}
