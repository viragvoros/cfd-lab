#include <chrono>
#include <iostream>
#include <mpi.h>
#include <string>

#include "Case.hpp"

int main(int argn, char **args) {

    if (argn > 1) {
        std::string file_name{args[1]};
        auto start = std::chrono::high_resolution_clock::now();
        MPI_Init(&argn, &args);
        Case problem(file_name, argn, args);
        problem.simulate();
        MPI_Finalize();
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "\nTotal computation time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s" << std::endl;
    } else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
}
