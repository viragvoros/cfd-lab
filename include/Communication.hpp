#pragma once

#include <mpi.h>

#include "Domain.hpp"
#include "Datastructures.hpp"

/**
  * @brief Class of container and modifier for the communication
 *
 */
class Communication {
  public:
    Communication() = default;

    /**
    * @brief Communicate between threads
    *
    * @param[in] arg from main.cpp
    *
    */
    static void communicate(Matrix<double> &A, Domain &domain, int &iproc, int &jproc);
};