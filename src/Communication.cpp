#include "Communication.hpp"

#include <mpi.h>

void Communication::communicate(Matrix<double> &A, Domain &domain, int &iproc, int &jproc) {
    int own_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &own_rank);

    double send_y[domain.size_y + 2];
    double receive_y[domain.size_y + 2];

    // TODO implement for other directions

    //----------------------------------------------------------------
    // SEND TO TOP, RECEIVE FROM BOTTOM
    //----------------------------------------------------------------
    // Copy data to send the TOP non-ghost layer of Matrix A
    for(int i = 0; i < domain.size_x + 2; i++){
        send_y[i] = A(i, domain.size_y);
    }
    // Find destination rank (TOP to current rank)
    int dest_rank = own_rank + iproc;
    // Find source rank (BOTTOM to current rank)
    int source_rank = own_rank - iproc;
    // Exchange data
    if ((own_rank + iproc) >= (iproc * jproc)){
        MPI_Sendrecv(&send_y, domain.size_y + 2, MPI_DOUBLE, MPI_PROC_NULL, 0,
            &receive_y, domain.size_y + 2, MPI_DOUBLE, source_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if (own_rank < iproc){
        MPI_Sendrecv(&send_y, domain.size_y + 2, MPI_DOUBLE, dest_rank, 0,
            &receive_y, domain.size_y + 2, MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Sendrecv(&send_y, domain.size_y + 2, MPI_DOUBLE, dest_rank, 0,
            &receive_y, domain.size_y + 2, MPI_DOUBLE, source_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // Copy data from revieve into Matrix A BOTTOM ghost layer
    for(int i = 0; i < domain.size_x + 2; i++){
        A(i, 0) = receive_y[i];
    }
}

