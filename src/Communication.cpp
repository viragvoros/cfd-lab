#include "Communication.hpp"

#include <mpi.h>

#include <iostream>

void Communication::communicate(Matrix<double> &A, Domain &domain, int &iproc, int &jproc) {
    int own_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &own_rank);

    if (jproc != 1) {
        //----------------------------------------------------------------
        // From BOTTOM rank to TOP rank
        //----------------------------------------------------------------
        // Initialize buffers to send and receive data
        double send_tb_x[domain.size_x + 2];
        double receive_tb_x[domain.size_x + 2];

        // Store topmost row of domain in send-buffer
        for (int i = 0; i < domain.size_x + 2; i++) {
            send_tb_x[i] = A(i, domain.size_y);
        }

        int dest_tb_rank, source_tb_rank;

        // Find destination rank (TOP to current rank)
        if (own_rank < iproc * jproc - iproc) { // Rank < MaxRank  is not existent
            dest_tb_rank = own_rank + iproc;
        }

        // Find source rank (BOTTOM to current rank)
        if (own_rank >= iproc) { // Negative ranks do not exist
            source_tb_rank = own_rank - iproc;
        }

        // Exchange data
        if (own_rank < iproc) { // SEND
            MPI_Sendrecv(&send_tb_x, domain.size_x + 2, MPI_DOUBLE, dest_tb_rank, 0, &receive_tb_x, domain.size_x + 2,
                         MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else if ((own_rank + iproc) >= (iproc * jproc)) { // RECEIVE
            MPI_Sendrecv(&send_tb_x, domain.size_x + 2, MPI_DOUBLE, MPI_PROC_NULL, 0, &receive_tb_x, domain.size_x + 2,
                         MPI_DOUBLE, source_tb_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Update matrix with data from receive-buffer
            for (int i = 0; i < domain.size_x + 2; i++) {
                A(i, 0) = receive_tb_x[i];
            }

        } else { // MIDDLE: SEND and RECEIVE
            MPI_Sendrecv(&send_tb_x, domain.size_x + 2, MPI_DOUBLE, dest_tb_rank, 0, &receive_tb_x, domain.size_x + 2,
                         MPI_DOUBLE, source_tb_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Update matrix with data from receive-buffer
            for (int i = 0; i < domain.size_x + 2; i++) {
                A(i, 0) = receive_tb_x[i];
            }
        }

        //----------------------------------------------------------------
        // From TOP rank to BOTTOM rank
        //----------------------------------------------------------------
        // Initialize buffers to send and receive data
        double send_bt_x[domain.size_x + 2];
        double receive_bt_x[domain.size_x + 2];

        // Store bottommost row of domain in send-buffer
        for (int i = 0; i < domain.size_x + 2; i++) {
            send_bt_x[i] = A(i, 1);
        }

        int dest_bt_rank, source_bt_rank;

        // Find destination rank (BOTTOM to current rank)
        if (own_rank >= iproc) { // Rank = -1 is not existent
            dest_bt_rank = own_rank - iproc;
        }

        // Find source rank (TOP to current rank)
        if (own_rank < iproc * jproc - iproc) { // Rank < MaxRank  is not existent
            source_bt_rank = own_rank + iproc;
        }

        // Exchange data
        if ((own_rank + iproc) >= (iproc * jproc)) { // SEND
            MPI_Sendrecv(&send_bt_x, domain.size_x + 2, MPI_DOUBLE, dest_bt_rank, 0, &receive_bt_x, domain.size_x + 2,
                         MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else if (own_rank < iproc) { // RECEIVE
            MPI_Sendrecv(&send_bt_x, domain.size_x + 2, MPI_DOUBLE, MPI_PROC_NULL, 0, &receive_bt_x, domain.size_x + 2,
                         MPI_DOUBLE, source_bt_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Update matrix with data from receive-buffer
            for (int i = 0; i < domain.size_x + 2; i++) {
                A(i, domain.size_y + 1) = receive_bt_x[i];
            }

        } else { // MIDDLE: SEND and RECEIVE
            MPI_Sendrecv(&send_bt_x, domain.size_x + 2, MPI_DOUBLE, dest_bt_rank, 0, &receive_bt_x, domain.size_x + 2,
                         MPI_DOUBLE, source_bt_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Update matrix with data from receive-buffer
            for (int i = 0; i < domain.size_x + 2; i++) {
                A(i, domain.size_y + 1) = receive_bt_x[i];
            }
        }
    }

    if (iproc != 1) {
        //----------------------------------------------------------------
        // From RIGHT rank to LEFT rank
        //----------------------------------------------------------------
        // Initialize buffers to send and receive data
        double send_lr_y[domain.size_y + 2];
        double receive_lr_y[domain.size_y + 2];

        // Store leftmost column of domain in send-buffer
        for (int i = 0; i < domain.size_y + 2; i++) {
            send_lr_y[i] = A(1, i);
        }

        int dest_lr_rank, source_lr_rank;

        // Find destination rank (LEFT to current rank)
        if (own_rank % iproc != 0) { // Ranks on LEFT edge have no LEFT neighbor
            dest_lr_rank = own_rank - 1;
        }

        // Find sender rank (RIGHT to current rank)
        if ((own_rank + 1) % iproc != 0) { // Rank < MaxRank  is not existent
            source_lr_rank = own_rank + 1;
        }

        // Exchange data
        if ((own_rank + 1) % iproc == 0) { // SEND
            MPI_Sendrecv(&send_lr_y, domain.size_y + 2, MPI_DOUBLE, dest_lr_rank, 0, &receive_lr_y, domain.size_y + 2,
                         MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else if (own_rank % iproc == 0) { // RECEIVE
            MPI_Sendrecv(&send_lr_y, domain.size_y + 2, MPI_DOUBLE, MPI_PROC_NULL, 0, &receive_lr_y, domain.size_y + 2,
                         MPI_DOUBLE, source_lr_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Update matrix with data from receive-buffer
            for (int i = 0; i < domain.size_y + 2; i++) {
                A(domain.size_x + 1, i) = receive_lr_y[i];
            }

        } else { // MIDDLE: SEND and RECEIVE
            MPI_Sendrecv(&send_lr_y, domain.size_y + 2, MPI_DOUBLE, dest_lr_rank, 0, &receive_lr_y, domain.size_y + 2,
                         MPI_DOUBLE, source_lr_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Update matrix with data from receive-buffer
            for (int i = 0; i < domain.size_y + 2; i++) {
                A(domain.size_x + 1, i) = receive_lr_y[i];
            }
        }

        //----------------------------------------------------------------
        // From LEFT rank to RIGHT rank
        //----------------------------------------------------------------
        // Initialize buffers to send and receive data
        double send_rl_y[domain.size_y + 2];
        double receive_rl_y[domain.size_y + 2];

        // Store rightmost column of domain in send-buffer
        for (int i = 0; i < domain.size_y + 2; i++) {
            send_rl_y[i] = A(domain.size_x, i);
        }

        int dest_rl_rank, source_rl_rank;

        // Find destination rank (RIGHT to current rank)
        if (own_rank != iproc * jproc - 1) { // Rank < MaxRank  is not existent
            dest_rl_rank = own_rank + 1;
        }

        // Find source rank (LEFT to current rank)
        if (own_rank != 0) { // Ranks on LEFT edge have no LEFT neighbor
            source_rl_rank = own_rank - 1;
        }

        // Exchange data
        if ((own_rank + 1) % iproc == 0) { // RECEIVE
            MPI_Sendrecv(&send_rl_y, domain.size_y + 2, MPI_DOUBLE, MPI_PROC_NULL, 0, &receive_rl_y, domain.size_y + 2,
                         MPI_DOUBLE, source_rl_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Update matrix with data from receive-buffer
            for (int i = 0; i < domain.size_y + 2; i++) {
                A(0, i) = receive_rl_y[i];
            }

        } else if (own_rank % iproc == 0) { // SEND
            MPI_Sendrecv(&send_rl_y, domain.size_y + 2, MPI_DOUBLE, dest_rl_rank, 0, &receive_rl_y, domain.size_y + 2,
                         MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else { // MIDDLE: SEND and RECEIVE
            MPI_Sendrecv(&send_rl_y, domain.size_y + 2, MPI_DOUBLE, dest_rl_rank, 0, &receive_rl_y, domain.size_y + 2,
                         MPI_DOUBLE, source_rl_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Update matrix with data from receive-buffer
            for (int i = 0; i < domain.size_y + 2; i++) {
                A(0, i) = receive_rl_y[i];
            }
        }
    }
}
