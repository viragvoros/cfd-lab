#include "Communication.hpp"

#include <mpi.h>

#include <iostream>

void Communication::communicate(Matrix<double> &A, Domain &domain, int &iproc, int &jproc) {
    int own_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &own_rank);

if (jproc != 1) {
    //----------------------------------------------------------------
    // Sender column: TOP. Receiver column: BOTTOM
    //----------------------------------------------------------------
    // Copy data to send the Top non-ghost layer of Matrix A
    double send_tb_x[domain.size_x + 2];
    double receive_tb_x[domain.size_x + 2];

    // store topmost row not being ghost layer
    for (int i = 0; i < domain.size_x + 2; i++) {
        send_tb_x[i] = A(i, domain.size_y);
    }

    int dest_tb_rank, source_tb_rank;

    // Find source rank (BOTTOM to current rank)
    if (own_rank >= iproc) { // Rank = -1 is not existent
        source_tb_rank = own_rank - iproc;
    }

    // Find destination rank (TOP to current rank)
    if (own_rank < iproc * jproc - iproc) { // Rank < MaxRank  is not existent
        dest_tb_rank = own_rank + iproc;
    }


    // Exchange data
    if (own_rank < iproc) {// sender
        MPI_Sendrecv(&send_tb_x, domain.size_x + 2, MPI_DOUBLE, dest_tb_rank, 0,
                     &receive_tb_x, domain.size_x + 2, MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if ((own_rank + iproc) >= (iproc * jproc)) {  // receiver
        MPI_Sendrecv(&send_tb_x, domain.size_x + 2, MPI_DOUBLE, MPI_PROC_NULL, 0,
                     &receive_tb_x, domain.size_x + 2, MPI_DOUBLE, source_tb_rank, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

        // TODO: MAYBE REMOVE LATER
        for (int i = 0; i < domain.size_x + 2; i++) {
            A(i, 0) = receive_tb_x[i];
        }
        
    } else { // MIDDLE: SEND and RECEIVE
        MPI_Sendrecv(&send_tb_x, domain.size_x + 2, MPI_DOUBLE, dest_tb_rank, 0,
                     &receive_tb_x, domain.size_x + 2, MPI_DOUBLE, source_tb_rank, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

        // TODO: MAYBE REMOVE LATER
        for (int i = 0; i < domain.size_x + 2; i++) {
            A(i, 0) = receive_tb_x[i];
        }
    }


    /*
    // Copy data from receive into Matrix A BOTTOM ghost layer
    if ((own_rank+1)%iproc != 0){    //prevent non receivers to get bullshit
        for (int i = 0; i < domain.size_x + 2; i++) {
            A(i, 0) = receive_tb_x[i];
        }
    }
    */

    //----------------------------------------------------------------
    // Sender column: BOTTOM. Receiver column: TOP
    //----------------------------------------------------------------
    // Copy data to send the Bottom non-ghost layer of Matrix A
    double send_bt_x[domain.size_x + 2];
    double receive_bt_x[domain.size_x + 2];

    // store bottommost row not being ghost layer
    for (int i = 0; i < domain.size_x + 2; i++) {
        send_bt_x[i] = A(i, 1);
    }

    int dest_bt_rank, source_bt_rank;

    // Find destination rank (BOTTOM to current rank)
    if (own_rank >= iproc) { // Rank = -1 is not existent
        dest_bt_rank = own_rank - iproc;
    }

    // Find sender rank (TOP to current rank)
    if (own_rank < iproc * jproc - iproc) { // Rank < MaxRank  is not existent
        source_bt_rank = own_rank + iproc;
    }


    // Exchange data
    if ((own_rank + iproc) >= (iproc * jproc)) { // sender
        MPI_Sendrecv(&send_bt_x, domain.size_x + 2, MPI_DOUBLE, dest_bt_rank, 0,
                     &receive_bt_x, domain.size_x + 2, MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if (own_rank < iproc) { // receiver
        MPI_Sendrecv(&send_bt_x, domain.size_x + 2, MPI_DOUBLE, MPI_PROC_NULL, 0,
                     &receive_bt_x, domain.size_x + 2, MPI_DOUBLE, source_bt_rank, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);


        // TODO: MAYBE REMOVE LATER
        for (int i = 0; i < domain.size_x + 2; i++) {
            A(i, domain.size_y + 1) = receive_bt_x[i];
        }

    } else { // MIDDLE: SEND and RECEIVE
        MPI_Sendrecv(&send_bt_x, domain.size_x + 2, MPI_DOUBLE, dest_bt_rank, 0,
                     &receive_bt_x, domain.size_x + 2, MPI_DOUBLE, source_bt_rank, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

        // TODO: MAYBE REMOVE LATER
        for (int i = 0; i < domain.size_x + 2; i++) {
            A(i, domain.size_y + 1) = receive_bt_x[i];
        }
    }


    /*
    // Copy data from receive into Matrix A TOP ghost layer
    if ((own_rank+1)%iproc != 0){    //prevent non receivers to get bullshit
        for (int i = 0; i < domain.size_x + 2; i++) {
            A(i, domain.size_y + 1) = receive_bt_x[i];
        }
    }
    */
}
//std::cout << "hello" << std::endl;

if ( iproc != 1){//(own_rank%iproc != 0) || ((own_rank+1)%iproc != 0 )){
    //----------------------------------------------------------------
    // Sender column: LEFT. Receiver column: RIGHT
    //----------------------------------------------------------------
    // Copy data to send the LEFT non-ghost layer of Matrix A
    double send_lr_y[domain.size_y + 2];
    double receive_lr_y[domain.size_y + 2];

    // store column leftmost, not ghost into send
    for(int i = 0; i < domain.size_y + 2; i++){
        send_lr_y[i] = A(1, i);
    }

    int dest_lr_rank, source_lr_rank;

    // Find destination rank (LEFT to current rank)
    if (own_rank%iproc != 0) { // Rank = -1 is not existent
        dest_lr_rank = own_rank - 1;
    }

    // Find sender rank (RIGHT to current rank)
    if ((own_rank+1)%iproc != 0 ) { // Rank < MaxRank  is not existent
      source_lr_rank = own_rank + 1;
    }


    // Exchange data
    if ((own_rank+1)%iproc == 0){ // sender
        MPI_Sendrecv(&send_lr_y, domain.size_y + 2, MPI_DOUBLE, dest_lr_rank, 0,
                     &receive_lr_y, domain.size_y + 2, MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if (own_rank%iproc == 0){ // receiver
        MPI_Sendrecv(&send_lr_y, domain.size_y + 2, MPI_DOUBLE, MPI_PROC_NULL, 0,
                     &receive_lr_y, domain.size_y + 2 , MPI_DOUBLE, source_lr_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // TODO: MAYBE REMOVE LATER
        for(int i = 0; i < domain.size_y + 2 ; i++){
            A(domain.size_x + 1, i ) = receive_lr_y[i];
        }

    } else { //MIDDLE: SEND and RECEIVE
        MPI_Sendrecv(&send_lr_y, domain.size_y + 2 , MPI_DOUBLE, dest_lr_rank, 0,
                     &receive_lr_y, domain.size_y + 2 , MPI_DOUBLE, source_lr_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // TODO: MAYBE REMOVE LATER
        for(int i = 0; i < domain.size_y + 2 ; i++){
            A(domain.size_x + 1, i ) = receive_lr_y[i];
        }
    }


    /*
    // Copy data from receive into Matrix A Right ghost layer
    if ((own_rank+1)%iproc != 0){//prevent non receivers to get bullshit
        for(int i = 0; i < domain.size_y + 2 ; i++){
            A(domain.size_x + 1, i ) = receive_lr_y[i];
        }
    }
    */

                //if (own_rank == 0){
                //
                //    std::cout << "------------------Array --------------------\n" << std::endl;
                //
                //    for (int j = 0 ; j <= domain.size_y + 1; j++ ){
                ////        for (int i = 0; i <= domain.size_x + 1; i++) {
                ////
                ////            std::cout << A(i,j) << " " ;
                ////
                ////        }
                ////        std::cout << "\n";
                //        std::cout << receive_lr_y[j] << " ";
                //
                //
                //    }
                //    std::cout << "End of array\n" << std::endl;
                //}



    //----------------------------------------------------------------
    // FROM RIGHT TO LEFT
    //----------------------------------------------------------------
    // Copy data to send the RIGHT non-ghost layer of Matrix A
    double send_rl_y[domain.size_y + 2];
    double receive_rl_y[domain.size_y + 2];


        for (int i = 0; i < domain.size_y + 2; i++) {
            send_rl_y[i] = A(domain.size_x, i);
        }

    int dest_rl_rank, source_rl_rank;

    // Find destination rank (RIGHT to current rank)
    if (own_rank != iproc*jproc -1 ) { // Rank < MaxRank  is not existent
        dest_rl_rank = own_rank + 1;
    }

    // Find source rank (LEFT to current rank)
    if (own_rank != 0 ) { // Rank = -1  does not existent
        source_rl_rank = own_rank - 1;
    }

    // Exchange data
    if ((own_rank+1)%iproc == 0){ // receiver
        MPI_Sendrecv(&send_rl_y, domain.size_y + 2 , MPI_DOUBLE, MPI_PROC_NULL, 0,
                     &receive_rl_y, domain.size_y + 2 , MPI_DOUBLE, source_rl_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // TODO: MAYBE REMOVE LATER
        for(int i = 0; i < domain.size_y + 2; i++){
            A(0, i) = receive_rl_y[i];
        }

    } else if (own_rank%iproc == 0){ // sender
        MPI_Sendrecv(&send_rl_y, domain.size_y + 2 , MPI_DOUBLE, dest_rl_rank, 0,
                     &receive_rl_y, domain.size_y + 2 , MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else { //MIDDLE: SEND and RECEIVE
        MPI_Sendrecv(&send_rl_y, domain.size_y + 2 , MPI_DOUBLE, dest_rl_rank, 0,
                     &receive_rl_y, domain.size_y + 2 , MPI_DOUBLE, source_rl_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // TODO: MAYBE REMOVE LATER
        for(int i = 0; i < domain.size_y + 2; i++){
            A(0, i) = receive_rl_y[i];
        }     
    }

    /*
    // Copy data from receive into Matrix A TOP ghost layer
    if (own_rank%iproc != 0){//prevent senders from receiving
        for(int i = 0; i < domain.size_y + 2; i++){
            A(0, i) = receive_rl_y[i];
        }
    }
    */
}




}

