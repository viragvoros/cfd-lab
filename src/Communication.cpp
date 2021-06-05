#include "Communication.hpp"

#include <mpi.h>

#include <iostream>

void Communication::communicate(Matrix<double> &A, Domain &domain, int &iproc, int &jproc) {
    int own_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &own_rank);

//    //----------------------------------------------------------------
//    // SEND TO TOP, RECEIVE INTO BOTTOM
//    //----------------------------------------------------------------
//    double send_tb_y[domain.size_y + 2];
//    double receive_tb_y[domain.size_y + 2];
//    // Copy data to send the TOP non-ghost layer of Matrix A
//    for(int i = 0; i < domain.size_x + 2; i++){
//        send_tb_y[i] = A(i, domain.size_y);
//    }
//    // Find destination rank (TOP to current rank)
//    int dest_rank = own_rank + iproc;
//    // Find source rank (BOTTOM to current rank)
//    int source_rank = own_rank - iproc;
//    // Exchange data
//    if ((own_rank + iproc) >= (iproc * jproc)){
//        MPI_Sendrecv(&send_tb_y, domain.size_y + 2, MPI_DOUBLE, MPI_PROC_NULL, 0,
//            &receive_tb_y, domain.size_y + 2, MPI_DOUBLE, source_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    } else if (own_rank < iproc){
//        MPI_Sendrecv(&send_tb_y, domain.size_y + 2, MPI_DOUBLE, dest_rank, 0,
//            &receive_tb_y, domain.size_y + 2, MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    } else {
//        MPI_Sendrecv(&send_tb_y, domain.size_y + 2, MPI_DOUBLE, dest_rank, 0,
//            &receive_tb_y, domain.size_y + 2, MPI_DOUBLE, source_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    }
//    // Copy data from receive into Matrix A BOTTOM ghost layer
//    for(int i = 0; i < domain.size_x + 2; i++){
//        A(i, 0) = receive_tb_y[i];
//    }
//
//    //----------------------------------------------------------------
//    // SEND TO BOTTOM, RECEIVE FROM TOP
//    //----------------------------------------------------------------
//    double send_bt_y[domain.size_y + 2];
//    double receive_bt_y[domain.size_y + 2];
//    // Copy data to send the BOTTOM non-ghost layer of Matrix A
//    for(int i = 0; i < domain.size_x + 2; i++){
//        send_bt_y[i] = A(i, 1);
//    }
//    // Find destination rank (TOP to current rank)
//    dest_rank = own_rank - iproc;
//    // Find source rank (BOTTOM to current rank)
//    source_rank = own_rank + iproc;
//    // Exchange data
//    if ((own_rank + iproc) >= (iproc * jproc)){ // TOP ROW: ONLY SEND
//        MPI_Sendrecv(&send_bt_y, domain.size_y + 2, MPI_DOUBLE, dest_rank, 0,
//                     &receive_bt_y, domain.size_y + 2, MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    } else if (own_rank < iproc){ //BOTTOM ROW: ONLY RECEIVE
//        MPI_Sendrecv(&send_bt_y, domain.size_y + 2, MPI_DOUBLE, MPI_PROC_NULL, 0,
//                     &receive_bt_y, domain.size_y + 2, MPI_DOUBLE, source_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    } else { //MIDDLE: SEND and RECEIVE
//        MPI_Sendrecv(&send_bt_y, domain.size_y + 2, MPI_DOUBLE, dest_rank, 0,
//                     &receive_bt_y, domain.size_y + 2, MPI_DOUBLE, source_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    }
//    // Copy data from receive into Matrix A TOP ghost layer
//    for(int i = 0; i < domain.size_x + 2; i++){
//        A(i, domain.size_y + 1) = receive_bt_y[i];
//    }

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
    if (own_rank != 0) { // Rank = -1 is not existent
        dest_lr_rank = own_rank - 1;
    }

    // Find sender rank (RIGHT to current rank)
    if (own_rank != iproc*jproc -1 ) { // Rank < MaxRank  is not existent
      source_lr_rank = own_rank + 1;
    }


    // Exchange data
    if (own_rank == 1){ // LEFT COL: ONLY SEND
        MPI_Sendrecv(&send_lr_y, domain.size_y + 2, MPI_DOUBLE, dest_lr_rank, 0,
                     &receive_lr_y, domain.size_y + 2, MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if (own_rank == 0){ //RIGHT COL: ONLY RECEIVE
        MPI_Sendrecv(&send_lr_y, domain.size_y + 2, MPI_DOUBLE, MPI_PROC_NULL, 0,
                     &receive_lr_y, domain.size_y + 2 , MPI_DOUBLE, source_lr_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else { //MIDDLE: SEND and RECEIVE
        MPI_Sendrecv(&send_lr_y, domain.size_y + 2 , MPI_DOUBLE, dest_lr_rank, 0,
                     &receive_lr_y, domain.size_y + 2 , MPI_DOUBLE, source_lr_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


//    std::cout << A.jmax() << std::endl;

    // Copy data from receive into Matrix A Right ghost layer
    if (own_rank == 0){
        for(int i = 0; i < domain.size_y + 2 ; i++){
            A(domain.size_x + 1, i ) = receive_lr_y[i];
            //std::cout << i << std::endl;
//        if (own_rank == 1) {
//            std::cout << i;
//        }
        }
    }

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
    // SEND TO RIGHT, RECEIVE FROM LEFT
    //----------------------------------------------------------------
    // Copy data to send the RIGHT non-ghost layer of Matrix A
    double send_rl_y[domain.size_y + 2];
    double receive_rl_y[domain.size_y + 2];


        for (int i = 0; i < domain.size_y + 2; i++) {
            send_rl_y[i] = A(domain.size_x, i);
        }

    //std::cout << "domain.size_y + 2 + 2 :" << domain.size_y + 2 + 2  << std::endl;


//    // Find destination rank (LEFT to current rank)
//        int dest_rank_rl = own_rank - 1;
//
//    // Find source rank (RIGHT to current rank)
//        int source_rank_rl = own_rank + 1;

    // Exchange data
    if (own_rank == 1){ // LEFT COL: ONLY RECEIVE
        MPI_Sendrecv(&send_rl_y, domain.size_y + 2 , MPI_DOUBLE, MPI_PROC_NULL, 0,
                     &receive_rl_y, domain.size_y + 2 , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //std::cout << "receiving" << std::endl;
    } else if (own_rank == 0){ //RIGHT COL: ONLY SEND
        MPI_Sendrecv(&send_rl_y, domain.size_y + 2 , MPI_DOUBLE, 1, 0,
                     &receive_rl_y, domain.size_y + 2 , MPI_DOUBLE, MPI_PROC_NULL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       // std::cout << "sending" << std::endl;;
    } else { //MIDDLE: SEND and RECEIVE
        MPI_Sendrecv(&send_rl_y, domain.size_y + 2 , MPI_DOUBLE, 1, 0,
                     &receive_rl_y, domain.size_y + 2 , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::cout << "both" << std::endl;;
    }
    // Copy data from receive into Matrix A TOP ghost layer
    if (own_rank == 1){
        for(int i = 0; i < domain.size_y + 2; i++){
            A(0, i) = receive_rl_y[i];
        }
    }





}

