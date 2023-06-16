#include <stdlib.h>
#include <mpi.h>

#include "../headers/network.h"
#include "../headers/seismic_sensor.h"
#include "../headers/base.h"
#include "../headers/struct_data.h"

// Global Variable
MPI_Datatype datetime_t;
MPI_Datatype data_t;
MPI_Datatype MPI_response_t;
MPI_Datatype MPI_report_t;
MPI_Datatype MPI_timespec_t;

int main(int argc, char * argv[]) 
{
    // Passing args 
    int m, n, delay;
    if (argc == 4) 
    {
        m = atoi(argv[1]);
        n = atoi(argv[2]);
        delay = atoi(argv[3]);
    }
    else 
    {
        printf("Usage : %s [m] [n] [delay]", argv[0]);
        return -1;
    }

    int rank, size, thread_level, rc = MPI_SUCCESS;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_level);

    if (thread_level != MPI_THREAD_MULTIPLE) {
            printf("Process Ended. Requirement is not satisfied.");
            fflush(stdout);
            rc = MPI_Finalize();
            return rc;
    }

    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Define the datatypes
    rc = commit_mpi_read_t(&datetime_t, &data_t);
    rc = commit_mpi_response_t(&MPI_response_t, &data_t);
    rc = commit_mpi_timespec(&MPI_timespec_t);
    rc = commit_mpi_report_t(&MPI_report_t, &MPI_response_t, &MPI_timespec_t);

    // Do seed the randomness
    srand(rank);

    /*
     *  Setting up Seismic Sensor Network in separate communicator
     *  We reserve rank 0 to be the base station. Others will be split to form
     *  Seismic sensor network.
     */
    MPI_Comm seismic_comm;
    rc = MPI_Comm_split( MPI_COMM_WORLD, rank == BASE_STATION, 0, &seismic_comm);

    if (rank != BASE_STATION) // Case where it's not the base station
    {
        rc = seismic_network(MPI_COMM_WORLD, seismic_comm, m, n, delay, &MPI_response_t, &MPI_report_t);
    }
    else if (rank == BASE_STATION ) // Base station code
    {
        rc = main_base(size, delay, m, n, &MPI_report_t);
    }

    // Free the datatype
    rc = MPI_Type_free(&datetime_t);
    rc = MPI_Type_free(&data_t);
    rc = MPI_Type_free(&MPI_report_t);
    rc = MPI_Type_free(&MPI_response_t);

    rc = MPI_Finalize();

    if (rc == MPI_SUCCESS) {
        printf("Programs Ended Successfully\n");
    } else {
        printf("Programs Crashes Somewhere!!!\n");
    }

    return rc;
}