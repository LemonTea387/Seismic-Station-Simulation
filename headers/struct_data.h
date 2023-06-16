#ifndef STRUCT_DATA_H
#define STRUCT_DATA_H

#include <mpi.h>
#include <time.h>

typedef struct {
    struct tm datetime;
    double latitude;
    double longitude;
    double magnitude;
    double depth;
} read_t;

typedef struct {
    int node_rank;
    int coord_x;
    int coord_y;
    read_t reading;
} response_t;

typedef struct {
    struct timespec start_report_time;
    response_t origin;
    response_t neighbours[4];
    int responses;
    int matches;
    int total_messages;
} report_t;

int commit_mpi_read_t(MPI_Datatype *datetime_t, MPI_Datatype *data_t);
int commit_mpi_response_t(MPI_Datatype *MPI_response_t, MPI_Datatype * data_t);
int commit_mpi_report_t(MPI_Datatype *MPI_report_t, MPI_Datatype * MPI_response_t, MPI_Datatype * MPI_timespec_t);
int commit_mpi_timespec(MPI_Datatype *MPI_timespec_t);

#endif //STRUCT_DATA_H