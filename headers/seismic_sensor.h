#ifndef SEISMIC_H
#define SEISMIC_H

#include <stdio.h>
#include <unistd.h>
#include <omp.h>
#include <mpi.h>
#include <time.h>

#include "network.h"
#include "util.h"
#include "struct_data.h"

#define SHIFT_ROW 0 
#define SHIFT_COL 1
#define DISP 1

#define LATITUDE_MIN -90.0
#define LATITUDE_MAX 90.0
#define LONGITUDE_MIN -180.0
#define LONGITUDE_MAX 180.0
#define DEPTH_MAX 700.0
#define MAGNITUDE_MAX 10.0

#define MAGNITUDE_THRESHOLD 2.5
#define LATITUDE_DIFF_THRESHOLD 3.0
#define LONGITUDE_DIFF_THRESHOLD 3.0
#define MAGNITUDE_DIFF_THRESHOLD 3.0
#define DEPTH_DIFF_THRESHOLD 250.0
#define DISTANCE_THRESHOLD 300.0

typedef struct {
    MPI_Comm grid_comm;
    int cart_rank;
    int coord_x;
    int coord_y;
    // Neighbours
    int neighbours[4];  // Top, bottom, left, right
} node_t;

int seismic_network(MPI_Comm master_comm, MPI_Comm nodes_comm, int m, int n, int delay, MPI_Datatype * MPI_response_t, MPI_Datatype * MPI_report_t);

int _init_node(MPI_Comm nodes_comm, int rank, int m, int n, node_t * node_data);

int _network_handler(int * exit_flag, int * request_flag, int * compare_flag, int * report_flag, int * request_counter, int * compare_counter,int * total_message_count, MPI_Comm master_comm, node_t * node_data, read_t * reading_data, response_t * neighbour_result, report_t * report, MPI_Datatype * MPI_response_t, MPI_Datatype * MPI_report_t);

void _simulate_sensor(read_t * results, node_t * node_data);

#endif // SEISMIC_H