#ifndef BASE_H
#define BASE_H

#include <mpi.h>
#include <omp.h>

#include "util.h"
#include "network.h"
#include "struct_data.h"

#define BASE_THREAD_NUM 4
#define END_PROGRAM 'e'
#define GLOBAL_ARRAY_SIZE 10

#define COLUMN_WIDTH 10

#define DEPTH_MAX 700.0
#define MAGNITUDE_MAX 10.0

#define MAGNITUDE_THRESHOLD 2.5
#define LATITUDE_DIFF_THRESHOLD 3.0
#define LONGITUDE_DIFF_THRESHOLD 3.0
#define MAGNITUDE_DIFF_THRESHOLD 3.0
#define DEPTH_DIFF_THRESHOLD 250.0
#define DISTANCE_THRESHOLD 300.0

#define SIMULATION_TIME 120
#define PING_DELAY 10
#define PING_TIME_THRESHOLD 2

int _sensors_handler(char *exit_flag, int delay, double *comm_time, int *compare_signal, report_t *recv_buffer, MPI_Datatype *MPI_report_t);
int _compare_and_write(char *exit_flag, int *compare_signal, report_t *recv_buffer, double *comm_time);
int main_base(int p, int delay, int m, int n, MPI_Datatype *MPI_report_t);
int _fault_detect_thread(char *exit_flag, int p);

#endif //BASE_H