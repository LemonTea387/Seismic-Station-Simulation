#include "../headers/base.h"
#include <ctype.h>

// Global array for base station and balloon
read_t *gReading;
int report_ind = -1;

/**
* @brief This function is called by thread handling the communications between the sensors with the base station
* @param exit_flag signal tells that the thread should be end
* @param delay period defined for checking the incoming request
* @param comm_time pointer to store the result of communication
* @param compare_signal pointer to store the signal whether comparison should be done
* @param recv_buffer pointer that used as buffer as receiving data from the sensors
* @param MPI_report_t pointer to defined datatype for communication purposes
* @return flag telling whether the behaviour is functioning correctly
*/
int _sensors_handler(char *exit_flag, int delay, double *comm_time, int *compare_signal, report_t *recv_buffer, MPI_Datatype *MPI_report_t)
{
    // Initialize the variables used in communication
    int msg_arrival = 0;
    MPI_Status status;
    
    // Loop until the exit signal is received
    while (tolower(*exit_flag) != END_PROGRAM) {
        // Check whether there is any report sent
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &msg_arrival, &status);

        // If report is sent and comparison has finished, get the new report
        if (msg_arrival == 1 && (*compare_signal) == 0) {
            MPI_Recv(recv_buffer, 1, *MPI_report_t, status.MPI_SOURCE, NODE_REPORT_SIGNAL, MPI_COMM_WORLD, &status);
            // printf("[Base] Received report from [Node %d]. %d(%d, %d) %d %d %d\n", 
            //     status.MPI_SOURCE, recv_buffer->origin.node_rank, recv_buffer->origin.coord_x, recv_buffer->origin.coord_y, recv_buffer->responses, recv_buffer->matches, recv_buffer->total_messages);
            // fflush(stdout);
            
            struct timespec curr_time;
            clock_gettime(CLOCK_REALTIME, &curr_time);
            
            #pragma omp critical
            *comm_time = (((curr_time.tv_sec - (recv_buffer -> start_report_time).tv_sec)*1e9) + (curr_time.tv_nsec - (recv_buffer->start_report_time).tv_nsec))*1e-9;

            // Update the signal
            #pragma omp critical
            *compare_signal = 1;
        }
        
        // To mimick the effect of periodically probing
        // sleep(delay);
    }

    return MPI_SUCCESS;
}

/**
* @brief This function is called by thread handling comparison and logging
* @param exit_flag signal tells the program to end
* @param compare_signal signal tells that the report is received and ready for logging
* @param recv_buffer buffer where the received report is stored
* @param comm_time communication time between sensor and base station
* @return flag telling whether the process is successful ended
*/
int _compare_and_write(char *exit_flag, int *compare_signal, report_t *recv_buffer, double *comm_time) 
{
    struct timespec start, end;
    double time_taken = 0;

    // Initialize variable used to produce a file
    FILE* pFile = fopen("./results/report.txt", "w");
    int report_counter = 0;

    // Only exit when the signal is received (input from the user)
    clock_gettime(CLOCK_MONOTONIC, &start);
    while (tolower(*exit_flag) != END_PROGRAM) {
        // Check whether report is received, if received, perform comparison and logging
        if ((*compare_signal) == 1) {
            // Initialize variables
            int conclusive_flag = 0;

            // Compare the differences against the threshold for conclusion
            while (report_ind < 0) {}

            printf("[Base] Logs the report\n");
            fflush(stdout);
            
            read_t reading;
            reading = gReading[report_ind];

            double lat_diff = fabs(recv_buffer -> origin.reading.latitude - reading.latitude);
            double long_diff = fabs(recv_buffer -> origin.reading.longitude - reading.longitude);
            double magnitude_diff = fabs(recv_buffer -> origin.reading.magnitude - reading.magnitude);
            double depth_diff = fabs(recv_buffer -> origin.reading.depth - reading.magnitude);
            double location = distance(reading.latitude, reading.longitude, recv_buffer -> origin.reading.latitude, recv_buffer -> origin.reading.longitude);
            
            if(lat_diff <= LATITUDE_DIFF_THRESHOLD && long_diff <= LONGITUDE_DIFF_THRESHOLD
                && magnitude_diff <= MAGNITUDE_DIFF_THRESHOLD && depth_diff <= DEPTH_DIFF_THRESHOLD
                && location <= DISTANCE_THRESHOLD)
            {
                conclusive_flag = 1;
            }

            // Logging the report to txt file
            char buffer[26];

            // Log the time
            strftime(buffer, 26, "%Y:%m:%d %H:%M:%S", generate_datetime());
            fprintf(pFile, "Report : %d\nLogged Time : %s\n", (++report_counter), buffer);
            
            // Log the reporting node
            strftime(buffer, 26, "%Y:%m:%d %H:%M:%S", &(recv_buffer -> origin.reading.datetime));
            fprintf(pFile, "Alert Report Time : %s\nAlert Type : %s\n\nReporting Node,Seismic Coord, Magnitude\n%d (%d, %d),(%lf, %lf),%lf\n\n", 
                buffer, conclusive_flag == 0 ? "inconclusive" : "conclusive", 
                recv_buffer -> origin.node_rank, recv_buffer -> origin.coord_x, recv_buffer -> origin.coord_y,
                recv_buffer -> origin.reading.latitude, recv_buffer -> origin.reading.longitude,
                recv_buffer -> origin.reading.magnitude);
            
            // Log the adjacent nodes
            fprintf(pFile, "Adjacent Nodes, Seismic Coord, Location_Diff, Magnitude, Mag_Diff\n");
            for(int i = 0; i < recv_buffer -> responses; i++) 
            {
                fprintf(pFile, "%d (%d, %d),(%lf, %lf),%lf,%lf,%lf\n", 
                    recv_buffer -> neighbours[i].node_rank, recv_buffer -> neighbours[i].coord_x, recv_buffer -> neighbours[i].coord_y,
                    recv_buffer -> neighbours[i].reading.latitude, recv_buffer -> neighbours[i].reading.longitude,
                    distance(recv_buffer -> origin.reading.latitude, recv_buffer -> origin.reading.longitude, recv_buffer -> neighbours[i].reading.latitude, recv_buffer -> neighbours[i].reading.longitude),
                    recv_buffer -> neighbours[i].reading.magnitude,
                    fabs(recv_buffer -> origin.reading.magnitude - recv_buffer -> neighbours[i].reading.magnitude));
            }
            
            // Log the balloon seismic node
            strftime(buffer, 26, "%Y:%m:%d %H:%M:%S", &reading.datetime);
            fprintf(pFile, "\nBalloon Seismic Report Time : %s\nBallon Seismic reporting Coord : (%lf, %lf)\nBalloon Seismic Location Different with reporting node : %lf\nBalloon Seismic Magnitude : %lf\nBalloon Seismic Magnitude Different with reporting node : %lf\n\n", 
                buffer, reading.latitude, reading.longitude, 
                distance(recv_buffer -> origin.reading.latitude, recv_buffer -> origin.reading.longitude, reading.latitude, reading.longitude),
                reading.magnitude, fabs(recv_buffer -> origin.reading.magnitude - reading.magnitude));

            // Log the performance metrics
            fprintf(pFile, "Communication Time (seconds) : %lf\nTotal Message between reporting nodes and adjacent nodes : %d\nNumber of adjacent matches with reporting node : %d\nLocation Threshold : %f\nMagnitude Threshold : %f\nMagnitude Difference Threshold : %f\n\n", 
            *comm_time, recv_buffer -> total_messages, recv_buffer -> matches, DISTANCE_THRESHOLD, MAGNITUDE_THRESHOLD, MAGNITUDE_DIFF_THRESHOLD);

            #pragma omp critical
            *compare_signal = 0;
        }

        // Check whether there is any input from the console to end the program
        struct timeval tv = { 0L, 0L };
        fd_set fds;
        FD_ZERO(&fds);
        FD_SET(0, &fds);

        if (select(1, &fds, NULL, NULL, &tv) > 0) {
            scanf("%c", exit_flag);
        }
        
        // for running on CAAS
        clock_gettime(CLOCK_MONOTONIC, &end);
        time_taken = (end.tv_sec - start.tv_sec) * 1e9;
        time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
        // printf("Time taken %lf\n", time_taken);
        
        if (time_taken > SIMULATION_TIME) {
            printf("Shutdown the program\n");
            *exit_flag = END_PROGRAM;
        }
    }

    fclose(pFile);
    return MPI_SUCCESS;
}

/**
* @brief Thread to ping other nodes to make sure they are functioning
* @param exit_flag flag tells the process to end
* @param p number of processes in the program
*/
int _fault_detect_thread(char *exit_flag, int p) 
{
    // Initialize variables
    int buffer, msg_arrival = 0;
    struct timespec start, end;
    char *flags = (char*)malloc(sizeof(char) * (p-1));
    int ping_flag = 1, response_counter = p-1;
    MPI_Status status;
    FILE* pFile = fopen("./results/ping.txt", "w");

    // Ping the nodes at fixed time interval
    while (tolower(*exit_flag) != END_PROGRAM) {
        if (ping_flag == 1) {
            for (int i = 1; i < p; i++) {
                MPI_Send(&buffer, 0, MPI_INT, i, BASE_PING_SIGNAL, MPI_COMM_WORLD);
            }
            clock_gettime(CLOCK_MONOTONIC, &start);
            ping_flag = 0;
        }

        // Check whether there is any report sent
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &msg_arrival, &status);

        // If report is sent and comparison has finished, get the new report
        if (msg_arrival == 1) {
            MPI_Recv(&buffer, 0, MPI_INT, status.MPI_SOURCE, NODE_ALIVE_SIGNAL, MPI_COMM_WORLD, &status);
            // printf("[Base] received ping response from [Node %d]\n", status.MPI_SOURCE);
            // fflush(stdout);
            response_counter--;
            flags[status.MPI_SOURCE-1] = 'a';

            if (response_counter == 0) {
                ping_flag = 1;
                response_counter = p-1;
                
                for (int i = 0; i < (p-1); i++) {
                    flags[i] = 's';
                }

                sleep(PING_DELAY);
            } else {
                clock_gettime(CLOCK_MONOTONIC, &end);
                double time_taken = (end.tv_sec - start.tv_sec) * 1e9;
                time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
                
                if (time_taken > PING_TIME_THRESHOLD) {
                    // Logging the report
                    for (int i = 0; i < (p-1); i++) {
                        if (flags[i] != 'a') {
                            fprintf(pFile, "[Node %d] does not response within time frame. [Node %d] has fault.\n", i+1, i+1);
                        }
                    }
                    printf("Some nodes are at fault. Shutdown the whole networks\n");
                    fflush(stdout);
                    *exit_flag = END_PROGRAM;
                }
            }
        }
    }

    fclose(pFile);
    free(flags);

    return 0;
}

/**
* @brief This is the main function that describes how base station should work
* @param p number of processes in this program
* @param delay determines period that the thread should check for any request
* @param m number of rows in the cartesian topology
* @param n number of columns in the cartesian topology
* @param MPI_report_t datatype defined for sending and receving the report
* @return flag telling whether the process is ended successfully
*/
int main_base(int p, int delay, int m, int n, MPI_Datatype *MPI_report_t) {
    // Initialize variables
    int rc = MPI_SUCCESS;
    gReading = (read_t*)malloc(GLOBAL_ARRAY_SIZE * sizeof(read_t));

    int compare_signal = 0;
    char exit_flag = 's';
    double comm_time = 0;
    report_t *recv_buffer = (report_t*)malloc(sizeof(report_t));

    // Set number of threads
    omp_set_num_threads(BASE_THREAD_NUM);

    // Run the threads according to their jobs
    #pragma omp parallel sections
    {
        #pragma omp section
        _compare_and_write(&exit_flag, &compare_signal, recv_buffer, &comm_time);

        #pragma omp section
        _sensors_handler(&exit_flag, delay, &comm_time, &compare_signal, recv_buffer, MPI_report_t);

        #pragma omp section
        {
            while (tolower(exit_flag) != END_PROGRAM) {
                #pragma omp critical
                {                
                    report_ind = (report_ind + 1) % GLOBAL_ARRAY_SIZE;
                    gReading[report_ind].datetime = *(generate_datetime()); 
                    gReading[report_ind].latitude = rand_double(-0.5, m+0.5); 
                    gReading[report_ind].longitude = rand_double(-0.5, n+0.5); 
                    gReading[report_ind].magnitude = rand_double(MAGNITUDE_THRESHOLD, MAGNITUDE_MAX);
                    gReading[report_ind].depth = rand_single(DEPTH_MAX, 1);
                }

                sleep(delay);
            }
        }
        
        // #pragma omp section
        // _fault_detect_thread(&exit_flag, p);
    }

    // Send signal to shutdown other nodes
    int termination_buffer;
    for (int i = 1; i < p; i++) {
        printf("[Base] Send termination signal to [Node %d]\n", i);
        fflush(stdout);
        rc = MPI_Send(&termination_buffer, 0, MPI_INT, i, EXIT_SIGNAL, MPI_COMM_WORLD);
    }

    printf("[Base] Exiting\n");
    fflush(stdout);

    free(recv_buffer);
    free(gReading);
    return rc;
}