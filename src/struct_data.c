#include "../headers/struct_data.h"

/**
* @brief functions used to create mpi datatype for readings
*/
int commit_mpi_read_t(MPI_Datatype *datetime_t, MPI_Datatype *data_t) 
{
    int rc = MPI_SUCCESS, length_time[] = {1, 1, 1, 1, 1, 1, 1, 1, 1}, length_read[] = {1, 1, 1, 1, 1};
    struct tm dummy_tm = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    MPI_Aint displacements_time[9];
    MPI_Aint base_address;

    MPI_Get_address(&dummy_tm, &base_address);
    MPI_Get_address(&dummy_tm.tm_sec, &displacements_time[0]);
    MPI_Get_address(&dummy_tm.tm_min, &displacements_time[1]);
    MPI_Get_address(&dummy_tm.tm_hour, &displacements_time[2]);
    MPI_Get_address(&dummy_tm.tm_mday, &displacements_time[3]);
    MPI_Get_address(&dummy_tm.tm_mon, &displacements_time[4]);
    MPI_Get_address(&dummy_tm.tm_year, &displacements_time[5]);
    MPI_Get_address(&dummy_tm.tm_wday, &displacements_time[6]);
    MPI_Get_address(&dummy_tm.tm_yday, &displacements_time[7]);
    MPI_Get_address(&dummy_tm.tm_isdst, &displacements_time[8]);

    for (int i = 0; i < 9; i++) {
        displacements_time[i] = MPI_Aint_diff(displacements_time[i], base_address);
    }

    MPI_Datatype types_time[] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    rc = MPI_Type_create_struct(9, length_time, displacements_time, types_time, datetime_t);
    rc = MPI_Type_commit(datetime_t);

    read_t dummy_read = {dummy_tm, 0.0, 0.0, 0.0, 0.0};
    MPI_Aint displacements_data[5];

    // Datetime struct may not be the same size as time_t struct... Need to check on it
    MPI_Get_address(&dummy_read, &base_address);
    MPI_Get_address(&dummy_read.datetime, &displacements_data[0]);
    MPI_Get_address(&dummy_read.latitude, &displacements_data[1]);
    MPI_Get_address(&dummy_read.longitude, &displacements_data[2]);
    MPI_Get_address(&dummy_read.magnitude, &displacements_data[3]);
    MPI_Get_address(&dummy_read.depth, &displacements_data[4]);

    for (int i = 0; i < 5; i++) {
        displacements_data[i] = MPI_Aint_diff(displacements_data[i], base_address);
    }

    MPI_Datatype types_data[] = {*datetime_t, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    rc = MPI_Type_create_struct(5, length_read, displacements_data, types_data, data_t);
    rc = MPI_Type_commit(data_t);

    return rc;
}

/**
* @brief function used to create mpi datatype for response
*/
int commit_mpi_response_t(MPI_Datatype *MPI_response_t, MPI_Datatype *data_t) 
{
    int rc = MPI_SUCCESS, length_response[] = {1, 1, 1, 1};

    response_t dummy_response = {
        0, 0, 0, 
        {
            {1, 1, 1, 1, 1, 1, 1, 1, 1}, 
            0.0, 0.0, 0.0, 0.0
        }
    };
    MPI_Aint displacements_response[4];
    MPI_Aint base_address;

    // Datetime struct may not be the same size as time_t struct... Need to check on it
    MPI_Get_address(&dummy_response, &base_address);
    MPI_Get_address(&dummy_response.node_rank, &displacements_response[0]);
    MPI_Get_address(&dummy_response.coord_x, &displacements_response[1]);
    MPI_Get_address(&dummy_response.coord_y, &displacements_response[2]);
    MPI_Get_address(&dummy_response.reading, &displacements_response[3]);

    for (int i = 0; i < 4; i++) {
        displacements_response[i] = MPI_Aint_diff(displacements_response[i], base_address);
    }

    MPI_Datatype types_data[] = {MPI_INT, MPI_INT, MPI_INT, *data_t};
    rc = MPI_Type_create_struct(4, length_response, displacements_response, types_data, MPI_response_t);
    rc = MPI_Type_commit(MPI_response_t);

    return rc;
}

/**
* function used to create mpi datatype for timespec
*/
int commit_mpi_timespec(MPI_Datatype *MPI_timespec_t) 
{
    int rc = MPI_SUCCESS, length_timespec[] = {1,1};
    MPI_Aint base_addr;

    struct timespec dummy_time = {1, 1};
    MPI_Aint displacements_timespec[2];
    MPI_Get_address(&dummy_time, &base_addr);
    MPI_Get_address(&dummy_time.tv_sec, &displacements_timespec[0]);
    MPI_Get_address(&dummy_time.tv_nsec, &displacements_timespec[1]);

    for (int i = 0; i < 2; i++) {
        displacements_timespec[i] = MPI_Aint_diff(displacements_timespec[i], base_addr);
    }

    MPI_Datatype types_timespec[] = {MPI_LONG, MPI_LONG};
    rc = MPI_Type_create_struct(2, length_timespec, displacements_timespec, types_timespec, MPI_timespec_t);
    rc = MPI_Type_commit(MPI_timespec_t);

    return rc;
}

/**
* @brief function used to create mpi datatype for report
*/
int commit_mpi_report_t(MPI_Datatype *MPI_report_t, MPI_Datatype *MPI_response_t, MPI_Datatype *MPI_timespec_t) 
{
    int rc = MPI_SUCCESS, length_report[] = {1,1,4,1,1,1};
    MPI_Aint base_address;

    response_t dummy_response = {
        0, 0, 0, 
        {
            {1, 1, 1, 1, 1, 1, 1, 1, 1}, 
            0.0, 0.0, 0.0, 0.0
        }
    };
    
    report_t dummy_report = {
        {1, 1},
        dummy_response,
        {dummy_response, dummy_response, dummy_response, dummy_response},
        0, 0, 0
    };
    
    MPI_Aint displacements_report[6];

    // Datetime struct may not be the same size as time_t struct... Need to check on it
    MPI_Get_address(&dummy_report, &base_address);
    MPI_Get_address(&dummy_report.start_report_time, &displacements_report[0]);
    MPI_Get_address(&dummy_report.origin, &displacements_report[1]);
    MPI_Get_address(&dummy_report.neighbours, &displacements_report[2]);
    MPI_Get_address(&dummy_report.responses, &displacements_report[3]);
    MPI_Get_address(&dummy_report.matches, &displacements_report[4]);
    MPI_Get_address(&dummy_report.total_messages, &displacements_report[5]);

    for (int i = 0; i < 6; i++) {
        displacements_report[i] = MPI_Aint_diff(displacements_report[i], base_address);
    }

    MPI_Datatype types_report[] = {*MPI_timespec_t, *MPI_response_t, *MPI_response_t, MPI_INT, MPI_INT, MPI_INT};
    rc = MPI_Type_create_struct(6, length_report, displacements_report, types_report, MPI_report_t);
    rc = MPI_Type_commit(MPI_report_t);

    return rc;
}