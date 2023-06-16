#include "../headers/seismic_sensor.h"

/**
* @brief This function is main function of the seismic sensors (slaves rank) forming the topology
* @param master_comm communicator that includes the master rank
* @param nodes_comm communicator that includes only the sensors
* @param m row of the cartesian topology
* @param n col of the cartesian topology
* @param delay period for checking the request / generating the readings
* @param MPI_response_t defined data type for sending the response / receiving response
* @param MPI_report_t defined data type for sending the report / receiving report
* @return flag tellings whether the functions runs successfully
*/
int seismic_network(MPI_Comm master_comm, MPI_Comm nodes_comm, int m, int n, int delay, MPI_Datatype * MPI_response_t, MPI_Datatype * MPI_report_t) 
{
    // Initialize variables
    int rank;
    int exit = 0;
    MPI_Comm_rank(nodes_comm, &rank);

    node_t node;
    if(_init_node(nodes_comm, rank, m, n, &node) != 0) 
    {
        printf("ERROR initializing seismic node network.\n");
        return -1;
    }

    // Reserve space for sensor readings
    read_t simulated_results;

    int request_flag = 0;
    int request_counter = 0;
    int compare_counter = 0;
    int comparison_flag = 0;
    int report_flag = 0;
    int total_messages_count = 0;
    read_t compare_result;
    response_t * neighbour_results = (response_t *) malloc(4 * sizeof(response_t));
    report_t report;
    
    // Run the threads as specified
    omp_set_num_threads(3);
    #pragma omp parallel sections
    {
        #pragma omp section
        _network_handler(&exit, &request_flag, &comparison_flag, &report_flag, &request_counter, &compare_counter ,&total_messages_count, 
            master_comm, &node, &simulated_results, neighbour_results, &report, MPI_response_t, MPI_report_t);

        #pragma omp section
        {
            while (!exit)
            {
                _simulate_sensor(&simulated_results, &node);
                // printf("Node : %d, Generated : %f %f %f %f\n", node.cart_rank, simulated_results.latitude, simulated_results.longitude, simulated_results.depth, simulated_results.magnitude);
                
                // If > Threshold, send request to neighbours
                if(simulated_results.magnitude >= MAGNITUDE_THRESHOLD && request_flag==0) 
                {
                    request_flag = 1;
                    total_messages_count = 0;
                    compare_result = simulated_results;
                    for(int i = 0; i<4; i++) 
                    {
                        if(node.neighbours[i] >= 0)
                        {
                            request_counter++;
                        }
                    }
                    compare_counter = request_counter;
                } 
                sleep(delay);
            }
        }

        #pragma omp section
        {
            // Handle comparison between the neighbours
            while (!exit) 
            {
                if (comparison_flag) {
                    int matches = 0;

                    for (int i = 0; i < compare_counter; i++) {
                        report.neighbours[i] = neighbour_results[i];

                        double lat_diff = fabs(neighbour_results[i].reading.latitude - compare_result.latitude);
                        double long_diff = fabs(neighbour_results[i].reading.longitude - compare_result.longitude);
                        double magnitude_diff = fabs(neighbour_results[i].reading.magnitude - compare_result.magnitude);
                        double depth_diff = fabs(neighbour_results[i].reading.depth - compare_result.depth);
                        double location = distance(compare_result.latitude, compare_result.longitude, neighbour_results[i].reading.latitude, neighbour_results[i].reading.longitude);
                        // printf("Node : %d, diff with neighbour %d, %lf %lf %lf %lf %lf\n",node.cart_rank, neighbour_results[i].node_rank ,lat_diff, long_diff, magnitude_diff, depth_diff, location);
                        if(lat_diff <= LATITUDE_DIFF_THRESHOLD && long_diff <= LONGITUDE_DIFF_THRESHOLD
                            && magnitude_diff <= MAGNITUDE_DIFF_THRESHOLD && depth_diff <= DEPTH_DIFF_THRESHOLD
                            && location <= DISTANCE_THRESHOLD)
                            {   
                                matches++;
                            }
                    }

                    // printf("Comparison done by %d, count : %d, matches : %d\n", node.cart_rank, compare_counter, matches);

                    if (matches > 1) {
                        clock_gettime(CLOCK_REALTIME, &report.start_report_time);
                        response_t origin = {node.cart_rank, node.coord_x, node.coord_y, compare_result};
                        report.origin = origin;
                        report.responses = compare_counter;
                        report.matches = matches;
                        report.total_messages = total_messages_count;
                        report_flag = 1;
                    }
                    request_flag = 0;
                    comparison_flag = 0;
                }
            }
        }
    }
    
    // Clear the memory before exit
    free(neighbour_results);
    return MPI_SUCCESS;
}

/**
* @brief This function is used to create the cartesian topology
* @param nodes_comm communicator involve only the sensors node
* @param rank rank of the sensor node
* @param m row of the cartesian topology
* @param n col of the cartesian topology
* @param node_data Struct storing the data about the node in the topology
* @return flags telling that the functions ended successfully
*/
int _init_node(MPI_Comm nodes_comm, int rank, int m, int n, node_t * node_data)
{
    // Set up Cartesian grid topology
    MPI_Comm new_comm;
    int size = m * n;
    int ndims=2;
    int dims[] = {m,n};
    int wrap_around[] = {0, 0}; // No wrap around
    int reorder = 1;
    int coord[ndims];
    int cart_rank;
    // Neighbours
    int nbr_i_lo, nbr_i_hi;
    int nbr_j_lo, nbr_j_hi;

    MPI_Dims_create(size, ndims, dims);
    int ierr;
    ierr = MPI_Cart_create(nodes_comm, ndims, dims, wrap_around, reorder, &new_comm);
    if(ierr != 0) {
        printf("ERROR[%d] creating Seismic Node Network\n",ierr);
        return -1;
    }
    
    MPI_Cart_coords(new_comm, rank, ndims, coord);
    /* use my cartesian coordinates to find my rank in cartesian group*/
    MPI_Cart_rank(new_comm, coord, &cart_rank);

    /* get my neighbors; axis is coordinate dimension of shift */
    /* axis=0 ==> shift along the rows: P[my_row-1]: P[me] : P[my_row+1] */
    /* axis=1 ==> shift along the columns P[my_col-1]: P[me] : P[my_col+1] */
    MPI_Cart_shift( new_comm, SHIFT_ROW, DISP, &nbr_i_lo, &nbr_i_hi );
    MPI_Cart_shift( new_comm, SHIFT_COL, DISP, &nbr_j_lo, &nbr_j_hi );

    // Populate Data
    node_data->grid_comm = new_comm;
    node_data->cart_rank = cart_rank;
    node_data->coord_x = coord[0];
    node_data->coord_y = coord[1];
    node_data->neighbours[0] = nbr_i_lo;
    node_data->neighbours[1] = nbr_i_hi;
    node_data->neighbours[2] = nbr_j_lo;
    node_data->neighbours[3] = nbr_j_hi;

    return 0;
}

/**
* @brief Function called by thread for handling communication between the nodes
* @param exit_flag flag to ends the thread
* @param request_flag flag to send request to neighbours, asking for data
* @param compare_flag flag that the buffer is filled and ready for comparison
* @param report_flag flag that tells the thread to send report to the base station
* @param request_counter counter telling that how many request yet received reply
* @param compare_counter counter telling that how many comparisons to make
* @param total_message_count count of total messages sent and received during communication
* @param master_comm communicator that includes the base station
* @param node_data struct containing node data in the topology
* @param reading_data struct containing readings obtained from the sensor
* @param neighbour_result arrays storing the readings related to the neighbours
* @param report report data to be sent
* @param MPI_response_t defined datatype for sending response / receiving response
* @param MPI_report_t defined datatype for sending report
*/
int _network_handler(int * exit_flag, int * request_flag, int * compare_flag, int * report_flag, int * request_counter, int * compare_counter,int * total_message_count, MPI_Comm master_comm, node_t * node_data, read_t * reading_data, response_t * neighbour_result, report_t * report, MPI_Datatype * MPI_response_t, MPI_Datatype * MPI_report_t)
{
    // Initialize the variables
    int master_msg_arrival = 0 , network_msg_arrival = 0;
    MPI_Status master_status, network_status, status;
    int buf;
    
    // Loop until the termination signal reached
    while(!(*exit_flag)) 
    {
        // Probe to check whether there is any request / message
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, master_comm, &master_msg_arrival, &master_status);
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, node_data->grid_comm, &network_msg_arrival, &network_status);

        if(network_msg_arrival) 
        {
            // Handler for Network messages
            switch (network_status.MPI_TAG)
            {
                case NODE_REQUEST_SIGNAL:
                    // Receive request from the source
                    //printf("[Node %d] SIG [REQ] Received from %d.\n", node_data ->cart_rank, network_status.MPI_SOURCE);
                    MPI_Recv(&buf, 0, MPI_INT, network_status.MPI_SOURCE, NODE_REQUEST_SIGNAL, node_data->grid_comm, &network_status);
                    
                    // Send the data to the source
                    response_t send_response = {
                        node_data -> cart_rank, 
                        node_data -> coord_x, 
                        node_data -> coord_y, 
                        {
                            (*reading_data).datetime,
                            (*reading_data).latitude,
                            (*reading_data).longitude,
                            (*reading_data).magnitude,
                            (*reading_data).depth
                            }};
                    MPI_Send(&send_response, 1, *MPI_response_t, network_status.MPI_SOURCE, NODE_RESPONSE_SIGNAL, node_data->grid_comm);
                    break;
                case NODE_RESPONSE_SIGNAL:
                    // Receive data from the source
                    if ((*request_counter) > 0) {
                        //printf("[Node %d] SIG [RESP] Received from %d.\n",node_data -> cart_rank, network_status.MPI_SOURCE);
                        // Receive data from the neighbours
                        MPI_Recv(&neighbour_result[(*compare_counter) - (*request_counter)], 1, *MPI_response_t, network_status.MPI_SOURCE, NODE_RESPONSE_SIGNAL, node_data -> grid_comm, &status);
                        (*request_counter)-=1;
                        (*total_message_count) += 1;
                        if ((*request_counter) == 0)
                        {
                            *compare_flag = 1;
                        }
                    }
                default:
                    break;
            }
        }
        
        // Handler for base station messages
        if(master_msg_arrival) 
        {
            switch (master_status.MPI_TAG)
            {
                case EXIT_SIGNAL:
                    *exit_flag = 1;
                    MPI_Recv(&buf, 0, MPI_INT, master_status.MPI_SOURCE, EXIT_SIGNAL, master_comm, &master_status);
                    printf("[Node %d] Exiting.\n", node_data->cart_rank);
                    break;
                case BASE_PING_SIGNAL:
                    // printf("[Node %d] Received Ping from BASE\n", node_data->cart_rank);
                    MPI_Recv(&buf, 0, MPI_INT, master_status.MPI_SOURCE, BASE_PING_SIGNAL, master_comm, &master_status);
                    MPI_Send(&buf, 0, MPI_INT, master_status.MPI_SOURCE, NODE_ALIVE_SIGNAL, master_comm);
                    break;
                default:
                    break;
            }
        }

        // Handler for sending request
        #pragma omp critical
        {
            if (*request_flag == 1 && !(*compare_flag)) 
            {
                // printf("[Node %d] Requesting data from neighbours\n", node_data ->cart_rank);
                for(int i = 0; i < 4; i++)
                {
                    if(node_data->neighbours[i] >= 0)
                    {
                        MPI_Send(&buf, 0, MPI_INT, node_data->neighbours[i], NODE_REQUEST_SIGNAL, node_data->grid_comm);
                        *total_message_count += 1;
                    }
                }
                *request_flag = 2;
            }
        }
        
        #pragma omp critical
        {
            // Handler for sending report
            if (*report_flag) 
            {
                printf("[Node %d] sending Report to BASE STATION\n", node_data->cart_rank);
                MPI_Send(report, 1, *MPI_report_t, BASE_STATION, NODE_REPORT_SIGNAL, master_comm);
                *report_flag = 0;
            }
        }
    }
    return 0;
}

/**
* @brief This would generate the random values
* @param results buffer for storing the readings
* @param node_data struct containing data of node in topology
*/
void _simulate_sensor(read_t * results, node_t * node_data)
{
    #pragma omp critical
    {
        // Result Generation
        results->datetime = *(generate_datetime());
        results->latitude = node_data->coord_x + rand_double(-0.5, 0.5);
        results->longitude = node_data->coord_y + rand_double(-0.5, 0.5);
        results->depth = rand_double(0.0, DEPTH_MAX);
        results->magnitude = rand_double(0.0, MAGNITUDE_MAX);
    }
}
