#ifndef NETWORK_H
#define NETWORK_H

#define BASE_STATION 0

/**
 * The signal to quit, upon receiving such signals, the system has to prepare
 * to shut down, and gracefully exit.
 */
#define EXIT_SIGNAL 666

#define NODE_REQUEST_SIGNAL 1
#define NODE_RESPONSE_SIGNAL 2
#define NODE_REPORT_SIGNAL 3
#define BASE_PING_SIGNAL 4
#define NODE_ALIVE_SIGNAL 5

#endif // NETWORK_H