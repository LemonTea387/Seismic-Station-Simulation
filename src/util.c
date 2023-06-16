#include "../headers/util.h"

/**
* @brief function for generating random floating point
*/
float get_random_float(float min, float max) 
{
    return (((float)rand()/(float)(RAND_MAX)) * max) + min;
}

/**
* @brief function for generating date and time
*/
struct tm *generate_datetime() 
{
    time_t current = time(NULL);

    if (current == -1) {
        printf("Error! Cannot get the time");
        fflush(stdout);
        return NULL;
    }

    return localtime(&current);
}

/**
* @brief function for generating random floating between low and high
*/
double rand_double(double low, double high) 
{
    return (rand_() * (high - low)) + low;
}

/**
* @brief function for generating random float between offset to offset + 1 or 0 to offset
*/
double rand_single(double offset, int flag) 
{
    return flag == 0 ? rand_() + offset : rand_() * offset;
}

/**
* @brief function for generating random floating point between 0 to 1
*/
double rand_() 
{
    return (double) rand() / (double) RAND_MAX;
}

/**
* @brief Haversine Law for estimating the location between 2 points on earth. The return value is in km
*/
double distance(double latitude1, double longitude1, double latitude2, double longitude2) 
{
    double phi_1 = radians(latitude1), phi_2 = radians(latitude2);
    double delta_phi = radians(latitude2 - latitude1);
    double delta_lambda = radians(longitude2 - longitude1);
    double a = pow(sin(delta_phi / 2.0), 2) + (cos(phi_1) * cos(phi_2) * pow(sin(delta_lambda / 2.0), 2));
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    return 6371 * c; // unit = km
}

/**
* @brief convert angle degree to radians
*/
double radians(double angle) 
{
    return angle / 180 * M_PI;
}