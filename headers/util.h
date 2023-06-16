#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <sys/select.h>

float get_random_float(float min, float max);

struct tm *generate_datetime();
double rand_double(double low, double high);
double rand_single(double offset, int flag);
double rand_();
double distance(double latitude1, double longitude1, double latitude2, double longitude2);
double radians(double angle);

#endif // UTIL_H