#include "lennard_jones.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SIGMA 1
#define EPSILON 1
#define RADIUS (2.5 * SIGMA)

double power(double value, int exponent) {
    if (value == 1 || exponent == 1){
        // don't need to calculate if base or exponent is 1
        return value;
    }
    double final_value = value;
    for (int i = 1; i < exponent; i++){
        final_value *= value;
    }
    return final_value;
}

void lennard_jones(double* particle1, double* particle2, double* distances, double* forces, double* ljforce){
    double abs_x_dist = distances[0];
    double abs_y_dist = distances[1];

    double distance = sqrtf((abs_x_dist*abs_x_dist) + (abs_y_dist*abs_y_dist));
    if (distance < 0.75) { // cap how close particles can get
        distance = 0.75;
    }
    //double lj_force = (double) (48 * EPSILON * power(SIGMA, 12) / power(distance, 13)) - (24 * EPSILON * power(SIGMA, 6) / power(distance, 7));
    double lj_force = (double) (48 / power(distance, 13)) - (24 / power(distance, 7));
    *ljforce += lj_force;

    // ensure we get negative or positive interactions where required, as distances passed in are absolute
    double x_dist = particle2[0] - particle1[0];
    double x_force = lj_force * x_dist;

    double y_dist = particle2[1] - particle1[1];
    double y_force = lj_force * y_dist;

    forces[0] = x_force;
    forces[1] = y_force;
}