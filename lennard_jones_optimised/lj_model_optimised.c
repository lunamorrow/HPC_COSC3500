#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "lennard_jones.h"

#define NUM_BOXES 9
#define RADIUS 2.5

int seed = 5;

int distance_checks(double* particle1, double* particle2, double* distances, double box_dimensions) {
    double x_dist;
    double y_dist;

    // x distance check
    if (fabs(particle2[0] - particle1[0]) > 5){
        // closest interaction will be across the PBC (if possible)
        if (particle2[0] - particle1[0] > 0){
            // particle2 on left and particle 1 on right of box
            if (((box_dimensions - particle2[0]) + particle1[0]) <= RADIUS){
                // particles can interact over PBC
                x_dist = (box_dimensions - particle2[0]) + particle1[0];
            } else {
                // particles cannot interact over PBC
                return 1;
            }
        } else {
            // particle1 on left and particle 2 on right of box
            if (((box_dimensions - particle1[0]) + particle2[0]) <= RADIUS){
                // particles can interact over PBC
                x_dist = (box_dimensions - particle1[0]) + particle2[0];
            } else {
                // particles cannot interact over PBC
                return 1;
            }
        }
    } else if (fabs(particle2[0] - particle1[0]) <= RADIUS) {
        // normal interaction inside box boundaries
        x_dist = (particle2[0] - particle1[0]);
    } else {
        // no interaction, particles too far away to interact within box or over PBC
        return 1;
    }

    if (fabs(particle2[1] - particle1[1]) > 5){
        // closest interaction will be across the PBC (if possible)
        if (particle2[1] - particle1[1] > 0){
            // particle2 above particle 1 in box
            if (((box_dimensions - particle2[1]) + particle1[1]) <= RADIUS){
                // particles can interact over PBC
                y_dist = (box_dimensions - particle2[1]) + particle1[1];
            } else {
                // particles cannot interact over PBC
                return 1;
            }
        } else {
            // particle1 above particle 2 in box
            if (((box_dimensions - particle1[1]) + particle2[1]) <= RADIUS){
                // particles can interact over PBC
                y_dist = (box_dimensions - particle1[1]) + particle2[1];
            } else {
                // particles cannot interact over PBC
                return 1;
            }
        }
    } else if (fabs(particle2[1] - particle1[1]) <= RADIUS) {
        // normal interaction inside box boundaries
        y_dist = (particle2[1] - particle1[1]);
    } else {
        // no interaction, particles too far away to interact within box or over PBC
        return 1;
    }

    distances[0] = x_dist;
    distances[1] = y_dist;

    return 0;
}

void make_particles(double* velocities, double* positions, double* accelerations, double* masses, int num_positions, double* output){
    double centre_of_mass_velocity_total_x = 0;
    double centre_of_mass_velocity_total_y = 0;
    double total_mass = 0;

    double x_position = 0.0;
    double y_position = 0.0;

    // assist with random double generation
    double velocity_rand = ((double) 4 / RAND_MAX);
    double acceleration_rand = ((double) 2 / RAND_MAX);

    for (int i = 0; i < (num_positions-8); i = i + 8) {
        // particle 1
        velocities[i] = (rand() * velocity_rand) + (-2); // x velocity 
        velocities[i+1] = (rand() * velocity_rand) + (-2); // y velocity
        positions[i] = x_position;
        positions[i+1] = y_position;
        accelerations[i] = (rand() * acceleration_rand) + 1; // x acceleration
        accelerations[i+1] = (rand() * acceleration_rand) + 1;// y acceleration
        masses[i] = 1;
        masses[i+1] = 1; // TO DO: refactor

        x_position += 2.5;
        if (x_position >= 20){
            x_position = 0;
            y_position += 2.5;
        }

        // particle 2
        velocities[i+2] = (rand() * velocity_rand) + (-2); // x velocity 
        velocities[i+3] = (rand() * velocity_rand) + (-2); // y velocity
        positions[i+2] = x_position;
        positions[i+3] = y_position;
        accelerations[i+2] = (rand() * acceleration_rand) + 1; // x acceleration
        accelerations[i+3] = (rand() * acceleration_rand) + 1; // y acceleration
        masses[i+2] = 1;
        masses[i+3] = 1; // TO DO: refactor

        x_position += 2.5;
        if (x_position >= 20){
            x_position = 0;
            y_position += 2.5;
        }

        // particle 3
        velocities[i+4] = (rand() * velocity_rand) + (-2); // x velocity
        velocities[i+5] = (rand() * velocity_rand) + (-2); // y velocity
        positions[i+4] = x_position;
        positions[i+5] = y_position;
        accelerations[i+4] = (rand() * acceleration_rand) + 1; // x acceleration
        accelerations[i+5] = (rand() * acceleration_rand) + 1; // y acceleration
        masses[i+4] = 1;
        masses[i+5] = 1; // TO DO: refactor

        x_position += 2.5;
        if (x_position >= 20){
            x_position = 0;
            y_position += 2.5;
        }

        // particle 4
        velocities[i+6] = (rand() * velocity_rand) + (-2); // x velocity
        velocities[i+7] = (rand() * velocity_rand) + (-2); // y velocity
        positions[i+6] = x_position;
        positions[i+7] = y_position;
        accelerations[i+6] = (rand() * acceleration_rand) + 1; // x acceleration
        accelerations[i+7] = (rand() * acceleration_rand) + 1; // y acceleration
        masses[i+6] = 1;
        masses[i+7] = 1; // TO DO: refactor

        x_position += 2.5;
        if (x_position >= 20){
            x_position = 0;
            y_position += 2.5;
        }

        // particle 5
        velocities[i+8] = (rand() * velocity_rand) + (-2); // x velocity
        velocities[i+9] = (rand() * velocity_rand) + (-2); // y velocity
        positions[i+8] = x_position;
        positions[i+9] = y_position;
        accelerations[i+8] = (rand() * acceleration_rand) + 1; // x acceleration
        accelerations[i+9] = (rand() * acceleration_rand) + 1; // y acceleration
        masses[i+8] = 1;
        masses[i+9] = 1; // TO DO: refactor

        x_position += 2.5;
        if (x_position >= 20){
            x_position = 0;
            y_position += 2.5;
        }

        // particle 6
        velocities[i+10] = (rand() * velocity_rand) + (-2); // x velocity
        velocities[i+11] = (rand() * velocity_rand) + (-2); // y velocity
        positions[i+10] = x_position;
        positions[i+11] = y_position;
        accelerations[i+10] = (rand() * acceleration_rand) + 1; // x acceleration
        accelerations[i+11] = (rand() * acceleration_rand) + 1; // y acceleration
        masses[i+10] = 1;
        masses[i+11] = 1; // TO DO: refactor

        x_position += 2.5;
        if (x_position >= 20){
            x_position = 0;
            y_position += 2.5;
        }

        // particle 7
        velocities[i+12] = (rand() * velocity_rand) + (-2); // x velocity
        velocities[i+13] = (rand() * velocity_rand) + (-2); // y velocity
        positions[i+12] = x_position;
        positions[i+13] = y_position;
        accelerations[i+12] = (rand() * acceleration_rand) + 1; // x acceleration
        accelerations[i+13] = (rand() * acceleration_rand) + 1; // y acceleration
        masses[i+12] = 1;
        masses[i+13] = 1; // TO DO: refactor

        x_position += 2.5;
        if (x_position >= 20){
            x_position = 0;
            y_position += 2.5;
        }

        // particle 8
        velocities[i+14] = (rand() * velocity_rand) + (-2); // x velocity
        velocities[i+15] = (rand() * velocity_rand) + (-2); // y velocity
        positions[i+14] = x_position;
        positions[i+15] = y_position;
        accelerations[i+14] = (rand() * acceleration_rand) + 1; // x acceleration
        accelerations[i+15] = (rand() * acceleration_rand) + 1; // y acceleration
        masses[i+14] = 1;
        masses[i+15] = 1; // TO DO: refactor

        x_position += 2.5;
        if (x_position >= 20){
            x_position = 0;
            y_position += 2.5;
        }

        centre_of_mass_velocity_total_y += ((masses[i+1] * velocities[i+1]) + 
                                            (masses[i+3] * velocities[i+3]) + 
                                            (masses[i+5] * velocities[i+5]) + 
                                            (masses[i+7] * velocities[i+7]) + 
                                            (masses[i+9] * velocities[i+9]) + 
                                            (masses[i+11] * velocities[i+11]) + 
                                            (masses[i+13] * velocities[i+13]) + 
                                            (masses[i+15] * velocities[i+15]));
        centre_of_mass_velocity_total_x += ((masses[i] * velocities[i]) + 
                                            (masses[i+2] * velocities[i+2]) + 
                                            (masses[i+4] * velocities[i+4]) + 
                                            (masses[i+6] * velocities[i+6]) + 
                                            (masses[i+8] * velocities[i+8]) + 
                                            (masses[i+10] * velocities[i+10]) + 
                                            (masses[i+12] * velocities[i+12]) + 
                                            (masses[i+14] * velocities[i+14]));
        total_mass += (masses[i] + masses[i+2] + masses[i+4] + masses[i+6] + 
                       masses[i+8] + masses[i+10] + masses[i+12] + masses[i+14]);

        output[0] = centre_of_mass_velocity_total_x;
        output[1] = centre_of_mass_velocity_total_y;
        output[2] = total_mass;
    }

    // for last particles that don't fit in the loop
    int rem = 8*(num_positions/8);
    for (int i = rem; i < num_positions; i = i + 2) {
        velocities[i] = (rand() * velocity_rand) + (-2); // x velocity
        velocities[i+1] = (rand() * velocity_rand) + (-2); // y velocity
        positions[i] = x_position;
        positions[i+1] = y_position;
        accelerations[i] = (rand() * acceleration_rand) + 1; // x acceleration
        accelerations[i+1] = (rand() * acceleration_rand) + 1; // y acceleration
        masses[i] = 1;
        masses[i+1] = 1; // TO DO: refactor

        x_position += 2.5;
        if (x_position >= 20){
            x_position = 0;
            y_position += 2.5;
        }
        centre_of_mass_velocity_total_x += (masses[i] * velocities[i]);
        centre_of_mass_velocity_total_y += (masses[i+1] * velocities[i+1]);
        total_mass += masses[i];
    }

    output[0] = centre_of_mass_velocity_total_x;
    output[1] = centre_of_mass_velocity_total_y;
    output[2] = total_mass;
}

double set_temperature(double* velocities, double* masses, double comvx, double comvy, int num_positions, int num_particles){
    // set temperature    
    int df = 2;
    int num_constraints = 2;
    double desired_temp = 0.75;
    double two_x_kinetic_energy = 0;
    int rem = 8*(num_positions/8);

    for (int i = 0; i < (num_positions-8); i = i + 8) {
        velocities[i] -= comvx;
        velocities[i+1] -= comvy;

        velocities[i+2] -= comvx;
        velocities[i+3] -= comvy;

        velocities[i+4] -= comvx;
        velocities[i+5] -= comvy;

        velocities[i+6] -= comvx;
        velocities[i+7] -= comvy;

        velocities[i+8] -= comvx;
        velocities[i+9] -= comvy;

        velocities[i+10] -= comvx;
        velocities[i+11] -= comvy;

        velocities[i+12] -= comvx;
        velocities[i+13] -= comvy;

        velocities[i+14] -= comvx;
        velocities[i+15] -= comvy;

        // sum of x and y for kinetic
        two_x_kinetic_energy += ((masses[i] * velocities[i] * velocities[i]) + 
                                 (masses[i+1] * velocities[i+1] * velocities[i+1]) +
                                 (masses[i+2] * velocities[i+2] * velocities[i+2]) +
                                 (masses[i+3] * velocities[i+3] * velocities[i+3]) +
                                 (masses[i+4] * velocities[i+4] * velocities[i+4]) +
                                 (masses[i+5] * velocities[i+5] * velocities[i+5]) +
                                 (masses[i+6] * velocities[i+6] * velocities[i+6]) +
                                 (masses[i+7] * velocities[i+7] * velocities[i+7]) +
                                 (masses[i+8] * velocities[i+8] * velocities[i+8]) +
                                 (masses[i+9] * velocities[i+9] * velocities[i+9]) +
                                 (masses[i+10] * velocities[i+10] * velocities[i+10]) +
                                 (masses[i+11] * velocities[i+11] * velocities[i+11]) +
                                 (masses[i+12] * velocities[i+12] * velocities[i+12]) +
                                 (masses[i+13] * velocities[i+13] * velocities[i+13]) +
                                 (masses[i+14] * velocities[i+14] * velocities[i+14]) +
                                 (masses[i+15] * velocities[i+15] * velocities[i+15]));
    }

    for (int i = rem; i < num_positions; i = i + 2) {
        velocities[i] -= comvx;
        velocities[i+1] -= comvy;

        velocities[i+2] -= comvx;
        velocities[i+3] -= comvy;

        velocities[i+4] -= comvx;
        velocities[i+5] -= comvy;

        velocities[i+6] -= comvx;
        velocities[i+7] -= comvy;

        // sum of x and y for kinetic
        two_x_kinetic_energy += ((masses[i] * velocities[i] * velocities[i]) + 
                                 (masses[i+1] * velocities[i+1] * velocities[i+1]));
    }
    
    double temperature = two_x_kinetic_energy/(df * num_particles - num_constraints);
    double frac_temp_difference = temperature/desired_temp;
    double modify_temp = 1/sqrtf(frac_temp_difference);

    // modify velocities to set temperature
    for (int i = 0; i < (num_positions-8); i = i + 8) {
        velocities[i] *= modify_temp;
        velocities[i+1] *= modify_temp;
        velocities[i+2] *= modify_temp;
        velocities[i+3] *= modify_temp;
        velocities[i+4] *= modify_temp;
        velocities[i+5] *= modify_temp;
        velocities[i+6] *= modify_temp;
        velocities[i+7] *= modify_temp;
        velocities[i+8] *= modify_temp;
        velocities[i+9] *= modify_temp;
        velocities[i+10] *= modify_temp;
        velocities[i+11] *= modify_temp;
        velocities[i+12] *= modify_temp;
        velocities[i+13] *= modify_temp;
        velocities[i+14] *= modify_temp;
        velocities[i+15] *= modify_temp;
    }
    for (int i = rem; i < num_positions; i = i + 2) {
        velocities[i] *= modify_temp;
        velocities[i+1] *= modify_temp;
    }

    return temperature;
}

int main(int argv, char** argc){
    // start time
    clock_t tic = clock();
    srand(seed);

    // TO DO: add error handling and enable passing in for num of steps
    const int num_particles = atoi(argc[1]);
    const float box_dimensions = atoi(argc[2]);

    // set values
    double dt = 0.000001;
    double simulation_length = dt*1000000;
    int num_positions = num_particles + num_particles;
    double num_frames = simulation_length/dt;

    // make particles
    double* velocities = (double*) malloc(sizeof(double) * num_positions);
    double* positions = (double*) malloc(sizeof(double) * num_positions);
    double* accelerations = (double*) malloc(sizeof(double) * num_positions);
    double* masses = (double*) malloc(sizeof(double) * num_positions);
    double output[3];

    make_particles(&velocities[0], &positions[0], &accelerations[0], &masses[0], num_positions, &output[0]);

    // set temperature
    double centre_of_mass_velocity_x = output[0]/output[2]; // centre_of_mass_velocity_total_x/total_mass
    double centre_of_mass_velocity_y = output[1]/output[2]; // centre_of_mass_velocity_total_y/total_mass

    double temperature = set_temperature(&velocities[0], &masses[0], centre_of_mass_velocity_x, centre_of_mass_velocity_y, num_positions, num_particles);

    // // create array to store trajectory and store initial values
    // double** trajectory = (double**) malloc(sizeof(double*) * num_frames); //MAKE 1D -> makes quick to send w CUDA
    // for (int ti = 0; ti < num_frames; ti++){
    //     trajectory[ti] = (double*) malloc(sizeof(double) * num_positions); //to store x and y positions
    // }

    double distances[2];
    double forces[2];
    double particle_one[2];
    double particle_two[2];
    double total_forces = 0;
    double* old_accelerations = (double*) malloc(sizeof(double) * num_positions);

    // run the simulation
    for (int t = 0; t < num_frames; t++){
        // // save positions to trajectory but skip some steps
        // if (t%1000 == 0){
        //     int t_index = (int) (t/1000);
        //     for (int i = 0; i < num_positions; i = i + 2){
        //         trajectory[t_index][i] = positions[i];
        //         trajectory[t_index][i+1] = positions[i+1];
        //     }
        // }

        // // check temp at end
        // if (t == 999999){
        //     double temp = set_temperature(&velocities[0], &masses[0], centre_of_mass_velocity_x, centre_of_mass_velocity_y, num_positions, num_particles);
        //     printf("temperature at t 999999 = %f\n", temp);
        //     fflush(stdout);
        // }
        
        // for each particle, update position according to velocity and store old accerations
        memset(old_accelerations, 0, (sizeof(double) * num_positions));
        for (int p = 0; p < (num_positions-8); p = p + 8){
            // particle 1
            positions[p] += (velocities[p] * dt) + (0.5 * accelerations[p] * (dt * dt));
            if (positions[p] > box_dimensions){
                positions[p] -= box_dimensions; // if out of top of box, reenter at bottom of box
            }
            if (positions[p] < 0){
                positions[p] += box_dimensions; // if out of bottom of box, reenter at top of box
            }
            positions[p+1] += (velocities[p+1] * dt) + (0.5 * accelerations[p+1] * (dt * dt));
            if (positions[p+1] > box_dimensions){
                positions[p+1] -= box_dimensions; // if out left side of box, reenter at right side of box
            }
            if (positions[p+1] < 0){
                positions[p+1] += box_dimensions; // if out right side of box, reenter at left side of box
            }
            old_accelerations[p] = accelerations[p];
            old_accelerations[p+1] = accelerations[p+1];

            // particle 2
            positions[p+2] += (velocities[p+2] * dt) + (0.5 * accelerations[p+2] * (dt * dt));
            if (positions[p+2] > box_dimensions){
                positions[p+2] -= box_dimensions; // if out of top of box, reenter at bottom of box
            }
            if (positions[p+2] < 0){
                positions[p+2] += box_dimensions; // if out of bottom of box, reenter at top of box
            }
            positions[p+3] += (velocities[p+3] * dt) + (0.5 * accelerations[p+3] * (dt * dt));
            if (positions[p+3] > box_dimensions){
                positions[p+3] -= box_dimensions; // if out left side of box, reenter at right side of box
            }
            if (positions[p+3] < 0){
                positions[p+3] += box_dimensions; // if out right side of box, reenter at left side of box
            }
            old_accelerations[p+2] = accelerations[p+2];
            old_accelerations[p+3] = accelerations[p+3];

            // particle 3
            positions[p+4] += (velocities[p+4] * dt) + (0.5 * accelerations[p+4] * (dt * dt));
            if (positions[p+4] > box_dimensions){
                positions[p+4] -= box_dimensions; // if out of top of box, reenter at bottom of box
            }
            if (positions[p+4] < 0){
                positions[p+4] += box_dimensions; // if out of bottom of box, reenter at top of box
            }
            positions[p+5] += (velocities[p+5] * dt) + (0.5 * accelerations[p+5] * (dt * dt));
            if (positions[p+5] > box_dimensions){
                positions[p+5] -= box_dimensions; // if out left side of box, reenter at right side of box
            }
            if (positions[p+5] < 0){
                positions[p+5] += box_dimensions; // if out right side of box, reenter at left side of box
            }
            old_accelerations[p+4] = accelerations[p+4];
            old_accelerations[p+5] = accelerations[p+5];

            // particle 4
            positions[p+6] += (velocities[p+6] * dt) + (0.5 * accelerations[p+6] * (dt * dt));
            if (positions[p+6] > box_dimensions){
                positions[p+6] -= box_dimensions; // if out of top of box, reenter at bottom of box
            }
            if (positions[p+6] < 0){
                positions[p+6] += box_dimensions; // if out of bottom of box, reenter at top of box
            }
            positions[p+7] += (velocities[p+7] * dt) + (0.5 * accelerations[p+7] * (dt * dt));
            if (positions[p+7] > box_dimensions){
                positions[p+7] -= box_dimensions; // if out left side of box, reenter at right side of box
            }
            if (positions[p+7] < 0){
                positions[p+7] += box_dimensions; // if out right side of box, reenter at left side of box
            }
            old_accelerations[p+6] = accelerations[p+6];
            old_accelerations[p+7] = accelerations[p+7];

            // particle 5
            positions[p+8] += (velocities[p+8] * dt) + (0.5 * accelerations[p+8] * (dt * dt));
            if (positions[p+8] > box_dimensions){
                positions[p+8] -= box_dimensions; // if out of top of box, reenter at bottom of box
            }
            if (positions[p+8] < 0){
                positions[p+8] += box_dimensions; // if out of bottom of box, reenter at top of box
            }
            positions[p+9] += (velocities[p+9] * dt) + (0.5 * accelerations[p+9] * (dt * dt));
            if (positions[p+9] > box_dimensions){
                positions[p+9] -= box_dimensions; // if out left side of box, reenter at right side of box
            }
            if (positions[p+9] < 0){
                positions[p+9] += box_dimensions; // if out right side of box, reenter at left side of box
            }
            old_accelerations[p+8] = accelerations[p+8];
            old_accelerations[p+9] = accelerations[p+9];

            // particle 6
            positions[p+10] += (velocities[p+10] * dt) + (0.5 * accelerations[p+10] * (dt * dt));
            if (positions[p+10] > box_dimensions){
                positions[p+10] -= box_dimensions; // if out of top of box, reenter at bottom of box
            }
            if (positions[p+10] < 0){
                positions[p+10] += box_dimensions; // if out of bottom of box, reenter at top of box
            }
            positions[p+11] += (velocities[p+11] * dt) + (0.5 * accelerations[p+11] * (dt * dt));
            if (positions[p+11] > box_dimensions){
                positions[p+11] -= box_dimensions; // if out left side of box, reenter at right side of box
            }
            if (positions[p+11] < 0){
                positions[p+11] += box_dimensions; // if out right side of box, reenter at left side of box
            }
            old_accelerations[p+10] = accelerations[p+10];
            old_accelerations[p+11] = accelerations[p+11];

            // particle 7
            positions[p+12] += (velocities[p+12] * dt) + (0.5 * accelerations[p+12] * (dt * dt));
            if (positions[p+12] > box_dimensions){
                positions[p+12] -= box_dimensions; // if out of top of box, reenter at bottom of box
            }
            if (positions[p+12] < 0){
                positions[p+12] += box_dimensions; // if out of bottom of box, reenter at top of box
            }
            positions[p+13] += (velocities[p+13] * dt) + (0.5 * accelerations[p+13] * (dt * dt));
            if (positions[p+13] > box_dimensions){
                positions[p+13] -= box_dimensions; // if out left side of box, reenter at right side of box
            }
            if (positions[p+13] < 0){
                positions[p+13] += box_dimensions; // if out right side of box, reenter at left side of box
            }
            old_accelerations[p+12] = accelerations[p+12];
            old_accelerations[p+13] = accelerations[p+13];

            // particle 8
            positions[p+14] += (velocities[p+14] * dt) + (0.5 * accelerations[p+14] * (dt * dt));
            if (positions[p+14] > box_dimensions){
                positions[p+14] -= box_dimensions; // if out of top of box, reenter at bottom of box
            }
            if (positions[p+14] < 0){
                positions[p+14] += box_dimensions; // if out of bottom of box, reenter at top of box
            }
            positions[p+15] += (velocities[p+15] * dt) + (0.5 * accelerations[p+15] * (dt * dt));
            if (positions[p+15] > box_dimensions){
                positions[p+15] -= box_dimensions; // if out left side of box, reenter at right side of box
            }
            if (positions[p+15] < 0){
                positions[p+15] += box_dimensions; // if out right side of box, reenter at left side of box
            }
            old_accelerations[p+14] = accelerations[p+14];
            old_accelerations[p+15] = accelerations[p+15];
            
        }
        // remaining values
        int rem = 8*(num_positions/8);
        for (int p = rem; p < num_positions; p = p + 2){
            positions[p] += (velocities[p] * dt) + (0.5 * accelerations[p] * (dt * dt));
            if (positions[p] > box_dimensions){
                positions[p] -= box_dimensions; // if out of top of box, reenter at bottom of box
            }
            if (positions[p] < 0){
                positions[p] += box_dimensions; // if out of bottom of box, reenter at top of box
            }
            positions[p+1] += (velocities[p+1] * dt) + (0.5 * accelerations[p+1] * (dt * dt));
            if (positions[p+1] > box_dimensions){
                positions[p+1] -= box_dimensions; // if out left side of box, reenter at right side of box
            }
            if (positions[p+1] < 0){
                positions[p+1] += box_dimensions; // if out right side of box, reenter at left side of box
            }
            old_accelerations[p] = accelerations[p];
            old_accelerations[p+1] = accelerations[p+1];    
        }


        // update accelerations according to current positions and all interactions
        for (int p = 0; p < num_positions; p = p + 2){
            double force_x = 0;
            double force_y = 0;
            particle_one[0] = positions[p];
            particle_one[1] = positions[p+1];
            for (int i = 0; i < num_positions - 8; i = i + 8){
                // particle 1
                if (i != p){
                    particle_two[0] = positions[i];
                    particle_two[1] = positions[i+1];
                    if (!distance_checks(&particle_one[0], &particle_two[0], &distances[0], box_dimensions)){
                        lennard_jones(&particle_one[0], &particle_two[0], &distances[0], &forces[0], &total_forces);
                        force_x += (forces[0]); // removed division by mass
                        force_y += (forces[1]); // removed division by mass
                    } 
                    distances[0] = 0;
                    distances[1] = 0;
                }

                // particle 2
                if (i + 2 != p){
                    particle_two[0] = positions[i+2];
                    particle_two[1] = positions[i+3];
                    if (!distance_checks(&particle_one[0], &particle_two[0], &distances[0], box_dimensions)){
                        lennard_jones(&particle_one[0], &particle_two[0], &distances[0], &forces[0], &total_forces);
                        force_x += (forces[0]); // removed division by mass
                        force_y += (forces[1]); // removed division by mass
                    } 
                    distances[0] = 0;
                    distances[1] = 0;
                }

                // particle 3
                if (i + 4 != p){
                    particle_two[0] = positions[i+4];
                    particle_two[1] = positions[i+5];
                    if (!distance_checks(&particle_one[0], &particle_two[0], &distances[0], box_dimensions)){
                        lennard_jones(&particle_one[0], &particle_two[0], &distances[0], &forces[0], &total_forces);
                        force_x += (forces[0]); // removed division by mass
                        force_y += (forces[1]); // removed division by mass
                    } 
                    distances[0] = 0;
                    distances[1] = 0;
                }

                // particle 4
                if (i + 6 != p){
                    particle_two[0] = positions[i+6];
                    particle_two[1] = positions[i+7];
                    if (!distance_checks(&particle_one[0], &particle_two[0], &distances[0], box_dimensions)){
                        lennard_jones(&particle_one[0], &particle_two[0], &distances[0], &forces[0], &total_forces);
                        force_x += (forces[0]); // removed division by mass
                        force_y += (forces[1]); // removed division by mass
                    } 
                    distances[0] = 0;
                    distances[1] = 0;
                }

                // particle 5
                if (i + 8 != p){
                    particle_two[0] = positions[i+8];
                    particle_two[1] = positions[i+9];
                    if (!distance_checks(&particle_one[0], &particle_two[0], &distances[0], box_dimensions)){
                        lennard_jones(&particle_one[0], &particle_two[0], &distances[0], &forces[0], &total_forces);
                        force_x += (forces[0]); // removed division by mass
                        force_y += (forces[1]); // removed division by mass
                    } 
                    distances[0] = 0;
                    distances[1] = 0;
                }

                // particle 6
                if (i + 10 != p){
                    particle_two[0] = positions[i+10];
                    particle_two[1] = positions[i+11];
                    if (!distance_checks(&particle_one[0], &particle_two[0], &distances[0], box_dimensions)){
                        lennard_jones(&particle_one[0], &particle_two[0], &distances[0], &forces[0], &total_forces);
                        force_x += (forces[0]); // removed division by mass
                        force_y += (forces[1]); // removed division by mass
                    } 
                    distances[0] = 0;
                    distances[1] = 0;
                }

                // particle 7
                if (i + 12 != p){
                    particle_two[0] = positions[i+12];
                    particle_two[1] = positions[i+13];
                    if (!distance_checks(&particle_one[0], &particle_two[0], &distances[0], box_dimensions)){
                        lennard_jones(&particle_one[0], &particle_two[0], &distances[0], &forces[0], &total_forces);
                        force_x += (forces[0]); // removed division by mass
                        force_y += (forces[1]); // removed division by mass
                    } 
                    distances[0] = 0;
                    distances[1] = 0;
                }

                // particle 8
                if (i + 14 != p){
                    particle_two[0] = positions[i+14];
                    particle_two[1] = positions[i+15];
                    if (!distance_checks(&particle_one[0], &particle_two[0], &distances[0], box_dimensions)){
                        lennard_jones(&particle_one[0], &particle_two[0], &distances[0], &forces[0], &total_forces);
                        force_x += (forces[0]); // removed division by mass
                        force_y += (forces[1]); // removed division by mass
                    } 
                    distances[0] = 0;
                    distances[1] = 0;
                }
            }
            // for particles that don't fit
            for (int i = rem; i < num_positions; i = i + 2){
                if (i != p){
                    particle_two[0] = positions[i];
                    particle_two[1] = positions[i+1];
                    if (!distance_checks(&particle_one[0], &particle_two[0], &distances[0], box_dimensions)){
                        lennard_jones(&particle_one[0], &particle_two[0], &distances[0], &forces[0], &total_forces);
                        force_x += (forces[0]); // removed division by mass
                        force_y += (forces[1]); // removed division by mass
                    } 
                    distances[0] = 0;
                    distances[1] = 0;
                }
            }
            accelerations[p] = force_x;
            accelerations[p+1] = force_y;
        }

        // update velocities
        for (int p = 0; p < (num_positions-8); p = p + 8){
            velocities[p] += 0.5 * (old_accelerations[p] + accelerations[p]) * dt;
            velocities[p+1] += 0.5 * (old_accelerations[p+1] + accelerations[p+1]) * dt;
            velocities[p+2] += 0.5 * (old_accelerations[p+2] + accelerations[p+2]) * dt;
            velocities[p+3] += 0.5 * (old_accelerations[p+3] + accelerations[p+3]) * dt;
            velocities[p+4] += 0.5 * (old_accelerations[p+4] + accelerations[p+4]) * dt;
            velocities[p+5] += 0.5 * (old_accelerations[p+5] + accelerations[p+5]) * dt;
            velocities[p+6] += 0.5 * (old_accelerations[p+6] + accelerations[p+6]) * dt;
            velocities[p+7] += 0.5 * (old_accelerations[p+7] + accelerations[p+7]) * dt;
            velocities[p+8] += 0.5 * (old_accelerations[p+8] + accelerations[p+8]) * dt;
            velocities[p+9] += 0.5 * (old_accelerations[p+9] + accelerations[p+9]) * dt;
            velocities[p+10] += 0.5 * (old_accelerations[p+10] + accelerations[p+10]) * dt;
            velocities[p+11] += 0.5 * (old_accelerations[p+11] + accelerations[p+11]) * dt;
            velocities[p+12] += 0.5 * (old_accelerations[p+12] + accelerations[p+12]) * dt;
            velocities[p+13] += 0.5 * (old_accelerations[p+13] + accelerations[p+13]) * dt;
            velocities[p+14] += 0.5 * (old_accelerations[p+14] + accelerations[p+14]) * dt;
            velocities[p+15] += 0.5 * (old_accelerations[p+15] + accelerations[p+15]) * dt;
        }
        // remaining particles
        for (int p = rem; p < num_positions; p = p + 2){
            velocities[p] += 0.5 * (old_accelerations[p] + accelerations[p]) * dt;
            velocities[p+1] += 0.5 * (old_accelerations[p+1] + accelerations[p+1]) * dt;
        }
        total_forces = 0;

    }

    // end time
    clock_t toc = clock();
    double time_spent = (double)(toc - tic) / CLOCKS_PER_SEC;
    printf("time taken for %d particles with optimised = %f\n", num_particles, time_spent);
    fflush(stdout);
    
    // //dump trajectory values to output
    // FILE* f = fopen("trajectory.data", "wb");
    // for (int i = 0; i < num_frames; i++){
    //     //printf("trj: %f", trajectory[0][i]);
    //     fwrite(trajectory[i], sizeof(double), num_positions, f);
    // }
    // fclose(f);

    // free(trajectory);
    free(velocities);
    free(positions);
    free(accelerations);
    free(masses);
    free(old_accelerations);

    return 0;
}