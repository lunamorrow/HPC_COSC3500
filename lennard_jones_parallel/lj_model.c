#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <omp.h> //openMP library

#include "lennard_jones.h"

#define NUM_BOXES 36
#define NUM_MPI 4
#define RADIUS 2.5
// #define NUM_THREADS 112

int seed = 5;

// divide simulation box in 4

//     0  1  2  | 9  10 11
//     3  4  5  | 12 13 14
//     6  7  8  | 15 16 17
//     -------------------
//     18 19 20 | 27 28 29
//     21 22 23 | 30 31 32
//     24 25 26 | 33 34 35 

// break particles into 4 arrays

// duplicate neightbours in border/buffer region -> each process considers
// particles in its grid + 'ghost' particles in the buffer

// update particles that move between processes or in/out of buffer regions as required

// if MPI process sees particles have left its grid and entered its buffer,
// send these particles to other processes as required

// recieving process adds incoming atoms to its arrays

int distance_checks(double* particle1, double* particle2, double* distances, double box_dimensions) {
    double x_dist;
    double y_dist;
    double x_force = 1;
    double y_force = 1;

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

void make_particles(double* velocities, double* positions, double* accelerations, double* masses, int num_positions, double* output, double dimension, double max_x, double max_y, int processNum){
    double centre_of_mass_velocity_total_x = 0;
    double centre_of_mass_velocity_total_y = 0;
    double total_mass = 0;
    /*
    p2 | p3
    -------
    p0 | p1

    p0 = 0,0
    p1 = dimension, 0
    p2 = 0, dimension
    p3 = dimenion, dimension
    */
    double x_position = 0.0;
    double y_position = 0.0;
    if (processNum == 1 || processNum == 3){
        x_position = dimension;
    } 
    if (processNum == 2 || processNum == 3){
        y_position = dimension;
    }

    double init_x_pos = *(&x_position); // copy initial x position
   
    // assist with random double generation
    double velocity_rand = ((double) 4 / RAND_MAX);
    double acceleration_rand = ((double) 2 / RAND_MAX);

    for (int i = 0; i < num_positions; i+=2) {
        velocities[i] = (rand() * velocity_rand) + (-2); // x velocity 
        velocities[i+1] = (rand() * velocity_rand) + (-2); // y velocity
        positions[i] = x_position;
        positions[i+1] = y_position;
        accelerations[i] = (rand() * acceleration_rand) + 1; // x acceleration
        accelerations[i+1] = (rand() * acceleration_rand) + 1;// y acceleration
        masses[i] = 1;
        masses[i+1] = 1;

        x_position += 2.5;
        if (x_position >= (max_x - 0.5)){
            x_position = init_x_pos;
            y_position += 2.5;
            if (y_position > (max_y - 1)){
                break;
            }
        }

        centre_of_mass_velocity_total_y += (masses[i+1] * velocities[i+1]);
        centre_of_mass_velocity_total_x += (masses[i] * velocities[i]);
        total_mass += (masses[i]);

        output[0] += centre_of_mass_velocity_total_x;
        output[1] += centre_of_mass_velocity_total_y;
        output[2] += total_mass;
    }
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

void share_border_particles(double* positions, double* myBorderParticles, int numPositions, double dimensions, int processNum){
    // each process checks which of its particles are along its outer 2.5 and sends them to the respective other processes to create their buffer regions
    double max_x = 0.0;
    double max_y = 0.0;
    double min_x = 0.0;
    double min_y = 0.0;

    double* forBufferHorizontal_positions = (double*) malloc(sizeof(double) * numPositions);
    double* forBufferVertical_positions = (double*) malloc(sizeof(double) * numPositions);

    memset(forBufferHorizontal_positions, 0, sizeof(double) * numPositions);
    memset(forBufferVertical_positions, 0, sizeof(double) * numPositions);

    if (processNum == 1 || processNum == 3){
        min_x = dimensions;
        max_x = dimensions + dimensions;
    } 
    if (processNum == 2 || processNum == 3){
        min_y = dimensions;
        max_y = dimensions + dimensions;
    }
    if (processNum == 0 || processNum == 2){
        max_x = dimensions;
    }
    if (processNum == 0 || processNum == 1){
        max_y = dimensions;
    }

    double x, y;
    int vInd = 0;
    int hInd = 0;
    int hAdded = 0;
    int vAdded = 0;
    for (int i = 0; i < numPositions; i+=2){
        x = positions[i];
        y = positions[i+1];
        
        if (x >= max_x - RADIUS || ((x <= min_x + RADIUS) && x > -90)){
            forBufferHorizontal_positions[hInd] = positions[i];
            forBufferHorizontal_positions[hInd+1] = positions[i+1];
            hAdded++;
            hInd+=2;
        }
        if (y >= max_y - RADIUS || ((y <= min_y + RADIUS) && y > -90)){
            forBufferVertical_positions[vInd] = positions[i];
            forBufferVertical_positions[vInd+1] = positions[i+1];
            vAdded++;
            vInd+=2;
        }
    }

    if (processNum == 0){
        // send sizes
        int nRecv0, nRecv1;
        MPI_Request requests[2];
        MPI_Irecv(&nRecv0, 1, MPI_INT, 2, 10, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(&nRecv1, 1, MPI_INT, 1, 10, MPI_COMM_WORLD, &requests[1]);
        MPI_Send(&vAdded, 1, MPI_INT, 2, 10, MPI_COMM_WORLD);
        MPI_Send(&hAdded, 1, MPI_INT, 1, 10, MPI_COMM_WORLD);
        MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

        // must recieve vertical from 2 and horizontal from 1
        // make recieve buffer size large to be safe and as don't know how many particles other process is sending
        MPI_Request rrequest[2];
        MPI_Status stat0, stat1, stat2, stat3;

        // FLAGS: 0 = positions; 1 = accelerations; 2 = velocities; 3 = masses
        // save positions
        MPI_Irecv(&myBorderParticles[0], nRecv0, MPI_DOUBLE, 2, 0,
              MPI_COMM_WORLD, &rrequest[0]);
        MPI_Irecv(&myBorderParticles[nRecv0], nRecv1, MPI_DOUBLE, 1, 0,
              MPI_COMM_WORLD, &rrequest[1]);
        
        // vertical is 2 and horizontal is 1
        // FLAGS: 0 = positions; 1 = accelerations; 2 = velocities; 3 = masses
        MPI_Send(forBufferVertical_positions, vAdded, MPI_DOUBLE, 2, 0,
              MPI_COMM_WORLD);
        
        MPI_Send(forBufferHorizontal_positions, hAdded, MPI_DOUBLE, 1, 0,
              MPI_COMM_WORLD);
        MPI_Waitall(2, rrequest, MPI_STATUSES_IGNORE);
    } else if (processNum == 1){
        // send sizes
        int nRecv0, nRecv1;
        MPI_Request requests[2];
        MPI_Irecv(&nRecv0, 1, MPI_INT, 3, 10, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(&nRecv1, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &requests[1]);
        MPI_Send(&vAdded, 1, MPI_INT, 3, 10, MPI_COMM_WORLD);
        MPI_Send(&hAdded, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
        MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

        // must recieve vertical from 3 and horizontal from 0
        // make recieve buffer size large to be safe and as don't know how many particles other process is sending
        MPI_Request rrequest[2];
        MPI_Status stat0, stat1, stat2, stat3;

        // FLAGS: 0 = positions; 1 = accelerations; 2 = velocities; 3 = masses
        // save positions
        MPI_Irecv(&myBorderParticles[0], nRecv0, MPI_DOUBLE, 3, 0,
              MPI_COMM_WORLD, &rrequest[0]);
        MPI_Irecv(&myBorderParticles[nRecv0], nRecv1, MPI_DOUBLE, 0, 0,
              MPI_COMM_WORLD, &rrequest[1]);
      
        // vertical is 3 and horizontal is 0
        // FLAGS: 0 = positions; 1 = accelerations; 2 = velocities; 3 = masses
        MPI_Send(forBufferVertical_positions, vAdded, MPI_DOUBLE, 3, 0,
              MPI_COMM_WORLD);

        MPI_Send(forBufferHorizontal_positions, hAdded, MPI_DOUBLE, 0, 0,
              MPI_COMM_WORLD);
        MPI_Waitall(2, rrequest, MPI_STATUSES_IGNORE);
    } else if (processNum == 2){
        // send sizes
        int nRecv0, nRecv1;
        MPI_Request sizereq_0, sizereq_1;
        MPI_Request requests[2];
        MPI_Irecv(&nRecv0, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(&nRecv1, 1, MPI_INT, 3, 10, MPI_COMM_WORLD, &requests[1]);
        MPI_Send(&vAdded, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
        MPI_Send(&hAdded, 1, MPI_INT, 3, 10, MPI_COMM_WORLD);
        MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

        // must recieve vertical from 0 and horizontal from 3
        // make recieve buffer size large to be safe and as don't know how many particles other process is sending
        MPI_Request rrequest[2];
        MPI_Status stat0, stat1, stat2, stat3;

        // FLAGS: 0 = positions; 1 = accelerations; 2 = velocities; 3 = masses
        // save positions
        MPI_Irecv(&myBorderParticles[0], nRecv0, MPI_DOUBLE, 0, 0,
              MPI_COMM_WORLD, &rrequest[0]);
        MPI_Irecv(&myBorderParticles[nRecv0], nRecv1, MPI_DOUBLE, 3, 0,
              MPI_COMM_WORLD, &rrequest[1]);
        
        // vertical is 0 and horizontal is 3
        // FLAGS: 0 = positions; 1 = accelerations; 2 = velocities; 3 = masses
        MPI_Send(forBufferVertical_positions, vAdded, MPI_DOUBLE, 0, 0,
              MPI_COMM_WORLD);
        
        MPI_Send(forBufferHorizontal_positions, hAdded, MPI_DOUBLE, 3, 0,
              MPI_COMM_WORLD);
        MPI_Waitall(2, rrequest, MPI_STATUSES_IGNORE);
    } else if (processNum == 3){
        // send sizes
        int nRecv0, nRecv1;
        MPI_Request sizereq_0, sizereq_1;
        MPI_Request requests[2];
        MPI_Irecv(&nRecv0, 1, MPI_INT, 1, 10, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(&nRecv1, 1, MPI_INT, 2, 10, MPI_COMM_WORLD, &requests[1]);
        MPI_Send(&vAdded, 1, MPI_INT, 1, 10, MPI_COMM_WORLD);
        MPI_Send(&hAdded, 1, MPI_INT, 2, 10, MPI_COMM_WORLD);
        MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

        // must recieve vertical from 1 and horizontal from 2
        // make recieve buffer size large to be safe and as don't know how many particles other process is sending
        MPI_Request rrequest[2];
        MPI_Status stat0, stat1, stat2, stat3;

        // FLAGS: 0 = positions; 1 = accelerations; 2 = velocities; 3 = masses
        // save positions
        MPI_Irecv(&myBorderParticles[0], nRecv0, MPI_DOUBLE, 1, 0,
              MPI_COMM_WORLD, &rrequest[0]);
        MPI_Irecv(&myBorderParticles[nRecv0], nRecv1, MPI_DOUBLE, 2, 0,
              MPI_COMM_WORLD, &rrequest[1]);

        // vertical is 1 and horizontal is 2
        // FLAGS: 0 = positions; 1 = accelerations; 2 = velocities; 3 = masses
        MPI_Send(forBufferVertical_positions, vAdded, MPI_DOUBLE, 1, 0,
              MPI_COMM_WORLD);
        
        MPI_Send(forBufferHorizontal_positions, hAdded, MPI_DOUBLE, 2, 0,
              MPI_COMM_WORLD);
        MPI_Waitall(2, rrequest, MPI_STATUSES_IGNORE);
    }
    // MPI_Barrier(MPI_COMM_WORLD);
    free(forBufferHorizontal_positions);
    free(forBufferVertical_positions);
    return;
}

int transfer_particles(int processNum, double max_x, double max_y, double min_x, double min_y, double* positions, double* velocities, double* accelerations, double* masses, int myNumParticles, int dimensions, double* forBufferHorizontal_positions, double* forBufferVertical_positions, double* forBufferHorizontal_velocities, double* forBufferVertical_velocities, double* forBufferHorizontal_accelerations, double* forBufferVertical_accelerations, double* forBufferHorizontal_masses, double* forBufferVertical_masses){
    // each process checks which of its particles are outside of its bounds and sends to other particles

    double x, y;
    int vInd = 0;
    int hInd = 0;
    int hAdded = 0;
    int vAdded = 0;
    for (int i = 0; i < myNumParticles; i+=2){
        x = positions[i];
        y = positions[i+1];
        if (x > max_x || (x < min_x && x > -90)){
            forBufferHorizontal_positions[hInd] = positions[i];
            forBufferHorizontal_positions[hInd+1] = positions[i+1];
            forBufferHorizontal_accelerations[hInd] = accelerations[i];
            forBufferHorizontal_accelerations[hInd+1] = accelerations[i+1];
            forBufferHorizontal_velocities[hInd] = velocities[i];
            forBufferHorizontal_velocities[hInd+1] = velocities[i+1];
            forBufferHorizontal_masses[hInd] = 1;
            forBufferHorizontal_masses[hInd+1] = 1;
            // remove particles
            positions[i] = (double)-100;
            positions[i+1] = (double)-100;
            hAdded+=2;
            hInd+=2;
        } else if (y > max_y || (y < min_y && y > -90)){
            forBufferVertical_positions[vInd] = positions[i];
            forBufferVertical_positions[vInd+1] = positions[i+1];
            forBufferVertical_accelerations[vInd] = accelerations[i];
            forBufferVertical_accelerations[vInd+1] = accelerations[i+1];
            forBufferVertical_velocities[vInd] = velocities[i];
            forBufferVertical_velocities[vInd+1] = velocities[i+1];
            forBufferVertical_masses[vInd] = 1;
            forBufferVertical_masses[vInd+1] = 1;
            // remove particles
            positions[i] = (double)-100;
            positions[i+1] = (double)-100;
            vAdded+=2;
            vInd+=2;
        }
    }

    int nRecv0 = 0;
    int nRecv1 = 0;

    if (processNum == 0){
        // send and recieve sizes
        MPI_Request requests[2];
        MPI_Irecv(&nRecv0, 1, MPI_INT, 2, 10, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(&nRecv1, 1, MPI_INT, 1, 10, MPI_COMM_WORLD, &requests[1]);
        MPI_Send(&vAdded, 1, MPI_INT, 2, 10, MPI_COMM_WORLD);
        MPI_Send(&hAdded, 1, MPI_INT, 1, 10, MPI_COMM_WORLD);
        MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

        // must recieve vertical from 2 and horizontal from 1
        MPI_Request rrequest_0, rrequest_1, rrequest_2, rrequest_3, rrequest_4, rrequest_5, rrequest_6, rrequest_7;
        // FLAGS: 0 = positions; 1 = accelerations; 2 = velocities; 3 = masses
        MPI_Irecv(&positions[myNumParticles], nRecv0, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, &rrequest_0);
        MPI_Irecv(&positions[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &rrequest_4);
        
        MPI_Irecv(&accelerations[myNumParticles], nRecv0, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD, &rrequest_1);
        MPI_Irecv(&accelerations[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &rrequest_5);
        
        MPI_Irecv(&velocities[myNumParticles], nRecv0, MPI_DOUBLE, 2, 2, MPI_COMM_WORLD, &rrequest_2);
        MPI_Irecv(&velocities[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD, &rrequest_6);
        
        MPI_Irecv(&masses[myNumParticles], nRecv0, MPI_DOUBLE, 2, 3, MPI_COMM_WORLD, &rrequest_3);
        MPI_Irecv(&masses[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 1, 3, MPI_COMM_WORLD, &rrequest_7);

        // send vertical to 2 and horizontal to 1
        MPI_Request request_0, request_1, request_2, request_3, request_4, request_5, request_6, request_7;
        MPI_Send(forBufferVertical_positions, vAdded, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_positions, hAdded, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);

        MPI_Send(forBufferVertical_accelerations, vAdded, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_accelerations, hAdded, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
        
        MPI_Send(forBufferVertical_velocities, vAdded, MPI_DOUBLE, 2, 2, MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_velocities, hAdded, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD);

        MPI_Send(forBufferVertical_masses, vAdded, MPI_DOUBLE, 2, 3, MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_masses, hAdded, MPI_DOUBLE, 1, 3, MPI_COMM_WORLD);
    } else if (processNum == 1){
        // send and recieve sizes
        MPI_Request requests[2];
        MPI_Irecv(&nRecv0, 1, MPI_INT, 3, 10, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(&nRecv1, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &requests[1]);
        MPI_Send(&vAdded, 1, MPI_INT, 3, 10, MPI_COMM_WORLD);
        MPI_Send(&hAdded, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
        MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

        // must recieve vertical from 3 and horizontal from 0
        // make recieve buffer size large to be safe and as don't know how many particles other process is sending
        MPI_Request rrequest_0, rrequest_1, rrequest_2, rrequest_3, rrequest_4, rrequest_5, rrequest_6, rrequest_7;

        // FLAGS: 0 = positions; 1 = accelerations; 2 = velocities; 3 = masses
        // save positions
        MPI_Irecv(&positions[myNumParticles], nRecv0, MPI_DOUBLE, 3, 0,
              MPI_COMM_WORLD, &rrequest_0);
        MPI_Irecv(&positions[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 0, 0,
              MPI_COMM_WORLD, &rrequest_4);
        
        MPI_Irecv(&accelerations[myNumParticles], nRecv0, MPI_DOUBLE, 3, 1,
              MPI_COMM_WORLD, &rrequest_1);
        MPI_Irecv(&accelerations[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 0, 1,
              MPI_COMM_WORLD, &rrequest_5);
        
        MPI_Irecv(&velocities[myNumParticles], nRecv0, MPI_DOUBLE, 3, 2,
              MPI_COMM_WORLD, &rrequest_2);
        MPI_Irecv(&velocities[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 0, 2,
              MPI_COMM_WORLD, &rrequest_6);
        
        MPI_Irecv(&masses[myNumParticles], nRecv0, MPI_DOUBLE, 3, 3,
              MPI_COMM_WORLD, &rrequest_3);
        MPI_Irecv(&masses[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 0, 3,
              MPI_COMM_WORLD, &rrequest_7);
        //MPI_Wait(&rrequest_7, MPI_STATUS_IGNORE);

        // vertical is 3 and horizontal is 0
        MPI_Request request_0, request_1, request_2, request_3, request_4, request_5, request_6, request_7;
        // FLAGS: 0 = positions; 1 = accelerations; 2 = velocities; 3 = masses
        MPI_Send(forBufferVertical_positions, vAdded, MPI_DOUBLE, 3, 0,
              MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_positions, hAdded, MPI_DOUBLE, 0, 0,
              MPI_COMM_WORLD);

        MPI_Send(forBufferVertical_accelerations, vAdded, MPI_DOUBLE, 3, 1,
              MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_accelerations, hAdded, MPI_DOUBLE, 0, 1,
              MPI_COMM_WORLD);
        
        MPI_Send(forBufferVertical_velocities, vAdded, MPI_DOUBLE, 3, 2,
              MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_velocities, hAdded, MPI_DOUBLE, 0, 2,
              MPI_COMM_WORLD);

        MPI_Send(forBufferVertical_masses, vAdded, MPI_DOUBLE, 3, 3,
              MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_masses, hAdded, MPI_DOUBLE, 0, 3,
              MPI_COMM_WORLD);
    } else if (processNum == 2){
        // send and recieve sizes
        MPI_Request requests[2];
        MPI_Irecv(&nRecv0, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(&nRecv1, 1, MPI_INT, 3, 10, MPI_COMM_WORLD, &requests[1]);
        MPI_Send(&vAdded, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
        MPI_Send(&hAdded, 1, MPI_INT, 3, 10, MPI_COMM_WORLD);
        MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

        // must recieve vertical from 0 and horizontal from 3
        // make recieve buffer size large to be safe and as don't know how many particles other process is sending
        MPI_Request rrequest_0, rrequest_1, rrequest_2, rrequest_3, rrequest_4, rrequest_5, rrequest_6, rrequest_7;

        // FLAGS: 0 = positions; 1 = accelerations; 2 = velocities; 3 = masses
        // save positions        
        MPI_Irecv(&positions[myNumParticles], nRecv0, MPI_DOUBLE, 0, 0,
              MPI_COMM_WORLD, &rrequest_0);
        MPI_Irecv(&positions[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 3, 0,
              MPI_COMM_WORLD, &rrequest_4);
        
        MPI_Irecv(&accelerations[myNumParticles], nRecv0, MPI_DOUBLE, 0, 1,
              MPI_COMM_WORLD, &rrequest_1);
        MPI_Irecv(&accelerations[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 3, 1,
              MPI_COMM_WORLD, &rrequest_5);
        
        MPI_Irecv(&velocities[myNumParticles], nRecv0, MPI_DOUBLE, 0, 2,
              MPI_COMM_WORLD, &rrequest_2);
        MPI_Irecv(&velocities[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 3, 2,
              MPI_COMM_WORLD, &rrequest_6);
        
        MPI_Irecv(&masses[myNumParticles], nRecv0, MPI_DOUBLE, 0, 3,
              MPI_COMM_WORLD, &rrequest_3);
        MPI_Irecv(&masses[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 3, 3,
              MPI_COMM_WORLD, &rrequest_7);
        //MPI_Wait(&rrequest_7, MPI_STATUS_IGNORE);

        // vertical is 0 and horizontal is 3
        MPI_Request request_0, request_1, request_2, request_3, request_4, request_5, request_6, request_7;
        // FLAGS: 0 = positions; 1 = accelerations; 2 = velocities; 3 = masses
        MPI_Send(forBufferVertical_positions, vAdded, MPI_DOUBLE, 0, 0,
              MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_positions, hAdded, MPI_DOUBLE, 3, 0,
              MPI_COMM_WORLD);

        MPI_Send(forBufferVertical_accelerations, vAdded, MPI_DOUBLE, 0, 1,
              MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_accelerations, hAdded, MPI_DOUBLE, 3, 1,
              MPI_COMM_WORLD);
        
        MPI_Send(forBufferVertical_velocities, vAdded, MPI_DOUBLE, 0, 2,
              MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_velocities, hAdded, MPI_DOUBLE, 3, 2,
              MPI_COMM_WORLD);

        MPI_Send(forBufferVertical_masses, vAdded, MPI_DOUBLE, 0, 3,
              MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_masses, hAdded, MPI_DOUBLE, 3, 3,
              MPI_COMM_WORLD);
    } else if (processNum == 3){
        // send and recieve sizes
        MPI_Request requests[2];
        MPI_Irecv(&nRecv0, 1, MPI_INT, 1, 10, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(&nRecv1, 1, MPI_INT, 2, 10, MPI_COMM_WORLD, &requests[1]);
        MPI_Send(&vAdded, 1, MPI_INT, 1, 10, MPI_COMM_WORLD);
        MPI_Send(&hAdded, 1, MPI_INT, 2, 10, MPI_COMM_WORLD);
        MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

        // must recieve vertical from 1 and horizontal from 2
        // make recieve buffer size large to be safe and as don't know how many particles other process is sending
        MPI_Request rrequest_0, rrequest_1, rrequest_2, rrequest_3, rrequest_4, rrequest_5, rrequest_6, rrequest_7;

        // FLAGS: 0 = positions; 1 = accelerations; 2 = velocities; 3 = masses
        // save positions
        MPI_Irecv(&positions[myNumParticles], nRecv0, MPI_DOUBLE, 1, 0,
              MPI_COMM_WORLD, &rrequest_0);
        MPI_Irecv(&positions[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 2, 0,
              MPI_COMM_WORLD, &rrequest_4);
        
        MPI_Irecv(&accelerations[myNumParticles], nRecv0, MPI_DOUBLE, 1, 1,
              MPI_COMM_WORLD, &rrequest_1);
        MPI_Irecv(&accelerations[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 2, 1,
              MPI_COMM_WORLD, &rrequest_5);
        
        MPI_Irecv(&velocities[myNumParticles], nRecv0, MPI_DOUBLE, 1, 2,
              MPI_COMM_WORLD, &rrequest_2);
        MPI_Irecv(&velocities[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 2, 2,
              MPI_COMM_WORLD, &rrequest_6);
        
        MPI_Irecv(&masses[myNumParticles], nRecv0, MPI_DOUBLE, 1, 3,
              MPI_COMM_WORLD, &rrequest_3);
        MPI_Irecv(&masses[myNumParticles + nRecv0], nRecv1, MPI_DOUBLE, 2, 3,
              MPI_COMM_WORLD, &rrequest_7);
        //MPI_Wait(&rrequest_7, MPI_STATUS_IGNORE);


        // vertical is 1 and horizontal is 2
        MPI_Request request_0, request_1, request_2, request_3, request_4, request_5, request_6, request_7;
        // FLAGS: 0 = positions; 1 = accelerations; 2 = velocities; 3 = masses
        MPI_Send(forBufferVertical_positions, vAdded, MPI_DOUBLE, 1, 0,
              MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_positions, hAdded, MPI_DOUBLE, 2, 0,
              MPI_COMM_WORLD);

        MPI_Send(forBufferVertical_accelerations, vAdded, MPI_DOUBLE, 1, 1,
              MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_accelerations, hAdded, MPI_DOUBLE, 2, 1,
              MPI_COMM_WORLD);
        
        MPI_Send(forBufferVertical_velocities, vAdded, MPI_DOUBLE, 1, 2,
              MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_velocities, hAdded, MPI_DOUBLE, 2, 2,
              MPI_COMM_WORLD);

        MPI_Send(forBufferVertical_masses, vAdded, MPI_DOUBLE, 1, 3,
              MPI_COMM_WORLD);
        MPI_Send(forBufferHorizontal_masses, hAdded, MPI_DOUBLE, 2, 3,
              MPI_COMM_WORLD);

    }

    // MPI_Barrier(MPI_COMM_WORLD);
    return (int)(myNumParticles + nRecv0 + nRecv1); // don't minus as buffer shuffles
}

int main(int argv, char** argc){
    MPI_Init(NULL, NULL);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    srand(seed);

    double share_time = 0;

    const int num_particles = atoi(argc[1]);
    const float box_dimensions = atof(argc[2]);

    // work out how many particles in each process and their grids
    const int num_particles_local = num_particles/world_size;
    const double gridDimension = box_dimensions/ (double)2;
    double myUpperX;
    double myUpperY;
    double myLowerX;
    double myLowerY;
    if (my_rank == 0){
        myUpperX = gridDimension;
        myUpperY = gridDimension;
        myLowerX = 0.0;
        myLowerY = 0.0;
    } else if (my_rank == 1){
        myUpperX = box_dimensions;
        myUpperY = gridDimension;
        myLowerX = gridDimension;
        myLowerY = 0.0;
    } else if (my_rank == 2){
        myUpperX = gridDimension;
        myUpperY = box_dimensions;
        myLowerX = 0.0;
        myLowerY = gridDimension;
    } else if (my_rank == 3){
        myUpperX = box_dimensions;
        myUpperY = box_dimensions;
        myLowerX = gridDimension;
        myLowerY = gridDimension;
    }

    // array to store border particle positions in
    double* myBorderParticles = (double*) malloc(sizeof(double) * 1024); // (num_particles_local + num_particles_local + num_particles_local)
    memset(myBorderParticles, 0, sizeof(double) * 1024);

    // set values
    double dt = 0.000001;
    double simulation_length = dt*100000;
    int num_positions = num_particles + num_particles;
    int num_positions_local = num_positions/world_size;
    double num_frames = simulation_length/dt;

    // make particles -> extra large to ensure buffer and mixing/crowding won't cause a sigfault
    double* velocities = (double*) malloc(sizeof(double) * num_positions);
    double* positions = (double*) malloc(sizeof(double) * num_positions);
    double* accelerations = (double*) malloc(sizeof(double) * num_positions);
    double* masses = (double*) malloc(sizeof(double) * num_positions);
    memset(&velocities[0], 0, sizeof(double) * num_positions); 
    memset(&positions[0], 0, sizeof(double) * num_positions); 
    memset(&accelerations[0], 0, sizeof(double) * num_positions); 
    memset(&masses[0], 0, sizeof(double) * num_positions); 
    double* output = (double*) calloc(3, sizeof(double));

    make_particles(&velocities[0], &positions[0], &accelerations[0], &masses[0], num_positions_local, &output[0], gridDimension, myUpperX, myUpperY, my_rank);

    int myNumParticles = num_positions_local;
    int bufferSize = num_positions_local + num_positions_local;

    double* forBufferHorizontal_positions = (double*) malloc(sizeof(double) * bufferSize);
    double* forBufferVertical_positions = (double*) malloc(sizeof(double) * bufferSize);
    double* forBufferHorizontal_velocities = (double*) malloc(sizeof(double) * bufferSize);
    double* forBufferVertical_velocities = (double*) malloc(sizeof(double) * bufferSize);
    double* forBufferHorizontal_accelerations = (double*) malloc(sizeof(double) * bufferSize);
    double* forBufferVertical_accelerations = (double*) malloc(sizeof(double) * bufferSize);
    double* forBufferHorizontal_masses = (double*) malloc(sizeof(double) * bufferSize);
    double* forBufferVertical_masses = (double*) malloc(sizeof(double) * bufferSize);
    memset(forBufferHorizontal_positions, 0, sizeof(double) * bufferSize);
    memset(forBufferVertical_positions, 0, sizeof(double) * bufferSize);
    memset(forBufferHorizontal_velocities, 0, sizeof(double) * bufferSize);
    memset(forBufferVertical_velocities, 0, sizeof(double) * bufferSize);
    memset(forBufferHorizontal_accelerations, 0, sizeof(double) * bufferSize);
    memset(forBufferVertical_accelerations, 0, sizeof(double) * bufferSize);
    memset(forBufferHorizontal_masses, 0, sizeof(double) * bufferSize);
    memset(forBufferVertical_masses, 0, sizeof(double) * bufferSize);

    // set temperature
    double centre_of_mass_velocity_x = output[0]/output[2]; // centre_of_mass_velocity_total_x/total_mass
    double centre_of_mass_velocity_y = output[1]/output[2]; // centre_of_mass_velocity_total_y/total_mass

    double temperature = set_temperature(&velocities[0], &masses[0], centre_of_mass_velocity_x, centre_of_mass_velocity_y, num_particles, num_particles);

    // // create array to store trajectory and store initial values
    // double** trajectory = (double**) malloc(sizeof(double*) * num_frames); //MAKE 1D -> makes quick to send w CUDA
    // for (int ti = 0; ti < num_frames; ti++){
    //     trajectory[ti] = (double*) malloc(sizeof(double) * num_positions); //to store x and y positions
    // }

    double distances[2];
    double forces[2];
    double particle_one[2];
    double particle_two[2];
    double particle_three[2];
    double total_forces = 0;
    double* old_accelerations = (double*) malloc(sizeof(double) * num_positions * 10);
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
        memset(old_accelerations, 0, (sizeof(double) * num_positions * 10));
        int vInd = 0;
        int hInd = 0;
        for (int p = 0; p < (myNumParticles); p = p + 2){
            if ((positions[p] < -90) || (positions[p+1] < -90)){
                continue;
                // skip 'deleted' particles
            }
            positions[p] += (velocities[p] * dt) + (0.5 * accelerations[p] * (dt * dt));
            if (positions[p] > myUpperX){
                if (positions[p] > box_dimensions){
                    positions[p] -= box_dimensions; // if out of top of box, reenter at bottom of box
                }
            } else if (positions[p] < myLowerX){
                if (positions[p] < 0){
                    positions[p] += box_dimensions; // if out of bottom of box, reenter at top of box
                }
            }
            positions[p+1] += (velocities[p+1] * dt) + (0.5 * accelerations[p+1] * (dt * dt));
            if (positions[p+1] > myUpperY){
                if (positions[p+1] > box_dimensions){
                    positions[p+1] -= box_dimensions; // if out of right of box, reenter at left of box
                }
            } else if (positions[p+1] < myLowerY){
                if (positions[p+1] < 0){
                    positions[p+1] += box_dimensions; // if out of left of box, reenter at right of box
                }
            }
            old_accelerations[p] = accelerations[p];
            old_accelerations[p+1] = accelerations[p+1];
        }

        // transfer any particles that are now outside of the grid
        MPI_Barrier(MPI_COMM_WORLD);

        myNumParticles = transfer_particles(my_rank, myUpperX, myUpperY, myLowerX, myLowerY, &positions[0], &velocities[0], &accelerations[0], &masses[0], myNumParticles, gridDimension,  forBufferHorizontal_positions, forBufferVertical_positions, forBufferHorizontal_velocities, forBufferVertical_velocities, forBufferHorizontal_accelerations, forBufferVertical_accelerations, forBufferHorizontal_masses, forBufferVertical_masses);
        
        memset(forBufferHorizontal_positions, 0, sizeof(double) * bufferSize);
        memset(forBufferVertical_positions, 0, sizeof(double) * bufferSize);
        memset(forBufferHorizontal_velocities, 0, sizeof(double) * bufferSize);
        memset(forBufferVertical_velocities, 0, sizeof(double) * bufferSize);
        memset(forBufferHorizontal_accelerations, 0, sizeof(double) * bufferSize);
        memset(forBufferVertical_accelerations, 0, sizeof(double) * bufferSize);
        memset(forBufferHorizontal_masses, 0, sizeof(double) * bufferSize);
        memset(forBufferVertical_masses, 0, sizeof(double) * bufferSize);

        share_border_particles(&positions[0], &myBorderParticles[0], num_positions_local, gridDimension, my_rank);

        // MPI_Barrier(MPI_COMM_WORLD);
        // update accelerations according to current positions and all interactions
        // #pragma omp parallel for num_threads(4) schedule(dynamic) private(distances, particle_one)
        for (int p = 0; p < myNumParticles; p = p + 2){
            // printf("p%d: I am thread %d\n", my_rank, omp_get_thread_num());
            // fflush(stdout);
            double force_x = 0;
            double force_y = 0;
            particle_one[0] = positions[p];
            particle_one[1] = positions[p+1];
            //#pragma omp parallel shared(force_x, force_y)
            // {
            if ((particle_one[0] > -90) && (particle_one[1] > -90)){
                // for particles in grid
                for (int i = 0; i < myNumParticles; i += 2){
                    if (i != p){
                        particle_two[0] = positions[i];
                        particle_two[1] = positions[i+1];
                        if (!(particle_two[0] < -90) && !(particle_two[1] < -90) && !distance_checks(&particle_one[0], &particle_two[0], &distances[0], box_dimensions)){
                            lennard_jones(&particle_one[0], &particle_two[0], &distances[0], &forces[0], &total_forces);
                            force_x += (forces[0]); // removed division by mass
                            force_y += (forces[1]); // removed division by mass
                        } 
                        distances[0] = 0;
                        distances[1] = 0;
                    }
                }
                // for border particles
                for (int j = 0; j < (num_particles_local); j += 2){
                    particle_three[0] = myBorderParticles[j];
                    particle_three[1] = myBorderParticles[j+1];
                    if (particle_three[0] != 0 && particle_three[1] != 0){
                        if (!distance_checks(&particle_one[0], &particle_three[0], &distances[0], box_dimensions)){
                            lennard_jones(&particle_one[0], &particle_three[0], &distances[0], &forces[0], &total_forces);
                            force_x += (forces[0]); // removed division by mass
                            force_y += (forces[1]); // removed division by mass
                        } 
                        distances[0] = 0;
                        distances[1] = 0;
                    } else {
                        // if particle with position {0, 0} found, you're at the end
                        break;
                    }
                }
            }
            accelerations[p] = force_x;
            accelerations[p+1] = force_y;
        }
        // clear border particles before next iteration
        memset(myBorderParticles, 0, sizeof(double) * 1024);

        // update velocities
        for (int p = 0; p < myNumParticles; p += 2){
            if ((positions[p] < -90) || (positions[p+1] < -90)){
                continue;
                // skip 'deleted' particles
            }
            velocities[p] += 0.5 * (old_accelerations[p] + accelerations[p]) * dt;
            velocities[p+1] += 0.5 * (old_accelerations[p+1] + accelerations[p+1]) * dt;
        }

        total_forces = 0;
       // MPI_Barrier(MPI_COMM_WORLD);
    }

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
    free(forBufferHorizontal_positions);
    free(forBufferVertical_positions);
    free(forBufferHorizontal_velocities);
    free(forBufferVertical_velocities);
    free(forBufferHorizontal_accelerations);
    free(forBufferVertical_accelerations);
    free(forBufferHorizontal_masses);
    free(forBufferVertical_masses);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}