/**
 * @file main.c
 * @brief Contains the main program.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "headers.h"

/**
 * @brief The main entry point of the program.
 * 
 * This program manages reading of the input, performing the MD simualtion and generation of output.
 * 
 * @return int Returns 0 upon successful execution.
 */

int main(int argc, char *argv[]) {
    // Default settings for number of simulation steps and time step  
    int n_steps = 1000; //! Number of simulation steps
    double dt = 0.2; //! Time step

    // Check which command line options are provided
    if (argc != 2) {
        for (int i = 1; i < argc; i++) { // Loop over command line arguments
            if (strcmp(argv[i], "-n") == 0) {
                if (i + 1 < argc) {
                    n_steps = atoi(argv[i + 1]);
                }
                else {
                    fprintf(stderr, "Option -n requires the specification of the number of steps, e.g. -n 2000");
                    return 1;
                }
            }
            if (strcmp(argv[i], "-t") == 0) {
                if (i + 1 < argc) {
                    dt = atoi(argv[i + 1]);
                }
                else {
                    fprintf(stderr, "Option -t requires the specification of the timestep, e.g. dt = 0.1"); 
                    return 1;
                }
            }
        }
        if (argc == 1) {
            fprintf(stderr, "Usage: %s <filename>\n", argv[0]);
            return 1;
        }
    }

    const char* filename = argv[1]; //! Name of the input file
    
    // Read number of atoms
    int n_atoms = read_natoms(filename); //! Number of atoms

    // Allocate arrays
    double** coords = allocate_2d_array(n_atoms, 3); //! 2D array of coordinates consisting of x, y, z for each atom
    double** distances = allocate_2d_array(n_atoms, n_atoms); //! 2D array of distances for each pair of atoms 
    double* masses = (double*)malloc(n_atoms * sizeof(double)); //! Array of masses of each atom
    if (masses == NULL) {
        fprintf(stderr, "Memory allocation failed for masses!\n");
        free_2d_array(coords, n_atoms);
        return 1;
    }

    // Read coordinates and masses
    read_coords_and_masses(filename, coords, masses, n_atoms);

    // Validate masses
    double epsilon; //! Epsilon parameter for the Lennard-Jones potential in j/mol
    double sigma; //! Sigma parameter for the Lennard-Jones potential in nm
    if (!validate_atoms(masses, &epsilon, &sigma, n_atoms)) {
        free_2d_array(coords, n_atoms);
        free(masses);
        return 1;
    }

    // Allocate array for velocities and initialize them to zero
    double** velocities = allocate_2d_array(n_atoms, 3); //! 2D array of velocities consisting of vx, vy, vz for each atom
    for (int i = 0; i < n_atoms; i++) {
        for (int j = 0; j < 3; j++) {
            velocities[i][j] = 0.0;
        }
    }
    // Allocate array for accelerations and initialize them to zero
    double** accelerations = allocate_2d_array(n_atoms, 3); //! 2D array of accelerations consiting of ax, ay, az for each atom
    for (int i = 0; i < n_atoms; i++) {
        for (int j = 0; j < 3; j++) {
            accelerations[i][j] = 0.0;
        }
    }

    // Open files for writing the output
    const char* trajectory_name = "trajectory.xyz"; //! Name of the file where the trajectory output is written
    FILE* trajectory_file = open_output(trajectory_name); //! File where the trajectory output is written
    const char* energy_name = "energies"; //! Name of the file where the energies are written
    FILE* energy_file = open_output(energy_name); //! File where the energies are written
    const char* extended_name = "trajectory_velocity.xyz"; //! Name of the file where the extended trajectory with velocities is written
    FILE* extended_file = open_output(extended_name); //! File where the extended trajectory with velocities is written
    const char* acceleration_name = "acceleration"; //! Name of the file where the accelerations are written
    FILE* acceleration_file = open_output(acceleration_name); //! File where the accelerations are written

    // Run 1000 steps of MD simulation
    
    double kinetic_energy; //! Variable for storing the kinetic energy
    double potential_energy; //! Variable for storing the potential energy
    double total_energy; //! Variable for storing the total energy
    double previous_energy; //! Variable for storing the total energy of the previous step

    // Initialize timing variables
    clock_t start_md, end_md;
    double total_md_time;

    start_md = clock(); // Start timing the MD simulation

    for (int i =0; i < n_steps; i++){

        // Update positions, velocities and accelerations
        update_positions(coords, velocities, accelerations, dt, n_atoms);
        update_velocities(velocities, accelerations, dt, n_atoms); // First velocity update with old accelerations
        calculate_accelerations(coords, masses, n_atoms, epsilon, sigma, distances, accelerations);
        update_velocities(velocities, accelerations, dt, n_atoms); // Second velocity update with new accelerations
        
        // Calculate energies
        kinetic_energy = calculate_kinetic_energy(velocities, masses, n_atoms);
        potential_energy = calculate_potential_energy(distances, n_atoms, epsilon, sigma);
        previous_energy = total_energy;
        total_energy = calculate_total_energy(kinetic_energy, potential_energy);
        
        // Check if the energy is conserved or varies by more than 10 %
        check_energy(previous_energy, total_energy, i);

        // Print output
        print_output(trajectory_file, energy_file, extended_file, acceleration_file, 
                     n_atoms, i, kinetic_energy, potential_energy, total_energy,
                     coords, velocities, accelerations);       
    }
    
    end_md = clock(); // End timing the MD simulation

    // Timing statistics
    total_md_time = ((double) (end_md - start_md)) / CLOCKS_PER_SEC;
    double average_step_time = total_md_time / n_steps;

    printf("\n################# Timing Information ################\n");
    printf("Total number of steps:          %d\n", n_steps);
    printf("Total MD simulation time:       %.6f seconds\n", total_md_time);
    printf("Average time per step:          %.6f seconds\n", average_step_time);

    // Close output files
    fclose(trajectory_file);  
    fclose(extended_file);
    fclose(acceleration_file);

    // Free the allocated memory
    free_2d_array(coords, n_atoms);
    free_2d_array(distances, n_atoms);
    free_2d_array(velocities, n_atoms);
    free_2d_array(accelerations, n_atoms);
    free(masses);

    return 0;
}
