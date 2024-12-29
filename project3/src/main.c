#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "headers.h"

int main(int argc, char *argv[]) {
    // Check if filename is provided
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <filename>\n", argv[0]);
        return 1;
    }

    const char* filename = argv[1];
    
    // Read number of atoms
    int n_atoms = read_natoms(filename);
    //printf("Number of atoms: %d\n", n_atoms);

    // Allocate arrays
    double** coords = allocate_2d_array(n_atoms, 3);
    double** distances = allocate_2d_array(n_atoms, n_atoms);
    double* masses = (double*)malloc(n_atoms * sizeof(double));
    if (masses == NULL) {
        fprintf(stderr, "Memory allocation failed for masses!\n");
        free_2d_array(coords, n_atoms);
        return 1;
    }

    // Read coordinates and masses
    read_coords_and_masses(filename, coords, masses, n_atoms);

    // Validate masses
    double epsilon, sigma;
    if (!validate_atoms(masses, &epsilon, &sigma, n_atoms)) {
        free_2d_array(coords, n_atoms);
        free(masses);
        return 1;
    }

    // Allocate array for velocities and initialize them to zero
    double** velocities = allocate_2d_array(n_atoms, 3);
    for (int i = 0; i < n_atoms; i++) {
        for (int j = 0; j < 3; j++) {
            velocities[i][j] = 0.0;
        }
    }
    // Allocate array for accelerations and initialize them to zero
    double** accelerations = allocate_2d_array(n_atoms, 3);
    for (int i = 0; i < n_atoms; i++) {
        for (int j = 0; j < 3; j++) {
            accelerations[i][j] = 0.0;
        }
    }

    // Open file for writing the trajectory
    const char* trajectory_name = "trajectory.xyz";
    FILE* trajectory_file = open_output(trajectory_name);  

    // Run 1000 steps of MD simulation
    int n_steps = 1000; //number of steps
    double dt = 0.2; //timestep in fs
    
    double kinetic_energy;
    double potential_energy;
    double total_energy;

    for (int i =0; i < n_steps; i++){

        // Update positions, velocities and accelerations
        update_positions(coords, velocities, accelerations, dt, n_atoms);
        update_velocities(velocities, accelerations, dt, n_atoms); // First velocity update with old accelerations
        calculate_accelerations(coords, masses, n_atoms, epsilon, sigma, distances, accelerations);
        update_velocities(velocities, accelerations, dt, n_atoms); // Second velocity update with new accelerations
        
        // Calculate and print energies
        kinetic_energy = calculate_kinetic_energy(velocities, masses, n_atoms);
        potential_energy = calculate_potential_energy(distances, n_atoms, epsilon, sigma);
        total_energy = calculate_total_energy(kinetic_energy, potential_energy);
        fprintf(trajectory_file, "\n%d\n#Step %d: %10.6e %10.6e %10.6e\n", 
                n_atoms, i, kinetic_energy, potential_energy, total_energy);
        // Print coordinates
        for (int i = 0; i < n_atoms; i++) {
            fprintf(trajectory_file, "Ar %8.6e %8.6e %8.6e\n", 
                    coords[i][0], coords[i][1], coords[i][2]); 
        }
    }
    
    // Close output file
    fclose(trajectory_file);  

    // Free all allocated memory
    free_2d_array(coords, n_atoms);
    free_2d_array(distances, n_atoms);
    free_2d_array(velocities, n_atoms);
    free_2d_array(accelerations, n_atoms);
    free(masses);

    return 0;
}
