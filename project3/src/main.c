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

    // Open files for writing the output
    const char* trajectory_name = "trajectory.xyz";
    FILE* trajectory_file = open_output(trajectory_name); 
    const char* extended_name = "trajectory_velocity.xyz";
    FILE* extended_file = open_output(extended_name); 
    const char* acceleration_name = "acceleration";
    FILE* acceleration_file = open_output(acceleration_name);

    // Run 1000 steps of MD simulation
    int n_steps = 1000; //number of steps
    double dt = 0.2; //timestep in fs
    
    double kinetic_energy;
    double potential_energy;
    double total_energy;
    double previous_energy;

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
        // Check if the energy is conserved or varies by more than 5 %
        check_energy(previous_energy, total_energy, i);

        // Print output
        print_output(trajectory_file, extended_file, acceleration_file, 
                     n_atoms, i, kinetic_energy, potential_energy, total_energy,
                     coords, velocities, accelerations);       
    }
    
    // Close output files
    fclose(trajectory_file);  
    fclose(extended_file);
    fclose(acceleration_file);

    // Free all allocated memory
    free_2d_array(coords, n_atoms);
    free_2d_array(distances, n_atoms);
    free_2d_array(velocities, n_atoms);
    free_2d_array(accelerations, n_atoms);
    free(masses);

    return 0;
}
