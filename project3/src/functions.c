#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"

// Function to allocate a 2D array
double** allocate_2d_array(int rows, int cols) {
    double** array = (double**)malloc(rows * sizeof(double*));
    if (array == NULL) {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(1);
    }
    
    for (int i = 0; i < rows; i++) {
        array[i] = (double*)malloc(cols * sizeof(double));
        if (array[i] == NULL) {
            fprintf(stderr, "Memory allocation failed!\n");
            for (int j = 0; j < i; j++) { // Free previously allocated memory
                free(array[j]);
            }
            free(array);
            exit(1);
        }
    }
    return array;
}

// Function to free a 2D array
void free_2d_array(double** array, int rows) {
    if (array == NULL) return;
    
    for (int i = 0; i < rows; i++) {
        free(array[i]);
    }
    free(array);
}

// Function to read number of atoms from input file
int read_natoms(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Could not open file %s, please check whether the filename is correct\n", filename);
        exit(1);
    }
    
    int n_atoms;
    fscanf(file, "%d", &n_atoms);
    fclose(file);
    return n_atoms;
}

// Function to read coordinates and masses
void read_coords_and_masses(const char* filename,
                            double** coords,
                            double* masses,
                            int n_atoms) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Could not open file %s\n", filename);
        exit(1);
    }
    
    // Skip the first line containing the number of atoms
    int dummy;
    fscanf(file, "%d", &dummy);
    
    // Read coordinates and masses
    for (int i = 0; i < n_atoms; i++) {
        fscanf(file, "%lf %lf %lf %lf", &coords[i][0], &coords[i][1], &coords[i][2], &masses[i]);
    }
    
    fclose(file);
}

// Function to validate atoms in the input file
int validate_atoms(double* masses, double* epsilon, double* sigma, int n_atoms) {
    for (int i = 0; i < n_atoms; i++) {
        if (masses[i] != ARGON_MASS) {
            fprintf(stderr, "Error: This MD engine isn't parametrized for some of the atoms in the input file\n");
            return 0;
        }
    }

    *epsilon = ARGON_EPSILON;
    *sigma = ARGON_SIGMA;
    return 1;
}

// Function to calculate internuclear distances
double** calculate_distances(double** coords, int n_atoms) {
    double** distances = allocate_2d_array(n_atoms, n_atoms);
    
    for (int i = 0; i < n_atoms; i++) {
        distances[i][i] = 0.0;  // Distance to the atom itself is 0.0
        for (int j = i + 1; j < n_atoms; j++) {
            double dx = coords[i][0] - coords[j][0];
            double dy = coords[i][1] - coords[j][1];
            double dz = coords[i][2] - coords[j][2];
            
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            // Pairwise distance matrix is symmetric
            distances[i][j] = r;
            distances[j][i] = r;
        }
    }
    
    return distances; // 2D array of size n_atoms x n_atoms
}

// Function to calculate the Lennard-Jones potential
static double lennard_jones_potential(double r,
                                      double epsilon,
                                      double sigma) {
    double sigma_r = sigma / r;
    double sigma_r_6 = pow(sigma_r, 6);
    double sigma_r_12 = sigma_r_6 * sigma_r_6;
    return 4.0 * epsilon * (sigma_r_12 - sigma_r_6);
}

// Function to calculate total potential energy
double calculate_potential_energy(double** distances,
                                  int n_atoms,
                                  double epsilon,
                                  double sigma) {
    double total_potential = 0.0;
    
    // Sum over all unique pairs of i and j where j > i
    for (int i = 0; i < n_atoms; i++) {
        for (int j = i + 1; j < n_atoms; j++) {
            double r = distances[i][j];
            total_potential += lennard_jones_potential(r, epsilon, sigma);
        }
    }

    return total_potential;
}

// Function to calculate kinetic energy
double calculate_kinetic_energy(double** velocities,
                                double* masses,
                                int n_atoms) {
    double total_kinetic = 0.0;
    
    for (int i = 0; i < n_atoms; i++) {
        double v_squared = velocities[i][0] * velocities[i][0] + 
                           velocities[i][1] * velocities[i][1] +
                           velocities[i][2] * velocities[i][2];
        total_kinetic += 0.5 * masses[i] * v_squared;
    }
    
    return total_kinetic;
}

// Function to calculate total energy
double calculate_total_energy(double kinetic_energy, double potential_energy) {
    return kinetic_energy + potential_energy;
}

// Helper function for acceleration calculation
static double calculate_U(double r,
                          double epsilon,
                          double sigma) {
    double sigma_r = sigma / r;
    double sigma_r_6 = pow(sigma_r, 6);
    double sigma_r_12 = sigma_r_6 * sigma_r_6;
    return 24.0 * (epsilon / r) * (sigma_r_6 - 2.0 * sigma_r_12);
}

// Function to calculate acceleration vectors
void calculate_accelerations(double** coords,
                             double* masses,
                             int n_atoms, 
                             double epsilon,
                             double sigma,
                             double** distances, 
                             double** accelerations) {
    for (int i = 0; i < n_atoms; i++) {
        for (int j = 0; j < n_atoms; j++) {
            if (i != j) {  // To avoid self-interaction
                double r = distances[i][j];
                double U = calculate_U(r, epsilon, sigma);
                
                double dx = coords[i][0] - coords[j][0];
                double dy = coords[i][1] - coords[j][1];
                double dz = coords[i][2] - coords[j][2];
                
                // Calculate accelerations
                accelerations[i][0] += - (1 / masses[i]) * U * (dx / r);
                accelerations[i][1] += - (1 / masses[i]) * U * (dy / r);
                accelerations[i][2] += - (1 / masses[i]) * U * (dz / r);
            }
        }
    }
}
