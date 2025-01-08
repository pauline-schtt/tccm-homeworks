#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"

/**
 * @brief Allocates a 2D array of doubles
 * @param rows Number of rows in the array
 * @param cols Number of columns in the array
 * @return Pointer to the allocated 2D array
 * @throws Exits with code 1 if memory allocation fails
 */
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

/**
 * @brief Frees a 2D array from memory
 * @param array Pointer to the 2D array to be freed
 * @param rows Number of rows in the array
 */
void free_2d_array(double** array, int rows) {
    if (array == NULL) return;
    
    for (int i = 0; i < rows; i++) {
        free(array[i]);
    }
    free(array);
}

/**
 * @brief Reads the number of atoms from an input file
 * @param filename Name of the input file
 * @return Number of atoms
 * @throws Exits with code 1 if file cannot be opened
 */
int read_natoms(const char* filename) {
    // Open file
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Could not open file %s, please check whether the filename is correct\n", filename);
        exit(1);
    }
    // Read number of atoms from first line of the file
    int n_atoms;
    fscanf(file, "%d", &n_atoms);
    fclose(file);
    return n_atoms;
}

/**
 * @brief Reads atomic coordinates and masses from an input file
 * @param filename Name of the input file
 * @param coords 2D array to store coordinates
 * @param masses Array to store atomic masses
 * @param n_atoms Number of atoms
 * @throws Exits with code 1 if file cannot be opened
 */
void read_coords_and_masses(const char* filename,
                            double** coords,
                            double* masses,
                            int n_atoms) {
    // Open file
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Could not open file %s\n", filename);
        exit(1);
    }
    
    // Skip the first line containing the number of atoms
    int dummy; //!< Dummy variable used by read_coords_and_masses for reading the input
    fscanf(file, "%d", &dummy);
    
    // Read coordinates and masses
    for (int i = 0; i < n_atoms; i++) {
        fscanf(file, "%lf %lf %lf %lf", &coords[i][0], &coords[i][1], &coords[i][2], &masses[i]);
    }
    
    fclose(file);
}

/**
 * @brief Validates that all atoms in the system are argon
 * @param masses Array of atomic masses
 * @param epsilon Pointer to store the epsilon parameter
 * @param sigma Pointer to store the sigma parameter
 * @param n_atoms Number of atoms
 * @return 1 if all atoms are valid, 0 otherwise
 */
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

/**
 * @brief Opens a file for writing output
 * @param filename Name of the output file
 * @return Pointer to the opened file
 * @throws Exits with code 1 if file cannot be opened
 */
FILE* open_output(const char* filename) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Could not open file %s for writing.\n", filename);
        exit(1);
    }
    return file;
}

/**
 * @brief Prints simulation output
 * @param trajectory_file File pointer for trajectory output
 * @param energy_file File pointer for energy output
 * @param extended_file File pointer for extended information
 * @param acceleration_file File pointer for acceleration data
 * @param n_atoms Number of atoms
 * @param step Current simulation step
 * @param kinetic_energy Current kinetic energy
 * @param potential_energy Current potential energy
 * @param total_energy Current total energy
 * @param coords Array of atomic coordinates
 * @param velocities Array of atomic velocities
 * @param accelerations Array of atomic accelerations
 */
void print_output(FILE* trajectory_file, FILE* energy_file, FILE* extended_file, FILE* acceleration_file, 
                 int n_atoms, int step, double kinetic_energy, double potential_energy, double total_energy,
                 double** coords, double** velocities, double** accelerations) {
    
    // Print comment line with number of atoms, step and energies
    fprintf(trajectory_file, "%d\nStep %d: E(kin) = %10.8f, E(pot) = %10.8f, E(tot) = %10.8f\n",
            n_atoms, step, kinetic_energy, potential_energy, total_energy);
    
    // Print energies to separate file
    fprintf(energy_file, "%10.8f %10.8f %10.8f\n",
            kinetic_energy, potential_energy, total_energy);    

    // Print coordinates, velocities and accelerations
    for (int i = 0; i < n_atoms; i++) {
        fprintf(trajectory_file, "Ar    %10.6f %10.6f %10.6f\n",
                coords[i][0], coords[i][1], coords[i][2]);
        fprintf(extended_file, "Ar     %10.6f %10.6f %10.6f     %10.6f %10.6f %10.6f\n",
                coords[i][0], coords[i][1], coords[i][2], velocities[i][0], velocities[i][1], velocities[i][2]);
        fprintf(acceleration_file, "Ar    %10.6f %10.6f %10.6f\n",
                accelerations[i][0], accelerations[i][1], accelerations[i][2]);
    }
}

/**
 * @brief Calculates distances between all pairs of atoms
 * @param coords Array of atomic coordinates
 * @param distances 2D array to store pairwise distances
 * @param n_atoms Number of atoms
 */
void calculate_distances(double** coords, double** distances, int n_atoms) {
    
    for (int i = 0; i < n_atoms; i++) {
        distances[i][i] = 0.0;  // Distance to the atom itself is 0.0
        for (int j = i + 1; j < n_atoms; j++) {
            double dx = coords[i][0] - coords[j][0]; //!< Variable used by calculate_distances() for the distance in x 
            double dy = coords[i][1] - coords[j][1]; //!< Variable used by calculate_distances() for the distance in y
            double dz = coords[i][2] - coords[j][2]; //!< Variable used by calculate_distances() for the distance in z
            
            double r = sqrt(dx*dx + dy*dy + dz*dz);  //!< Varible used by calculate_distances() for the distance in 3D
            // Pairwise distance matrix is symmetric
            distances[i][j] = r;
            distances[j][i] = r;
        }
    }
}

/**
 * @brief Calculates the Lennard-Jones potential between two atoms
 * @param r Distance between atoms
 * @param epsilon Epsilon parameter for LJ potential
 * @param sigma Sigma parameter for LJ potential
 * @return Value of the Lennard-Jones potential
 */
static double lennard_jones_potential(double r,
                                      double epsilon,
                                      double sigma) {
    double sigma_r = sigma / r; //!< Variable used by lennard_jones_potential() for storing sigma/r
    double sigma_r_6 = pow(sigma_r, 6); //!< Variable used by lennard_jones_potential() for storing the sixth power of sigma_r
    double sigma_r_12 = sigma_r_6 * sigma_r_6; //!< Variable used by lennard_jones_potential() for storing the 12th power of sigma:r
    return 4.0 * epsilon * (sigma_r_12 - sigma_r_6);
}

/**
 * @brief Calculates total potential energy of the system
 * @param distances 2D array of pairwise distances
 * @param n_atoms Number of atoms
 * @param epsilon Epsilon parameter for LJ potential
 * @param sigma Sigma parameter for LJ potential
 * @return Total potential energy of the system
 */
double calculate_potential_energy(double** distances,
                                  int n_atoms,
                                  double epsilon,
                                  double sigma) {
    double total_potential = 0.0; //!< Variable used by calculate_potential_energy() for updating the potential energy
    
    // Sum over all unique pairs of i and j where j > i
    for (int i = 0; i < n_atoms; i++) {
        for (int j = i + 1; j < n_atoms; j++) {
            double r = distances[i][j];
            total_potential += lennard_jones_potential(r, epsilon, sigma);
        }
    }

    return total_potential;
}

/**
 * @brief Calculates total kinetic energy of the system
 * @param velocities Array of atomic velocities
 * @param masses Array of atomic masses
 * @param n_atoms Number of atoms
 * @return Total kinetic energy of the system
 */
double calculate_kinetic_energy(double** velocities,
                                double* masses,
                                int n_atoms) {
    double total_kinetic = 0.0; //!< Variable used by calculate_kinetic_energy() for updating the kinetic energy
    
    for (int i = 0; i < n_atoms; i++) {
        double v_squared = velocities[i][0] * velocities[i][0] + 
                           velocities[i][1] * velocities[i][1] +
                           velocities[i][2] * velocities[i][2]; //!< Varibale used by calculate_kinetic_energy() for the sum of velocity squares
        total_kinetic += 0.5 * masses[i] * v_squared;
    }
    
    return total_kinetic;
}

/**
 * @brief Calculates total energy of the system
 * @param kinetic_energy Current kinetic energy
 * @param potential_energy Current potential energy
 * @return Total energy of the system
 */
double calculate_total_energy(double kinetic_energy, double potential_energy) {
    return kinetic_energy + potential_energy;
}

/**
 * @brief Checks energy conservation between simulation steps
 * @param previous_energy Energy from previous step
 * @param total_energy Current total energy
 * @param step Current simulation step
 * @warning Prints warning if energy varies by more than 10%
 */
void check_energy(double previous_energy, double total_energy, int step) {
    double difference = abs(total_energy - previous_energy); //!< Variable used by check_energy() for the difference in total energy between subsequent steps
    if (difference > 0.10 * abs(previous_energy)) {
        printf("WARNING: The total energy is varying by more than 10 %% in step %5d.\n", step);
    }
}

/**
 * @brief Helper function for acceleration calculation
 * @param r Distance between atoms
 * @param epsilon Epsilon parameter for LJ potential
 * @param sigma Sigma parameter for LJ potential
 * @return U value
 */
static double calculate_U(double r,
                          double epsilon,
                          double sigma) {
    double sigma_r = sigma / r; //!< Variable used by calculate_U() for storing sigma/r
    double sigma_r_6 = pow(sigma_r, 6); //!< Variable used by calculate_U() for the sixth power of sigma_r
    double sigma_r_12 = sigma_r_6 * sigma_r_6; //!< Variable used by calculate_U() for the 12th power of sigma_r
    return 24.0 * (epsilon / r) * (sigma_r_6 - 2.0 * sigma_r_12); 
}

/**
 * @brief Calculates acceleration vectors for all atoms
 * @param coords Array of atomic coordinates
 * @param masses Array of atomic masses
 * @param n_atoms Number of atoms
 * @param epsilon Epsilon parameter for LJ potential
 * @param sigma Sigma parameter for LJ potential
 * @param distances 2D array of pairwise distances
 * @param accelerations Array to store calculated accelerations
 */
void calculate_accelerations(double** coords,
                             double* masses,
                             int n_atoms, 
                             double epsilon,
                             double sigma,
                             double** distances, 
                             double** accelerations) {

    for (int i = 0; i < n_atoms; i++) { // Initialize accelerations to zero
        accelerations[i][0] = 0.0;
        accelerations[i][1] = 0.0;
        accelerations[i][2] = 0.0;
    }

    calculate_distances(coords, distances, n_atoms); // Update the 2D array of distances    

    for (int i = 0; i < n_atoms; i++) {
        for (int j = 0; j < n_atoms; j++) {
            if (i != j) {  // To avoid self-interaction
                double r = distances[i][j]; //!< Variable used by calculate_accelerations() for the distance
                double U = calculate_U(r, epsilon, sigma); //!< Variable used by calculate_accelerations() for the scaled potential U
                
                double dx = coords[i][0] - coords[j][0]; //!< Variable used by calculate_accelerations() for the distance in x
                double dy = coords[i][1] - coords[j][1]; //!< Variable used by calculate_accelerations() for the distance in y
                double dz = coords[i][2] - coords[j][2]; //!< Variable used by calculate_accelerations() for the distance in z
                
                // Calculate accelerations
                accelerations[i][0] += - (1 / masses[i]) * U * (dx / r);
                accelerations[i][1] += - (1 / masses[i]) * U * (dy / r);
                accelerations[i][2] += - (1 / masses[i]) * U * (dz / r);
            }
        }
    }
}

/**
 * @brief Updates atomic positions using Verlet algorithm
 * @param coords Array of atomic coordinates
 * @param velocities Array of atomic velocities
 * @param accelerations Array of atomic accelerations
 * @param dt Time step
 * @param n_atoms Number of atoms
 */
void update_positions(double** coords,
                      double** velocities,
                      double** accelerations,
                      double dt,
                      int n_atoms) {
    double dt_2 = 0.5 * dt * dt; //!< Variable used by update_positions for the square of dt
    for (int i = 0; i < n_atoms; i++) {
        // Calculate new positions
        coords[i][0] += velocities[i][0] * dt + accelerations[i][0] * dt_2;
        coords[i][1] += velocities[i][1] * dt + accelerations[i][1] * dt_2;
        coords[i][2] += velocities[i][2] * dt + accelerations[i][2] * dt_2;
    }
}

/**
 * @brief Updates atomic velocities using Verlet algorithm
 * @param velocities Array of atomic velocities
 * @param accelerations Array of atomic accelerations
 * @param dt Time step
 * @param n_atoms Number of atoms
 */
void update_velocities(double** velocities,
                       double** accelerations,
                       double dt,
                       int n_atoms) {
    for (int i = 0; i < n_atoms; i++) {
        // Calculate new velocities
        velocities[i][0] += 0.5 * accelerations[i][0] * dt;
        velocities[i][1] += 0.5 * accelerations[i][1] * dt;
        velocities[i][2] += 0.5 * accelerations[i][2] * dt;
    }
}
