/**
 * @file headers.h
 * @brief Contains the function headers.
 */

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

double** allocate_2d_array(int rows, int cols);
void free_2d_array(double** array, int rows);
int read_natoms(const char* filename);
void read_coords_and_masses(const char* filename, double** coords, double* masses, int n_atoms);
int validate_atoms(double* masses, double* epsilon, double* sigma, int n_atoms);
void initialize_velocities(double** velocities, double* masses, double temperature, int n_atoms);
void calculate_distances(double** coords, double** distances, int n_atoms);
double calculate_potential_energy(double** distances, int n_atoms, double epsilon, double sigma);
double calculate_kinetic_energy(double** velocities, double* masses, int n_atoms);
double calculate_total_energy(double kinetic_energy, double potential_energy);
void thermostat(double kinetic_energy, double temperature, double** velocities, int n_atoms);
void check_energy(double previous_energy, double total_energy, int step);
void calculate_accelerations(double** coords, double* masses, int n_atoms, double epsilon, double sigma, double** distances, double** accelerations);
void update_positions(double** coords, double** velocities, double** accelerations, double dt, int n_atoms);
void update_velocities(double** velocities, double** accelerations, double dt, int n_atoms);
FILE* open_output(const char* filename);
void print_output(FILE* trajectory_file, FILE* energy_file, FILE* extended_file, FILE* acceleration_file, int n_atoms, int step, double kinetic_energy, double potential_energy, double total_energy, double** coords, double** velocities, double** accelerations); 
#endif

// Constants
#ifndef CONSTANTS_H
#define CONSTANTS_H

#define ARGON_MASS 39.948
#define ARGON_EPSILON 0.0661 // j/mol
#define ARGON_SIGMA 0.3345 // nm
#define PI 3.14159265358979323846
#define R 8.31446261815324 // Ideal gas constant in J/K/mol

#endif
