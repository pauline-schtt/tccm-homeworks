// Function prototypes
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

double** allocate_2d_array(int rows, int cols);
void free_2d_array(double** array, int rows);
int read_natoms(const char* filename);
void read_coords_and_masses(const char* filename, double** coords, double* masses, int n_atoms);
int validate_atoms(double* masses, double* epsilon, double* sigma, int n_atoms);
void calculate_distances(double** coords, double** distances, int n_atoms);
double calculate_potential_energy(double** distances, int n_atoms, double epsilon, double sigma);
double calculate_kinetic_energy(double** velocities, double* masses, int n_atoms);
double calculate_total_energy(double kinetic_energy, double potential_energy);
void calculate_accelerations(double** coords, double* masses, int n_atoms, double epsilon, double sigma, double** distances, double** accelerations);
void update_positions(double** coords, double** velocities, double** accelerations, double dt, int n_atoms);
void update_velocities(double** velocities, double** accelerations, double dt, int n_atoms);
FILE* open_output(const char* filename);

#endif

// Constants
#ifndef CONSTANTS_H
#define CONSTANTS_H

#define ARGON_MASS 39.948
#define ARGON_EPSILON 0.0661 // j/mol
#define ARGON_SIGMA 0.3345 // nm

#endif
