// Function headers
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

double one_electron_energy(double* data, int32_t n_up, int32_t mo_num);
double two_electron_energy(int32_t* index, double* value, int32_t n_up, int64_t n_integrals);
double hartree_fock_energy(double nuc_repul, double one_el_energy, double two_el_energy);

double get_integral(int i, int j, int k, int l, const int32_t* index, const double* value, int64_t n_integrals);
double MP2_energy_correction(int32_t* index, double* value, double* mo_energy, int32_t n_up, int32_t mo_num, int64_t n_integrals);

#endif