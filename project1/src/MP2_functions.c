// Functions for calculating the MP2 energy correction#

#include <stdint.h>

// Function to get integral values considering 8-fold symmetry
double get_integral(int i, int j, int k, int l, const int32_t* index, const double* value, int64_t n_integrals) {
    for (int64_t n = 0; n < n_integrals; n++) { // Try all possible permutations
        if ((index[4*n] == i && index[4*n+1] == j && index[4*n+2] == k && index[4*n+3] == l) ||
            (index[4*n] == i && index[4*n+1] == l && index[4*n+2] == k && index[4*n+3] == j) ||
            (index[4*n] == k && index[4*n+1] == l && index[4*n+2] == i && index[4*n+3] == j) ||
            (index[4*n] == k && index[4*n+1] == j && index[4*n+2] == i && index[4*n+3] == l) ||
            (index[4*n] == j && index[4*n+1] == i && index[4*n+2] == l && index[4*n+3] == k) ||
            (index[4*n] == l && index[4*n+1] == i && index[4*n+2] == j && index[4*n+3] == k) ||
            (index[4*n] == l && index[4*n+1] == k && index[4*n+2] == j && index[4*n+3] == i) ||
            (index[4*n] == j && index[4*n+1] == k && index[4*n+2] == l && index[4*n+3] == i)) {
            return value[n];
        }
    }
    return 0.0;  // Return 0 if integral not found
}

// Function to calculate MP2 energy
double MP2_energy_correction(int32_t* index, double* value, double* mo_energy, int32_t n_up, int32_t mo_num, int64_t n_integrals) {
    double MP2_energy = 0.0;
    for (int i = 0; i < n_up; i++) {
        for (int j = 0; j < n_up; j++) {
            for (int a = n_up; a < mo_num; a++) {
                for (int b = n_up; b < mo_num; b++) {
                    double ijab = get_integral(i, j, a, b, index, value, n_integrals);
                    double ijba = get_integral(i, j, b, a, index, value, n_integrals);
                    
                    MP2_energy += ijab * (2.0 * ijab - ijba) / (mo_energy[i] + mo_energy[j] - mo_energy[a] - mo_energy[b]);
                }
            }
        }
    }
    return MP2_energy;
}
