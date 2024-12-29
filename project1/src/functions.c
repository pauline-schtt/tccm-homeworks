#include <stdint.h>

//
// Functions associated with calculating the Hartree-Fock energy
//

// Function to calculate one-electron energy contribution

double one_electron_energy(double* data, int32_t n_up, int32_t mo_num){
    // printf("\nOne-electron integrals:\n i    value\n");
    double one_el_energy = 0;
    for (int i=0; i < n_up; i++) { //! Iterate over the occupied orbitals
        one_el_energy += 2 * data[i*mo_num+i]; //! Get the value of the integral and add
        // printf("%2d    %9.6lf\n", i+1, data[i*mo_num+i]);
    }
    return one_el_energy;
}


// Function to calculate the two-electron energy contribution

double two_electron_energy(int32_t* index, double* value, int32_t n_up, int64_t n_integrals) {
    // printf("\nTwo-electron integrals:\n       i  j  k  l     value\n");
    double two_el_energy = 0; //! Variable to store the two-electron interaction energy
    for (int n=0; n < n_integrals; n++) { //! Iterate over the stored integrals
        //Get the indices
        int i = index[4*n];
        int j = index[4*n+1];
        int k = index[4*n+2];
        int l = index[4*n+3];
        if (i < n_up && j < n_up) { //! Check if the first two indices belong to occupied orbitals
            if (i == j && j == k && k == l) {
                two_el_energy += value[n]; //! if all indices are the same, only add once
                // printf("2J-K: %2d %2d %2d %2d    %9.6lf\n", i+1, j+1, k+1, l+1, value[n]);
            }
            else if (k == i && l ==j) {
                two_el_energy += (2 * 2 * value[n]); //! Add x2 the Coulomb integral, *2 for permutational symmetry
                // printf("J:    %2d %2d %2d %2d    %9.6lf\n", i+1, j+1, k+1, l+1, value[n]);
            }
            else if (i == j && k == l) {
                two_el_energy -= 2 * value[n]; //! Substract the exchange integral, *2 for permutational symmetry
                // printf("K:    %2d %2d %2d %2d    %9.6lf\n", i+1, j+1, k+1, l+1, value[n]);
            }
        }
        else if (l >= n_up) {
            break; //! Break the loop once there are only integrals with virtual orbitals in the list
        }
    }
    return two_el_energy;
}

// Function to calculate Hartree-Fock energy

double hartree_fock_energy(double nuc_repul, double one_el_energy, double two_el_energy) {
    return nuc_repul + one_el_energy + two_el_energy;
}

//
// Functions for calculating the MP2 energy correction
//

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
