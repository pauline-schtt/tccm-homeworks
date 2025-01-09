/**
 * @file functions.c
 * @brief Contains the functions associated with the HF and MP2 energy calculation.
 */

#include <stdio.h>

#include <stdint.h>
#include <stdio.h>

/**
 * @brief Calculates the one-electron energy contribution to the Hartree-Fock energy
 * @param data Array containing one-electron integrals
 * @param n_up Number of occupied orbitals
 * @param mo_num Total number of molecular orbitals
 * @return One-electron energy contribution
 */
double one_electron_energy(double* data, int32_t n_up, int32_t mo_num){
    // printf("\nOne-electron integrals:\n i    value\n");
    /**
    * @var double one_el_energy
    * @brief Variable used by one_electron_energy() while calculating the one-electron energy
    */
    double one_el_energy = 0; 
    for (int i=0; i < n_up; i++) { // Iterate over the occupied orbitals
        one_el_energy += 2 * data[i*mo_num+i]; // Get the value of the integral and add
        // printf("%2d    %9.6lf\n", i+1, data[i*mo_num+i]);
    }
    return one_el_energy;
}

/**
 * @brief Calculates the two-electron energy contribution to the Hartree-Fock energy
 * @param index Array containing four-index combinations for two-electron integrals
 * @param value Array containing values of two-electron integrals
 * @param n_up Number of occupied orbitals
 * @param n_integrals Total number of two-electron integrals
 * @return Two-electron energy contribution
 */
double two_electron_energy(int32_t* index, double* value, int32_t n_up, int64_t n_integrals) {
    // printf("\nTwo-electron integrals:\n       i  j  k  l     value\n");
    /**
     * @var double two_el_energy
     * @brief Variable used by two_electron_energy() while calculating the two-electron interaction energy
     */
    double two_el_energy = 0;
    for (int n=0; n < n_integrals; n++) { //! Iterate over the stored integrals
        // Get the indices
        /**
         * @var int i
         * @brief Variable used by two_electron_energy() for storing the index i
         */
        int i = index[4*n]; 
        /**
         * @var int j
         * @brief Variable used by two_electron_energy() for storing the index j
         */
        int j = index[4*n+1];
        /**
         * @var int k
         * @brief Variable used by two_electron_energy() for storing the index k
         */
        int k = index[4*n+2]; 
        /**
         * @var int l
         * @brief Variable used by two_electron_energy() for storing the index l
         */
        int l = index[4*n+3];
        if (i < n_up && j < n_up) { // Check if the first two indices belong to occupied orbitals
            if (i == j && j == k && k == l) {
                two_el_energy += value[n]; // if all indices are the same, only add once
                // printf("2J-K: %2d %2d %2d %2d    %9.6lf\n", i+1, j+1, k+1, l+1, value[n]);
            }
            else if (k == i && l ==j) {
                two_el_energy += (2 * 2 * value[n]); // Add x2 the Coulomb integral, *2 for permutational symmetry
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

/**
 * @brief Calculates the total Hartree-Fock energy
 * @param nuc_repul Nuclear repulsion energy
 * @param one_el_energy One-electron energy contribution
 * @param two_el_energy Two-electron energy contribution
 * @return Total Hartree-Fock energy
 */
double hartree_fock_energy(double nuc_repul, double one_el_energy, double two_el_energy) {
    return nuc_repul + one_el_energy + two_el_energy;
}

/**
 * @brief Retrieves integral value
 * @param i First index
 * @param j Second index
 * @param k Third index
 * @param l Fourth index
 * @param index Array containing four-index combinations
 * @param value Array containing integral values
 * @param n_integrals Total number of integrals
 * @return Value of the requested integral or 0.0 if not found
 */
double get_integral(int i, int j, int k, int l, const int32_t* index, const double* value, int64_t n_integrals) {
    for (int64_t n = 0; n < n_integrals; n++) { // Try both possible permutations       
        if ((index[4*n] == i && index[4*n+1] == j && index[4*n+2] == k && index[4*n+3] == l) ||
            (index[4*n] == j && index[4*n+1] == i && index[4*n+2] == l && index[4*n+3] == k)) {
            return value[n];
        }
    }
    return 0.0;
}

/**
 * @brief MP2 energy correction calculation
 * @param index Array containing four-index combinations
 * @param value Array containing integral values
 * @param mo_energy Array of molecular orbital energies
 * @param n_up Number of occupied orbitals
 * @param n_integrals Total number of integrals
 * @return MP2 energy correction
 */
double MP2_energy_correction(int32_t* index, double* value, double* mo_energy, int32_t n_up, int64_t n_integrals) { 
    /**
     * @var double MP2_energy
     * @brief Variable used by MP2_energy_correction() while calculating the MP2 energy
     */
    double MP2_energy = 0; 
    /**
     * @var double ijab
     * @brief Variable used by MP2_alter() to store the integral with indices ijab
     */
    double ijab = 0; 
    /**
     * @var double ijba
     * @brief Variable used by MP2_alter() to store the integral with indices ijba
     */
    double ijba = 0;
    /**
     * @var double denominator
     * @brief Variable used by MP2_alter() to store the denominator
     */
    double denominator = 0; 
    /**
     * @var double symmetry
     * @brief Variable used by MP2_alter to account for permutational symmetry of the integrals
     */
    double symmetry = 0; 
    for (int n=0; n < n_integrals; n++) { //Iterate over the stored integrals
        // Get the indices
        /**
         * @var int i
         * @brief Variable used by MP2_alter() for storing index i
         */
        int i = index[4*n]; 
        /**
         * @var int j
         * @brief Variable used by MP2_alter() for storing index j
         */
        int j = index[4*n+1]; 
        /**
         * @var int a
         * @brief Variable used by MP2_alter() for storing index a
         */
        int a = index[4*n+2];
        /**
         * @var int b
         * @brief Variable used by MP2_alter() for storing index b
         */
        int b = index[4*n+3];
        // Use <ij|ab> = <ab|ij> 
        // Check if the first two indices belong to virtual orbitals and the last two to occupied
        if (i >= n_up && j >= n_up && a < n_up && b < n_up) { 
            ijab = value[n];
            //printf("ijab: %2d %2d %2d %2d    %9.6lf\n", i+1, j+1, a+1, b+1, ijab);
            denominator = mo_energy[a] + mo_energy[b] - mo_energy[i] - mo_energy[j];
            ijba = get_integral(i, j, b, a, index, value, n_integrals);      
            // Account for permutational symmetry
	    if (i == j && a == b) {
                symmetry = 1.0;
            }
            else {
                symmetry = 2.0;
            }
            // Add integrals to MP2 correction
            MP2_energy += symmetry * ((ijab * (2.0 * ijab - ijba))/denominator);
        }  
    }
    return MP2_energy;
}
