// Imports
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <trexio.h>

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

int main(int argc, char *argv[]) {
    char* filename; //! Variable for the name of the HDF5 file

    // Check if a HDF5 file was specified as argument
    if (argc != 2) {
        fprintf(stderr, "No HDF5 file containing the data was specified. You can do so by using the program as follows: ./HF 'path/to/hdf5' or by providing the path for your file below:\nPath to HDF5 file: ");
        scanf("%s", filename);
    }
    else if (argc == 2) {
        filename = argv[1];
    }

    // Greet the user
    // Start with an ASCII art of the program name
    printf(" ___  ___  ________      _____ ______   ________    _______     \n");
    printf("|\\  \\|\\  \\|\\  _____\\    |\\   _ \\  _   \\|\\   __  \\  /  ___  \\    \n");
    printf("\\ \\  \\\\\\  \\ \\  \\__/     \\ \\  \\\\\\__\\ \\  \\ \\  \\|\\  \\/__/|_/  /|   \n");
    printf(" \\ \\   __  \\ \\   __\\     \\ \\  \\\\|__| \\  \\ \\   ____\\__|//  / /   \n");
    printf("  \\ \\  \\ \\  \\ \\  \\_|      \\ \\  \\    \\ \\  \\ \\  \\___|   /  /_/__  \n");
    printf("   \\ \\__\\ \\__\\ \\__\\        \\ \\__\\    \\ \\__\\ \\__\\     |\\________\\\n");
    printf("    \\|__|\\|__|\\|__|         \\|__|     \\|__|\\|__|      \\|_______|\n");
    printf("\nWelcome to the Hartree-Fock and MP2 energy calculation program.\n");

    // Open TREXIO file for reading the data
    trexio_exit_code rc;
    trexio_t* trexio_file = trexio_open(filename, 'r', TREXIO_AUTO, &rc);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error: %s\n", trexio_string_of_error(rc));
        exit(1);
    }

    // Read the nuclear repulsion energy
    double nuc_repul; //! Variable where the nuclear repulsion energy is read
    rc = trexio_read_nucleus_repulsion(trexio_file, &nuc_repul);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading nuclear repulsion energy:\n%s\n",
        trexio_string_of_error(rc));
        exit(1);
    }

    // Obtain the number of occupied orbitals
    int32_t n_up; //! Variable where the number of spin-up electrons is read
    rc = trexio_read_electron_up_num(trexio_file, &n_up);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading number of spin-up electrons: n%s\n",
        trexio_string_of_error(rc));
        exit(1);
    }

    // Obtain the number of molecular orbitals
    int32_t mo_num; //! Variable where the number of molecular orbitals is read
    rc = trexio_read_mo_num(trexio_file, &mo_num);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading number of molecular orbitals: n%s\n",
        trexio_string_of_error(rc));
        exit(1);
    }

    // Read in the one-electron integrals
    int64_t elements = mo_num*mo_num; //! Variable for the number of elements to store in data array
    // Allocate data array
    double* data = malloc(elements*sizeof(double)); //! Request mo_num x mo_num doubles for storing the one-electron integrals
    if (data == NULL) { // Check that the allocation was OK
        fprintf(stderr, "Allocation of data array for one-electron integrals failed\n");
        exit(-1);
    }
    // Use TREXIO function to fill data array with one-electron integrals
    rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, data);
    // Check the return code to be sure reading was OK
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading one-electron integrals: n%s\n",
        trexio_string_of_error(rc));
        exit(1);
    }

    // Read in the two-electron integrals

    // Get number of non-zero integrals
    int64_t n_integrals; //! Variable for number of non-zero two-electron integrals
    rc = trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);
    // Check the return code to be sure reading was OK
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading number of non-zero two-electron integrals: n%s\n",
        trexio_string_of_error(rc));
        exit(1);
    }

    // Allocate memory for storing the indices
    int32_t* index = malloc(4 * n_integrals * sizeof(int32_t)); 
    if (index == NULL) { // Check that the allocation was OK
        fprintf(stderr, "Allocation of index array for two-electron integrals failed\n");
        exit(-1);
    }

    // Allocate memory for storing the values of the integrals
    double* value = malloc(n_integrals * sizeof(double)); //! Array of values of the two-electron integrals
    if (value == NULL) { // Check that the allocation was OK 
        fprintf(stderr, "Allocation of value array for two-electron integrals failed\n");
        exit(-1);
    }

    // Read in the integrals
    int64_t buffer_size = n_integrals;
    rc = trexio_read_mo_2e_int_eri(trexio_file, 0, &buffer_size, index, value);
    // Check the return code to be sure reading was OK
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading two-electron integrals: n%s\n",
        trexio_string_of_error(rc));
        exit(1);
    }

    // Check if buffer size is equal to number of integrals
    if (buffer_size != n_integrals) {
        fprintf(stderr, "Not all two-electron integrals were read correctly\n");
        exit(1);
    }

    //Calculate one-electron energy contribution
    printf("\nCalculating the Hartree-Fock energy...\n");
    // printf("\nOne-electron integrals:\n i    value\n");
    double one_el_energy = 0;
    for (int i=0; i < n_up; i++) { //! Iterate over the occupied orbitals
        one_el_energy = one_el_energy + 2 * data[i*mo_num+i]; //! Get the value of the integral and add
        // printf("%2d    %9.6lf\n", i+1, data[i*mo_num+i]);
    }

    // Calculate the two-electron energy contribution
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
                two_el_energy = two_el_energy + value[n]; //! if all indices are the same, only add once
                // printf("2J-K: %2d %2d %2d %2d    %9.6lf\n", i+1, j+1, k+1, l+1, value[n]);
            }
            else if (k == i && l ==j) {
                two_el_energy = two_el_energy + (2 * 2 * value[n]); //! Add x2 the Coulomb integral, *2 for permutational symmetry
                // printf("J:    %2d %2d %2d %2d    %9.6lf\n", i+1, j+1, k+1, l+1, value[n]);
            }
            else if (i == j && k == l) {
                two_el_energy = two_el_energy - (2 * value[n]); //! Substract the exchange integral, *2 for permutational symmetry
                // printf("K:    %2d %2d %2d %2d    %9.6lf\n", i+1, j+1, k+1, l+1, value[n]);
            }
        }
        else if (l >= n_up) {
            break; //! Break the loop once there are only integrals with virtual orbitals in the list
        }
    }

    // Calculate the Hartree-Fock energy
    double HF_energy = nuc_repul + one_el_energy + two_el_energy;
    printf("Done!\n");

    // Calculate MP2 energy correction
    double* mo_energy = malloc(mo_num * sizeof(double)); // First, read the orbital energies
    if (mo_energy == NULL) {
        fprintf(stderr, "Failed to allocate memory for orbital energies\n");
        exit(1);
    }
    rc = trexio_read_mo_energy(trexio_file, mo_energy);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading orbital energies: %s\n", 
                trexio_string_of_error(rc));
        free(mo_energy);
        exit(1);
    }

    // Calculate MP2 energy
    double MP2_energy = 0.0;
    printf("\nCalculating the MP2 energy correction...\n");
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
    printf("Done!\n"); 

    // Print a summary
    printf("\n################## Energy Summary ##################\n");
    printf("\nNuclear repulsion energy:             %9.6lf\n", nuc_repul);
    printf("One-electron energy:                  %9.6lf\n", one_el_energy);
    printf("Two-electron energy:                  %9.6lf\n", two_el_energy);
    printf("Hartree-Fock energy:                  %9.6lf\n", HF_energy);
    printf("MP2 energy correction:                %9.6lf\n", MP2_energy);
    printf("Total energy (HF + MP2):              %9.6lf\n", HF_energy + MP2_energy);
    printf("\n################ System Information ################\n");
    printf("\nNumber of occupied orbitals:          %d\n", n_up);
    printf("Number of molecular orbitals:         %d\n", mo_num);
    printf("Number of two-electron integrals:     %ld\n", n_integrals);

    // Finalize the calculation
    printf("\nYour calculation is done.\n");
    printf("Thank you for using the program!\n");

    // Free the allocated arrays
    free(mo_energy);
    mo_energy = NULL;
    free(value);
    value = NULL;
    free(index); // Reset pointer? Either skip this add "const" in variable defintion
    index = NULL;
    free(data);
    data = NULL;

    // Close the TREXIO file
    rc = trexio_close(trexio_file);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error: %s\n", trexio_string_of_error(rc));
        exit(1);
    }
    trexio_file = NULL;

    /*Print some indexes to get used to the order of the integrals
    for (int n=0; n < 200; n++) {
        printf("%2d %2d %2d %2d\n", index[4*n], index[4*n+1], index[4*n+2], index[4*n+3]);
    }*/

    return 0;
}