/**
 * @file main.c
 * @brief Contains the main program.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <trexio.h>
#include <time.h>
#include "headers.h" // Function headers

int main(int argc, char *argv[]) {
    // Start timing the entire program
    clock_t start_total = clock();
    clock_t start_hf, end_hf, start_mp2, end_mp2;

    /**
     * @var char* filename
     * @brief Name of the HDF5 file
     * @details This variable stores the path to the HDF5 file provided as input.
     */
    char* filename;

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
    /**
     * @var trexio_exit_code rc
     * @brief TREXIO output code
     * @details This variable is used to store the output code returned by TREXIO operations.
     */
    trexio_exit_code rc;

    /**
     * @var trexio_t* trexio_file
     * @brief TREXIO file handle
     * @details This variable represents the handle for the TREXIO file.
     */
    trexio_t* trexio_file = trexio_open(filename, 'r', TREXIO_AUTO, &rc);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error: %s\n", trexio_string_of_error(rc));
        exit(1);
    }

    // Read the nuclear repulsion energy
    /**
     * @var double nuc_repul
     * @brief Nuclear repulsion energy
     * @details This variable stores the nuclear repulsion energy read from the TREXIO file.
     */
    double nuc_repul;
    rc = trexio_read_nucleus_repulsion(trexio_file, &nuc_repul);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading nuclear repulsion energy:\n%s\n",
                trexio_string_of_error(rc));
        exit(1);
    }

    // Obtain the number of occupied orbitals
    /**
     * @var int32_t n_up
     * @brief Number of spin-up electrons
     * @details This variable stores the number of spin-up electrons read from the TREXIO file.
     */
    int32_t n_up;
    rc = trexio_read_electron_up_num(trexio_file, &n_up);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading number of spin-up electrons: n%s\n",
                trexio_string_of_error(rc));
        exit(1);
    }

    // Obtain the number of molecular orbitals
    /**
     * @var int32_t mo_num
     * @brief Number of molecular orbitals
     * @details This variable stores the number of molecular orbitals read from the TREXIO file.
     */
    int32_t mo_num;
    rc = trexio_read_mo_num(trexio_file, &mo_num);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading number of molecular orbitals: n%s\n",
                trexio_string_of_error(rc));
        exit(1);
    }

    // Allocate array for orbital energies
    /**
     * @var double* mo_energy
     * @brief Array of orbital energies
     * @details This array stores the orbital energies read from the TREXIO file.
     */
    double* mo_energy = malloc(mo_num * sizeof(double));
    if (mo_energy == NULL) {
        fprintf(stderr, "Failed to allocate memory for orbital energies\n");
        exit(1);
    }
    // Read in the orbital energies
    rc = trexio_read_mo_energy(trexio_file, mo_energy);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading orbital energies: %s\n", 
                trexio_string_of_error(rc));
        free(mo_energy);
        exit(1);
    }

    // Read in the one-electron integrals

    /**
     * @var int64_t elements
     * @brief Number of elements in the data array
     * @details This variable represents the total number of elements to store in the data array.
     */
    int64_t elements = mo_num * mo_num;

    /**
     * @var double* data
     * @brief Array for one-electron integrals
     * @details This array stores the one-electron integrals calculated in the Hartree-Fock method.
     */
    double* data = malloc(elements * sizeof(double));
    if (data == NULL) {
        fprintf(stderr, "Allocation of data array for one-electron integrals failed\n");
        free(mo_energy);
        exit(-1);
    }
    // Use TREXIO function to fill data array with one-electron integrals
    rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, data);
    // Check the return code to be sure reading was OK
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading one-electron integrals: n%s\n",
                trexio_string_of_error(rc));
        free(mo_energy);
        free(data);
        exit(1);
    }

    // Read in the two-electron integrals

    // Get number of non-zero integrals
    /**
     * @var int64_t n_integrals
     * @brief Number of non-zero two-electron integrals
     * @details This variable stores the total count of non-zero two-electron integrals read from the TREXIO file.
     */
    int64_t n_integrals;
    rc = trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);
    // Check the return code to be sure reading was OK
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading number of non-zero two-electron integrals: n%s\n",
                trexio_string_of_error(rc));
        free(mo_energy);
        free(data);
        exit(1);
    }

    // Allocate memory for storing the indices
    /**
     * @var int32_t* index
     * @brief Array for two-electron integral indices
     * @details This array stores the indices of the two-electron integrals.
     */
    int32_t* index = malloc(4 * n_integrals * sizeof(int32_t));
    if (index == NULL) {
        fprintf(stderr, "Allocation of index array for two-electron integrals failed\n");
        free(mo_energy);
        free(data);
        exit(-1);
    }

    // Allocate memory for storing the values of the integrals
    /**
     * @var double* value
     * @brief Array for two-electron integral values
     * @details This array stores the values of the two-electron integrals.
     */
    double* value = malloc(n_integrals * sizeof(double));
    if (value == NULL) {
        fprintf(stderr, "Allocation of value array for two-electron integrals failed\n");
        free(mo_energy);
        free(data);
        free(index);
        exit(-1);
    }

    // Read in the integrals
    /**
     * @var int64_t buffer_size
     * @brief Buffer size for two-electron integrals
     * @details This variable is used as a buffer size for reading the two-electron integrals, initially set to n_integrals.
     */
    int64_t buffer_size = n_integrals;
    rc = trexio_read_mo_2e_int_eri(trexio_file, 0, &buffer_size, index, value);
    // Check the return code to be sure reading was OK
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading two-electron integrals: n%s\n",
                trexio_string_of_error(rc));
        free(mo_energy);
        free(data);
        free(index);
        free(value);
        exit(1);
    }

    // Check if buffer size is equal to number of integrals
    if (buffer_size != n_integrals) {
        fprintf(stderr, "Not all two-electron integrals were read correctly\n");
        free(mo_energy);
        free(data);
        free(index);
        free(value);
        exit(1);
    }

    // Close the TREXIO file
    rc = trexio_close(trexio_file);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error: %s\n", trexio_string_of_error(rc));
        exit(1);
    }
    trexio_file = NULL;

    printf("\nCalculating the Hartree-Fock energy...\n");
    
    start_hf = clock(); // Start timing the Hartree-Fock energy calculation

    //Calculate one-electron energy contribution
    /**
     * @var double one_el_energy
     * @brief One-electron energy contribution
     * @details This variable stores the one-electron energy contribution to the Hartree-Fock energy.
     */
    double one_el_energy = one_electron_energy(data, n_up, mo_num);

    // Calculate the two-electron energy contribution
    /**
     * @var double two_el_energy
     * @brief Two-electron energy contribution
     * @details This variable stores the two-electron energy contribution to the Hartree-Fock energy.
     */
    double two_el_energy = two_electron_energy(index, value, n_up, n_integrals);
    
    // Calculate the Hartree-Fock energy
    /**
     * @var double HF_energy
     * @brief Overall Hartree-Fock energy
     * @details This variable stores the overall Hartree-Fock energy, including nuclear, one-electron, and two-electron contributions.
     */
    double HF_energy = hartree_fock_energy(nuc_repul, one_el_energy, two_el_energy);

    end_hf = clock(); // End timing the Hartree-Fock energy calculation
    
    printf("Done!\n");

    printf("\nCalculating the MP2 energy correction...\n");
    
    start_mp2 = clock(); // Start timing the MP2 energy correction calculation

    // Calculate MP2 energy
    /**
     * @var double MP2_energy
     * @brief MP2 energy correction
     * @details This variable stores the MP2 energy correction calculated during the post-Hartree-Fock method.
     */
    double MP2_energy = MP2_energy_correction(index, value, mo_energy, n_up, n_integrals);
 
    end_mp2 = clock(); // End timing the MP2 energy correction calculation

    printf("Done!\n"); 

    clock_t end_total = clock(); // End timing the entire program

    // Calculate times in seconds
    double time_total = ((double) (end_total - start_total)) / CLOCKS_PER_SEC;
    double time_hf = ((double) (end_hf - start_hf)) / CLOCKS_PER_SEC;
    double time_mp2 = ((double) (end_mp2 - start_mp2)) / CLOCKS_PER_SEC;
    double time_other = time_total - (time_hf + time_mp2);

    // Print a summary
    printf("\n################## Energy Summary ##################\n");
    printf("\nNuclear repulsion energy:            %9.6lf\n", nuc_repul);
    printf("One-electron energy:                 %9.6lf\n", one_el_energy);
    printf("Two-electron energy:                 %9.6lf\n", two_el_energy);
    printf("Hartree-Fock energy:                 %9.6lf\n", HF_energy);
    printf("MP2 energy correction:               %9.6lf\n", MP2_energy);
    printf("Total energy (HF + MP2):             %9.6lf\n", HF_energy + MP2_energy);
    printf("\n################ System Information ################\n");
    printf("\nNumber of occupied orbitals:         %d\n", n_up);
    printf("Number of molecular orbitals:        %d\n", mo_num);
    printf("Number of two-electron integrals:    %ld\n", n_integrals);
    printf("\n################# Timing Information ################\n");
    printf("HF calculation time:                 %.6f seconds\n", time_hf);
    printf("MP2 calculation time:                %.6f seconds\n", time_mp2);
    printf("I/O and setup time:                  %.6f seconds\n", time_other);
    printf("Total execution time:                %.6f seconds\n", time_total);

    // Finalize the calculation
    printf("\nYour calculation is done.\n");
    printf("Thank you for using the program!\n");

    // Free the allocated arrays
    free(mo_energy);
    mo_energy = NULL;
    free(value);
    value = NULL;
    free(index);
    index = NULL;
    free(data);
    data = NULL;

    return 0;
}

