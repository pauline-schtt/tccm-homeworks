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

    char* filename; //! Name of the HDF5 file

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
    trexio_exit_code rc; //! TREXIO output
    trexio_t* trexio_file = trexio_open(filename, 'r', TREXIO_AUTO, &rc); //! TREXIO file handler
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error: %s\n", trexio_string_of_error(rc));
        exit(1);
    }

    // Read the nuclear repulsion energy
    double nuc_repul; //! Nuclear repulsion energy
    rc = trexio_read_nucleus_repulsion(trexio_file, &nuc_repul);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading nuclear repulsion energy:\n%s\n",
                trexio_string_of_error(rc));
        exit(1);
    }

    // Obtain the number of occupied orbitals
    int32_t n_up; //! Number of spin-up electrons
    rc = trexio_read_electron_up_num(trexio_file, &n_up);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading number of spin-up electrons: n%s\n",
                trexio_string_of_error(rc));
        exit(1);
    }

    // Obtain the number of molecular orbitals
    int32_t mo_num; //! Number of molecular orbitals
    rc = trexio_read_mo_num(trexio_file, &mo_num);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading number of molecular orbitals: n%s\n",
                trexio_string_of_error(rc));
        exit(1);
    }

    // Allocate array for orbital energies
    double* mo_energy = malloc(mo_num * sizeof(double)); //! Array of molecular orbital energies
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

    int64_t elements = mo_num * mo_num; //! Number of elements in the data array

    double* data = malloc(elements * sizeof(double)); //! Array of one-electron integrals
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
    int64_t n_integrals; //! Number of non-zero two-electron integrals
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
    int32_t* index = malloc(4 * n_integrals * sizeof(int32_t)); //! Array of indices of the two-electron integrals
    if (index == NULL) {
        fprintf(stderr, "Allocation of index array for two-electron integrals failed\n");
        free(mo_energy);
        free(data);
        exit(-1);
    }

    // Allocate memory for storing the values of the integrals
    double* value = malloc(n_integrals * sizeof(double)); //! Array of values of the two-electron integrals
    if (value == NULL) {
        fprintf(stderr, "Allocation of value array for two-electron integrals failed\n");
        free(mo_energy);
        free(data);
        free(index);
        exit(-1);
    }

    // Read in the integrals
    int64_t buffer_size = n_integrals; //! Buffer size for reading, initially equal to number of two-electron integrals
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
    double one_el_energy = one_electron_energy(data, n_up, mo_num); //! One-electron energy contribution

    // Calculate the two-electron energy contribution
    double two_el_energy = two_electron_energy(index, value, n_up, n_integrals); //! Two-electron energy contribution
    
    // Calculate the Hartree-Fock energy
    double HF_energy = hartree_fock_energy(nuc_repul, one_el_energy, two_el_energy); //! Hartree-Fock energy

    end_hf = clock(); // End timing the Hartree-Fock energy calculation
    
    printf("Done!\n");

    printf("\nCalculating the MP2 energy correction...\n");
    
    start_mp2 = clock(); // Start timing the MP2 energy correction calculation

    // Calculate MP2 energy
    double MP2_energy = MP2_energy_correction(index, value, mo_energy, n_up, n_integrals); //! MP2 energy
 
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

