#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <trexio.h>
#include "headers.h" // Function headers

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
    trexio_exit_code rc; //! Varibale for storing the TREXIO output
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

    // Allocate array for orbital energies
    double* mo_energy = malloc(mo_num * sizeof(double)); //! Array for storing the orbital energies
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

    int64_t elements = mo_num*mo_num; //! Variable for the number of elements to store in data array
    // Allocate data array
    double* data = malloc(elements*sizeof(double)); //! Request mo_num x mo_num doubles for storing the one-electron integrals
    if (data == NULL) { // Check that the allocation was OK
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
    int64_t n_integrals; //! Variable for number of non-zero two-electron integrals
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
    int32_t* index = malloc(4 * n_integrals * sizeof(int32_t)); //! Array for storing the indices of the two-electron integrals
    if (index == NULL) { // Check that the allocation was OK
        fprintf(stderr, "Allocation of index array for two-electron integrals failed\n");
        free(mo_energy);
        free(data);
        exit(-1);
    }

    // Allocate memory for storing the values of the integrals
    double* value = malloc(n_integrals * sizeof(double)); //! Array for storing the values of the two-electron integrals
    if (value == NULL) { // Check that the allocation was OK 
        fprintf(stderr, "Allocation of value array for two-electron integrals failed\n");
        free(mo_energy);
        free(data);
        free(index);
        exit(-1);
    }

    // Read in the integrals
    int64_t buffer_size = n_integrals; //! Variable buffer size for reading the two-electron integrals, equal to n_integrals
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
    
    //Calculate one-electron energy contribution
    double one_el_energy = one_electron_energy(data, n_up, mo_num); //! Variable for the one-electron energy contribution

    // Calculate the two-electron energy contribution
    double two_el_energy = two_electron_energy(index, value, n_up, n_integrals); //! Variable for the two-electron energy contribution
    
    // Calculate the Hartree-Fock energy
    double HF_energy = hartree_fock_energy(nuc_repul, one_el_energy, two_el_energy); //! Variable for the overall Hartree-Fock energy
    
    printf("Done!\n");

    printf("\nCalculating the MP2 energy correction...\n");
    
    // Calculate MP2 energy
    double MP2_energy = MP2_energy_correction(index, value, mo_energy, n_up, mo_num, n_integrals); //! Variable for the MP2 energy
    
    printf("Done!\n"); 

    // Print a summary
    printf("\n################## Energy Summary ##################\n");
    printf("\nNuclear repulsion energy:           %9.6lf\n", nuc_repul);
    printf("One-electron energy:                 %9.6lf\n", one_el_energy);
    printf("Two-electron energy:                 %9.6lf\n", two_el_energy);
    printf("Hartree-Fock energy:                 %9.6lf\n", HF_energy);
    printf("MP2 energy correction:               %9.6lf\n", MP2_energy);
    printf("Total energy (HF + MP2):             %9.6lf\n", HF_energy + MP2_energy);
    printf("\n################ System Information ################\n");
    printf("\nNumber of occupied orbitals:         %d\n", n_up);
    printf("Number of molecular orbitals:        %d\n", mo_num);
    printf("Number of two-electron integrals:    %ld\n", n_integrals);

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
