//Imports

#include <stdio.h>
#include <stdlib.h>
#include <trexio.h>

int main(){

//Open a TREXIO file

const char* filename = "./data/h2o.h5";

trexio_exit_code rc;
trexio_t* trexio_file = trexio_open(filename, 'r', TREXIO_AUTO, &rc);
if (rc != TREXIO_SUCCESS) {
    fprintf(stderr, "TREXIO Error: %s\n", trexio_string_of_error(rc));
    exit(1);
}

//Read the nuclear repulsion energy

double nuc_repul; // Variable where the nuclear repulsion energy is read
rc = trexio_read_nucleus_repulsion(trexio_file, &nuc_repul);
// Check the return code to be sure reading was OK
if (rc != TREXIO_SUCCESS) {
    fprintf(stderr, "TREXIO Error reading nuclear repulsion energy:\n%s\n",
    trexio_string_of_error(rc));
    exit(1);
}

//Obtain the number of occupied orbitals

int32_t n_up; //Variable where the number of spin-up electrons is read
rc = trexio_read_electron_up_num(trexio_file, &n_up);
// Check the return code to be sure reading was OK
if (rc != TREXIO_SUCCESS) {
    fprintf(stderr, "TREXIO Error reading number of spin-up electrons: n%s\n",
    trexio_string_of_error(rc));
    exit(1);
}

//Obtain the number of molecular orbitals

int32_t mo_num; //Variable where the number of molecular orbitals is read
rc = trexio_read_mo_num(trexio_file, &mo_num);
// Check the return code to be sure reading was OK
if (rc != TREXIO_SUCCESS) {
    fprintf(stderr, "TREXIO Error reading number of molecular orbitals: n%s\n",
    trexio_string_of_error(rc));
    exit(1);
}

//Read in the one-electron integrals

int64_t elements = mo_num*mo_num; // Calculate number of elements to store in data array
//Allocate data array
double* data = malloc(elements*sizeof(double)); // Request mo_num x mo_num doubles
if (data == NULL) { // Check that the allocation was OK
    fprintf(stderr, "Allocation of data array for one-electron integrals failed\n");
    exit(-1);
}
//Use TREXIO function to fill data array with one-electron integrals
rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, data);
// Check the return code to be sure reading was OK
if (rc != TREXIO_SUCCESS) {
    fprintf(stderr, "TREXIO Error reading one-electron integrals: n%s\n",
    trexio_string_of_error(rc));
    exit(1);
}

//Calculate one-electron energy contribution

double one_el_energy = 0;
for (int i=0; i < n_up; i++){ // Iterate over the occupied orbitals
    one_el_energy = one_el_energy + 2 * data[i*mo_num+i]; // Get the value the integral and add
    printf("%9.6lf\n", data[i*mo_num+i]);
} 

//Free the data array
free(data);
data = NULL; //Reset pointer

//Read in the two-electron integrals

//Get number of non-zero integrals
int64_t n_integrals;
rc = trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);
// Check the return code to be sure reading was OK
if (rc != TREXIO_SUCCESS) {
    fprintf(stderr, "TREXIO Error reading number of non-zero two-electron integrals: n%s\n",
    trexio_string_of_error(rc));
    exit(1);
}
//Allocate memory for storing the indices
int32_t* index = malloc(4 * n_integrals * sizeof(int32_t));
if (index == NULL) { //Check that the allocation was OK
    fprintf(stderr, "Allocation of index array for two-electron integrals failed\n");
    exit-(1);
}
//Allocate memory for storing the values of the integrals
double* value = malloc(n_integrals * sizeof(double));
if (value == NULL) { // Check that the allocation was OK 
    fprintf(stderr, "Allocation of value array for two-electron integrals failed\n");
    exit(-1);
}
//Read in the integrals
int64_t buffer_size = n_integrals;
rc = trexio_read_mo_2e_int_eri(trexio_file, 0, &buffer_size, index, value);
// Check the return code to be sure reading was OK
if (rc != TREXIO_SUCCESS) {
    fprintf(stderr, "TREXIO Error reading two-electron integrals: n%s\n",
    trexio_string_of_error(rc));
    exit(1);
}
//Check if buffer size is equal to number of integrals
if (buffer_size != n_integrals) {
    fprintf(stderr, "Not all two-electron integrals were read correctly\n");
    exit(1);
}

//Close the TREXIO file

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

// Calculate the electronic interaction energy from the two-electron integrals

double el_inter = 0; //Variable to store the electronic interaction energy
for (int n=0; n < n_integrals; n++) { //Iterate over the stored integrals
    //Get the indices
    int i = index[4*n];
    int j = index[4*n+1];
    int k = index[4*n+2];
    int l = index[4*n+3];
    if (i < n_up && j < n_up) { //Check if the first two indices belong to occupied orbitals
        if (i == j && j == k && k == l) {
            el_inter = el_inter + value[n]; //all indices are the same, only add once
            printf("S: %2d, %2d, %2d, %2d    %9.6lf\n", i, j, k, l, value[n]);
        }
        else if (k == i && l ==j) { 
            el_inter = el_inter + (2 * 2 * value[n]); //Add x2 the Coulomb integral, *2 for permutational symmetry
            printf("J: %2d, %2d, %2d, %2d    %9.6lf\n", i, j, k, l, value[n]);
        }
        else if (i == j && k == l) {
            el_inter = el_inter - (2 * value[n]); //Substract the exchange integral, *2 for permutational symmetry
            printf("K: %2d, %2d, %2d, %2d    %9.6lf\n", i, j, k, l, value[n]);
        }
    }
    else if (l >= n_up) {
        break; //Break the loop once there are only integrals with virtual orbitals in the list
    }    
}

// Calculate the Hartree-Fock energy

double HF_energy = nuc_repul + one_el_energy + el_inter;

//Print a summary
printf("Nuclear repulsion energy:    %9.6lf\n", nuc_repul);
printf("Core Hamiltonian: %9.6lf\n", one_el_energy);
printf("Electronic interaction energy: %9.6lf\n", el_inter);
printf("Hartree-Fock energy: %9.6lf\n", HF_energy);
printf("Number of occupied orbitals:  %d\n", n_up);
printf("Number of molecular orbitals: %d\n", mo_num);
printf("Number of two-electron integrals: %ld\n", n_integrals);

//Free the allocated arrays
free(index);
index = NULL; //Reset pointer? Either skip this add "const" in variable defintion
free(value);
value = NULL;

return 0;
}
