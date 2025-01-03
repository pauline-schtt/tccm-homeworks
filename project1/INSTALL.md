## Required software

Ensure you have the following installed on your system:
- `gcc` (version 14.2.0 and ... were tested)
- `make` (version 3.81 and ... were tested)
- `hdf5`
- `trexio` (version 2.5.0 was tested)

## Installation steps

1. First, clone the repository:
    ```sh
    git clone https://github.com/pauline-schtt/tccm-homeworks/tree/master/project1
    cd project1
    ```

2. Compile the program using `make`:
    ```sh
    make
    ```

3. And finally, run the program:
    ```sh
    ./HF_and_MP2 <path_to_the_input_file>
    ```

## Test

The program was tested on MacOS Sequoia 15.2 and Ubuntu ... . To test whether the installation was succesful, you could calculate the Hartree-Fock and MP2 energies on the set of molecules provided in the `data` folder. The reference values are given in the `README.org` file.

## Cleaning up

To clean up the compiled files, run:
```sh
make clean
```

## Troubleshooting

If you encounter any issues during the installation, ensure that your versions of the prerequisites meet the required versions mentioned above. In case you're still facing issues, feel free to reach out to any of the contributors or opening an issue on GitHub.