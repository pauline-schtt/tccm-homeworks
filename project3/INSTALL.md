## Required software

Ensure you have the following installed on your system:
- `gcc` (version 14.2.0 and 13.3.0 were tested)
- `make` (version 3.81 and 4.3 were tested)

## Installation steps

1. First, clone the repository:
    ```sh
    git clone https://github.com/pauline-schtt/tccm-homeworks/tree/master/project3
    cd project3
    ```

2. Compile the program using `make`:
    ```sh
    make
    ```

3. And finally, run the program:
    ```sh
    ./MD <path_to_the_input_file>
    ```

## Test

The program was tested on MacOS Sequoia 15.2 and Ubuntu 24.04. You can find an example input file in the `data` folder, while the corresponding output files can be found in the `test` folder. As no random velocity initialization is performed, you should obtain qualitatively similar results when running the program with the test input.

## Cleaning up

To clean up the compiled files, run:
```sh
make clean
```

## Troubleshooting

If you encounter any issues during the installation, ensure that your versions of the prerequisites meet the required versions mentioned above. In case you're still facing issues, feel free to reach out to any of the contributors or opening an issue on GitHub.
