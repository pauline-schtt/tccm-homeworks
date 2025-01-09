# Molecular Dynamics Simulation Program 

This project contains a program for running molecular dynamics (MD) simulations. The program reads atomic coordinates and masses from an input file, performs MD simulations, and outputs the trajectory file in XYZ format. The trajectory file contains atomic coordinates, as well as the values of kinetic, potential, and total energies.

## Project Structure

This project is composed of several folders. For instance, the `data` folder has the files for testing purposes, the `src` folder contains the source code, and the rest consists of files like `AUTHORS`, `INSTALL.md` with the installation instructions, `LICENSE`, as well as the `Makefile` used for the compilation process.

```
└── 📁project3
    └── 📁data
        └── inp.txt
        └── trajectory.xyz
    └── 📁src
        └── functions.c
        └── headers.h
        └── main.c
    └── AUTHORS
    └── dynamics.pdf
    └── INSTALL.md
    └── LICENSE
    └── Makefile
    └── README.md
```

## Installation

The installation steps can be found in the INSTALL.md file.

## Usage

To run the molecular dynamics simulation, provide the full path to the input file containing the atomic coordinates and masses as an argument to the program. The example of the input file can be found in `data/inp.txt`.

Example:
```sh
./MD data/inp.txt
```