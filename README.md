# fix-image-flags
When trying to simulate molecules in periodic boundary conditions using the molecular dynamics package, LAMMPS, it is useful to keep track of which "periodic image" a given part of the molecule belongs to when unwrapped - especially for visualisation purposes. This is done using image flags.

This code aims to correct the problem of image flags that have been set incorrectly, a problem that can go unnoticed as it typically does not affect dynamics. 


Method: Reads in LAMMPS data file, looks for bonds that span periodic boundaries, adds image flags and writes new LAMMPS data file

NB: completely ignores existing image flags, assumes 'atom_style full' and requires pizza.py (pizza.sandia.gov)

Usage: python fix_image_flags.py initial_data_file new_data_file
