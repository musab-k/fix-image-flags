# fix-image-flags
Generates image flags for atoms that are part of bonds that span periodic boundaries in a LAMMPS data file. 

Method: Reads in LAMMPS data file, looks for bonds that span periodic boundaries, adds image flags and writes new LAMMPS data file

Note: completely ignores existing image flags, assumes 'atom_style full' and requires pizza.py (pizza.sandia.gov)

Usage: python fix_image_flags.py initial_data_file new_data_file
