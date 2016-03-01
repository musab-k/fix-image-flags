# fix-image-flags
Generates image flags for atoms that are part of bonds that span periodic boundaries in a LAMMPS data file. 

method:
    - reads in LAMMPS data file
    - looks for bonds that span periodic boundaries
    - adds image flags
    - writes new LAMMPS data file

note: 
    - completely ignores existing image flags
    - assumes 'atom_style full' 
    - requires pizza.py (pizza.sandia.gov)

usage:
    - python fix_image_flags.py <initial_data_file> <new_data_file>
