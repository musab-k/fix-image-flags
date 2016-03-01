from data import *      
import numpy as np
import sys, math

""" fix_image_flags.py
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
"""


def sgn(a):
    return int((a<0)*-1 + (a>=0))


# Check that all bonds are no more than 'threshold' units long
def distance_check(i, j, flags, positions, box, threshold=3):
    r = 0.
    for ii in range(3):
        r += ((positions[i][ii]+flags[i][ii]*box[ii])-(positions[j][ii]+flags[j][ii]*box[ii]))**2
    if math.sqrt(r) > threshold:
        print 'bond %s:%s, |r| = %s' % (i, j, math.sqrt(r))
    return


# computes r_parent - r_child, in unwrapped co-ordinates, and sets flag of child
def correct_flag(parent, child, flags, positions, box):
    for i in range(3):
        #r = (positions[parent][i]+1*flags[parent][i]*box[i]) - (positions[child][i]+1*flags[child][i]*box[i])
        r = (positions[parent][i] - positions[child][i])
        if abs(r) > 0.5*box[i]:
            flags[child][i] = flags[parent][i] + sgn(r)
            #print 'modified: ', parent, child
        else:
            flags[child][i] = flags[parent][i]

    return



def main(infile, outfile):
    d = data(infile)

    # Get box co-ordinates and lengths
    box = [d.headers['xlo xhi'], d.headers['ylo yhi'], d.headers['zlo zhi']]
    box_len = [b[1]-b[0] for b in box]

    # Get atomic positions (assumes 
    positions = [np.array([a for a in ali.strip().split()[4:7]],dtype=float) for ali in d.sections['Atoms']]
    for pi in positions:
        for i in range(3):
            # Assume all initial co-ordinates are unwrapped
            try:
                assert pi[i] >= box[i][0] and pi[i] <= box[i][1]
            except:
                print 'atom outside box:', pi, i, box[i]
                sys.exit()
    natoms = len(positions)

    # Convert bonds to list
    bonds = [[int(b) for b in bli.strip().split()[2:4]] for bli in d.sections['Bonds']]

    # Compute neighbours for every atom
    atom_neighbours = [[] for i in range(natoms)]
    for bi, bj in bonds:
        atom_neighbours[bi-1].append(bj-1)
        atom_neighbours[bj-1].append(bi-1)

    # Set initial flags to 0
    flags = [[0,0,0] for i in range(natoms)]

    mol_ids = {}
    for i in range(len(d.sections['Atoms'])):
        li = d.sections['Atoms'][i]
        mol_no = int(li.strip().split()[1])
        if mol_no not in mol_ids.keys():
            mol_ids[mol_no] = i

    #parent = [0]
    parent = mol_ids.values()
    steps = 0
    while True:
        temp = []
        for p in parent:
            #ptype = d.sections['Atoms'][p].strip().split()[-1]
            for child in atom_neighbours[p]:
                # fix the flag
                correct_flag(p, child, flags, positions, box_len)

                # append child to temp for next iteration
                temp.append(child)

            # kill atom_neighbours
            atom_neighbours[p] = []
        parent = temp
        steps +=1
        if np.size(temp) == 0:
            break

    print 'doing distance check...'
    for bi, bj in bonds:
        distance_check(bi-1, bj-1, flags, positions, box_len)


    # Append flags to atom section in lammps data object
    for i in range(natoms):
        line = d.sections['Atoms'][i]
        l = [li for li in line.strip().split() if '#' not in li]
        assert len(l) == 7  or len(l) == 10 # atom-num, mol-num, atom-type, charge, x, y, z, (image flags)
        l = l[:7] + [str(f) for f in flags[i]]
        nl = "%-10s %8s %8s %10s %15s %15s %15s %8s %8s %8s\n" % tuple(l)
        d.sections['Atoms'][i] = nl

    d.write(outfile)


if __name__ == "__main__":
    assert len(sys.argv) == 3
    infile = sys.argv[1]
    outfile  = sys.argv[2]
    main(infile, outfile)
