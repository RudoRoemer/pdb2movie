'''
remove_rotamers.py: this is a pymol python script designed to get rid of all rotamers other than the first one at each region with rotamers
'''

import sys
from pymol import cmd
# print(sys.argv)


# the python script is called with an extra command-line argument which is the filename of the PDB structure from which rotamers will be removed
output_filename=sys.argv[3]

# we load the file in question into pymol
cmd.load(output_filename,"pro")

# then, we remove all atoms which have anything but A as "alt" (which is how PDB files track rotamers) - tst only saves the return status of this command
tst=cmd.remove("not alt A+""",quiet=1)
# print(tst)

# since there are no more alternative rotamers, what was left in the PDB file has its "alt" parameter erased
cmd.alter(all, "alt=''")

# finally, we just save the file
cmd.save(output_filename,"pro")
