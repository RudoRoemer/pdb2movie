import sys
from pymol import cmd
print(sys.argv)
output_filename=sys.argv[3]
cmd.load(output_filename,"pro")
tst=cmd.remove("not alt A+""",quiet=1)
print(tst)
cmd.alter(all, "alt=''")
cmd.save(output_filename,"pro")
