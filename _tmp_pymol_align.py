
from pymol import cmd
import sys

ref, mob, outpdb, outtxt = sys.argv[1:5]

cmd.load(ref, "ref")
cmd.load(mob, "mob")

# backbone atoms only (paper wording)
r = cmd.align("mob and polymer.protein and backbone",
              "ref and polymer.protein and backbone")

cmd.save(outpdb, "mob")

with open(outtxt, "w") as f:
    f.write(str(r[0]))   # RMSD only

cmd.quit()
