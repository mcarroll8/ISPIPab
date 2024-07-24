# Evan Edelstein
import pymol2
import sys
"""
1) color TPs as deep blue 

"""

with pymol2.PyMOL() as p1:
    cmd = p1.cmd
    output_path_dir = sys.argv[1]
    pdbid = sys.argv[2]
    tp_residues = sys.argv[3]
    fp_residues = sys.argv[4]
    annotated_resiues = sys.argv[5]
    pred = sys.argv[6]

    cmd.fetch(f"{pdbid}")
    cmd.orient(f"{pdbid}")

    cmd.color("orange")
    cmd.set("cartoon_transparency", "0.75")
    cmd.select("ann", f"resi {annotated_resiues}")
    cmd.indicate("bycalpha ann")
    cmd.create("annotated", "indicate")

    cmd.show("sphere", "annotated")
    cmd.color("pink", "annotated")
    cmd.set("sphere_transparency", "0.4", "annotated")

    cmd.select("tp_residues", f"resi {tp_residues}")
    cmd.indicate("bycalpha tp_residues")
    cmd.create("tp_residues", "indicate")

    cmd.select("fp_residues", f"resi {fp_residues}")
    cmd.indicate("bycalpha fp_residues")
    cmd.create("fp_residues", "indicate")

    cmd.show("sphere", "tp_residues")
    cmd.set("sphere_scale", "0.5", "tp_residues")
    cmd.color("green", "tp_residues")
    cmd.set("sphere_transparency", "0", "tp_residues")
    cmd.set("cartoon_transparency", "1", "tp_residues")

    cmd.show("sphere", "fp_residues")
    cmd.set("sphere_scale", "0.5", "fp_residues")
    cmd.color("deepblue", "fp_residues")
    cmd.set("sphere_transparency", "0.45", "fp_residues")
    cmd.set("cartoon_transparency", "1", "fp_residues")

    cmd.remove("resn hoh and not polymer.protein")
    cmd.bg_color("white")
    cmd.zoom(complete=1)
    cmd.save(f"{output_path_dir}/proteins/{pdbid}/pymol_{pred}.pse")
    # cmd.png(f"{output_path_dir}/{pdbid}/pymol_viz_{pred}.png")
