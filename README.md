Input file;
1. PDB file (make sure that your PDB have separate chain name for protein and ligand. For instance, protein = chain name "A", "B", ligand = chain name "UNK")
2. GROMACS xtc file (MD simulation output file)

In case if use of different ligand chain name, you need to optimized this part:
unk_chain_indices = traj.topology.select("resname UNK")
Change "resname UNK" into your ligand chain name.
