import os
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# Load topology and trajectory
topology_file = "step3_input.pdb"
trajectory_file = "step5_production.xtc"
traj = md.load(trajectory_file, top=topology_file)

# Select the UNK chain
unk_chain_indices = traj.topology.select("resname UNK")
unk_chain_traj = traj.atom_slice(unk_chain_indices)

# Select the rest of the protein
protein_indices = traj.topology.select("protein")
protein_traj = traj.atom_slice(protein_indices)

# Create pairs of atoms for which to compute distances
pairs = np.array([[i, j] for i in unk_chain_indices for j in protein_indices])

# Compute distances
distances = md.compute_distances(traj, pairs)

# Threshold for interaction (in nanometers)
threshold = 0.5

# Initialize interaction counter
interaction_counter = {}

# Iterate over each frame
for i in range(traj.n_frames):
    # Get the distances for this frame
    frame_distances = distances[i]

    # Find interactions (distances less than threshold)
    interactions = np.where(frame_distances < threshold)[0]

    # Iterate over each interaction
    for interaction in interactions:
        # Get the residue names
        res1 = traj.topology.atom(pairs[interaction][0]).residue
        res2 = traj.topology.atom(pairs[interaction][1]).residue

        # Skip if both residues are UNK
        if res1.name == 'UNK' and res2.name == 'UNK':
            continue

        # Create a key for the interaction
        key = (res1, res2)

        # Increment the counter for this interaction
        if key in interaction_counter:
            interaction_counter[key] += 1
        else:
            interaction_counter[key] = 1

# Sort interactions by count
sorted_interactions = sorted(interaction_counter.items(), key=lambda x: x[1], reverse=True)

# Create a figure and axis for the plot
fig, ax = plt.subplots(figsize=(18, 12))  # Adjust the size of the plot

# Get interactions and counts
interactions = [f"{interaction[0].name} ({interaction[0].index}) - {interaction[1].name} ({interaction[1].index+1})" for interaction, count in sorted_interactions]
# Divide counts by 200
counts = [count / 200 for interaction, count in sorted_interactions]

# Create a bar plot with increased font size
bar_plot = ax.bar(interactions, counts, color='red')  # Change the color of the bars

# Set labels and title with increased font size
ax.set_xlabel("Interaction", fontsize=16)
ax.set_ylabel("Simulation time (ns)", fontsize=16)
ax.set_title("Interactions over Time", fontsize=16)

# Rotate x-axis labels for better visibility with increased font size
plt.xticks(rotation=90, fontsize=12)

# Set font size for y-axis ticks
plt.yticks(fontsize=16)

# Set font size for bar labels
for rect in bar_plot:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width() / 2, height, f'{height:.2f}', ha='center', va='bottom', rotation=45, fontsize=10)

# Save the plot as a PNG file
plt.tight_layout()
plt.savefig('Simulation_time_plot.png')  # Provide the desired file name and extension    
    
# Show the plot
plt.tight_layout()
plt.show()