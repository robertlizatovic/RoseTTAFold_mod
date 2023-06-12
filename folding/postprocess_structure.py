import pyrosetta
import multiprocessing as mp
import argparse

# Initialize PyRosetta with the desired options
pyrosetta.init("""
    -mute all
    -relax:constrain_relax_to_start_coords
    -ex1
    -ex2
    -use_input_sc
    -detect_disulf true
    -detect_disulf_tolerance 3.5
    -preserve_crystinfo true
""")

# Create a score function
scorefxn = pyrosetta.get_fa_scorefxn()

# Set up the relax protocol
def run_relax_trajectory(pdb_file:str):
    """
    Run the FastRelax protocol to pack SC and minimize
    the BB and SCs
    """
    # Create a pose object and load the PDB file
    pose = pyrosetta.pose_from_pdb(pdb_file)
    relax = pyrosetta.rosetta.protocols.relax.FastRelax()
    relax.set_scorefxn(scorefxn)
    relax.constrain_relax_to_start_coords(True)  # keep close to the starting coordinates
    #relax.set_auto_relax_disulf(True) # Enable automatic disulfide bond formation

    # Run the relax protocol and return the minimized pose
    relax.apply(pose)
    return pose

def propagate_confidences(pose) -> None:
    """
    Propagates the pLDDT scores of each CA atom to all
    the corresponding SC atoms (like in AF2) in the passed pose
    """
    # Iterate over each residue in the pose
    for i in range(1, pose.total_residue() + 1):
        # Get the B-factor value of the C-alpha atom
        residue = pose.residue(i)
        #print(residue.name3())
        if residue.name3() != "XXX":
            ca_index = residue.atom_index("CA")
            b_factor = round(pose.pdb_info().bfactor(i, ca_index) * 100, 2)
            #print(i, b_factor)

            # Assign the B-factor value to all atoms in the residue
            for j in range(1, residue.natoms() + 1):
                pose.pdb_info().bfactor(i, j, b_factor)
                #print(i, j, pose.pdb_info().bfactor(i, j))

def relax(pdb_file:str, nstruct:int, nproc:int):
    """
    Performs several MC FastRelax runs in parallel and
    returns the lowest energy structure obtained
    """
    assert nstruct > 0, "nstruct must be at least 1!"
    assert nproc > 0, "nproc must be at least 1!"

    # Run relax trajectories in parallel
    if nproc > 1 and nstruct > 1: 
        workers = min(nstruct, nproc)
        pool = mp.Pool(processes=workers)
        poses = pool.map(run_relax_trajectory, [pdb_file for _ in range(nstruct)])
    else:
        poses = [run_relax_trajectory(pdb_file) for _ in range(nstruct)]

    # Find the pose with the lowest energy
    best_pose = min(poses, key=lambda pose: scorefxn(pose))
    return best_pose

def main(args) -> None:
    """
    Performs post-processing of RoseTTAFold end2end structure prediction pipeline by:
    1. Performing several FastRelax MC trajectories in parallel
    2. Taking the lowest energy structure
    3. Propagating the BB pLDDT scores to all SC atoms
    """
    pose = relax(args.input, args.nstruct, args.nproc)
    propagate_confidences(pose)
    # Output the lowest energy pose to a PDB file
    pose.dump_pdb(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, required=True, dest="input", help="path to input pdb file")
    parser.add_argument("-o", type=str, required=True, dest="output", help="path to the output pdb file")
    parser.add_argument("-n", type=int, required=False, default=10, dest="nstruct", help="number of relax trajectories to run")
    parser.add_argument("-p", type=int, required=False, default=10, dest="nproc", help="number of processes to use")
    args = parser.parse_args()
    main(args)
