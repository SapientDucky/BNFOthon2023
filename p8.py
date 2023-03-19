'''
Prompt 8
Takes fasta file, alignment file, and Newick phylogenetic tree file.
Replaces the Leaf names with the accession number and organism name
'''


import re
import argparse
import os
import subprocess
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
import stat
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo

def get_phylo(aligned_sequences):
    # Calculate the distance matrix
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(aligned_sequences)

    # Build the tree
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)

    # Print the tree in Newick format
    Phylo.write(tree, "tree.nwk", "newick")

    # Calculate the distance matrix
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(aligned_sequences)

    # Build the tree
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)

    # Print the tree in Newick format
    Phylo.write(tree, "tree.nwk", "newick")
    #print(str(tree))
    return str(tree)


def get_inputs() -> tuple:
    
    fastapath = input("Please enter the file path location of your fasta file: ")
    alignpath = input("Please enter the file path location of your alignment file: ")
    treepath = input("Please enter the file path location of your Newick tree file: ")
    output_path = input("Please enter the file path location of your output: ")
    return (fastapath, output_path, alignpath, treepath)
    
def get_muscle() -> os.PathLike:
    muscle_path = os.path.normpath(os.getcwd() + '/MUSCLE')
    if not os.path.exists(muscle_path):
        os.mkdir(muscle_path)
    os.chdir(muscle_path)
    command = ''
    if os.uname().sysname == 'Linux':
        muscle_name = 'muscle5.1.linux_intel64'
        muscle_path = os.path.join(muscle_path, muscle_name)
        if not os.path.exists(muscle_path):
            command = 'wget "https://github.com/rcedgar/muscle/releases/download/5.1.0/muscle5.1.linux_intel64"'
        
    elif os.uname().sysname == 'Darwin':
        muscle_name = 'muscle5.1.macos_arm64'
        muscle_path = os.path.join(muscle_path, muscle_name)
        if not os.path.exists(muscle_path):
            command = 'wget "https://github.com/rcedgar/muscle/releases/download/5.1.0/muscle5.1.macos_arm64"'
        
    elif os.uname().sysname == 'Windows':
        muscle_name = 'muscle5.1.win64.exe'
        muscle_path = os.path.join(muscle_path, muscle_name)
        if not os.path.exists(muscle_path):
            command = 'wget "https://github.com/rcedgar/muscle/releases/download/5.1.0/muscle5.1.win64.exe"'
    
    else: 
        raise SystemError(f'Unsupperted system {os.uname().sysname}')

    subprocess.run(command, shell=True)
    return muscle_path

def parse_args() -> tuple:
    '''Parse args when file run. (args.fasta, args.output_dir, args.clustal, args.tree)'''

    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description='Phylogenetic Tree Generator', add_help='test')

    # Add arguments
    parser.add_argument(
        '--fasta',
        type=str,
        help='Input a .fasta file from directory.'
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        help='Output newick formatted tree file to directory.'
    )
    parser.add_argument(
        '--clustal',
        type=str,
        help='Input a msa in clustal format.'
    )
    parser.add_argument(
        '--tree',
        type=str,
        help='Input a phylogenetic tree file.'
    )


    # Parse the arguments
    args = parser.parse_args()

    # Return the arguments
    return tuple([args.fasta, args.output_dir, args.clustal, args.tree])

def read_fastafile(filepath):
    """A function to read a fasta file. Returns a list of lines"""
    with open(filepath) as infile:
        fasta_file = infile.readlines()
    return fasta_file

def read_alignmentfile(filepath):
    """A function to read an alignment file. Returns a list of lines"""
    with open(filepath) as infile:
        alignment_file = infile.readlines()
    return alignment_file

def read_treefile(filepath):
    """Read a Newick file in, split by commas. Returns a list."""
    with open(filepath) as infile:
        tree_file = infile.read().split(",")
    return tree_file

def get_seqIds(fasta_file):
    """A function to obtain the sequence IDs from the fasta file lines. Returns a list of sequence IDs."""
    seqIds = []
    for line in fasta_file:
        if line.startswith(">"):
            seqIds.append(line.rstrip().lstrip("> "))
    return seqIds

def get_identpair(seqIds):
    """Accepts sequence ID lines, returns a dictionary where the key:value pairs are accession#:Organism"""
    identpairs = {}  # empty starting dictionary
    orgpattern = "OS=(.*?)\sOX="  # pattern for obtaining organism name
    for i in seqIds:
        key = re.split("\|", i)  # split line around |
        key = key[1]  # accession number is at index 1 of line split around |
        org = re.search(orgpattern, i).group(1)
        identpairs[key] = org
    return identpairs

def replace_Ids(tree_file, identpairs):
    """Accepts a list made from Newick file split around commas, and a dictionary
    with accession:organism pairs. Returns a rejoined string in Newick format."""
    newtree = []
    for line in tree_file:
        for key in identpairs:
            if key in line:
                # join the accession and organism with undercore
                replacement_text = f"{key}_{identpairs[key]}"
                # Replace the text between '(' and ':' with replacement text
                new_str = re.sub(r"(\()?(\w{2}\|).*?:", f"({replacement_text}:", line)
                newtree.append(new_str)  # add to list
    newtree = ",".join(newtree)  # joins list with commas to restore Newick formatting
    return newtree

def run_muscle(input_file:os.PathLike, muscle_path:os.PathLike, output_file:os.PathLike):
    os.chmod(muscle_path, stat.S_IRWXU)
    # Set MUSCLE commandline parameters
    #muscle_cline = MuscleCommandline(muscle_path, input=input_file, out=output_file)

    # Run MUSCLE commandline
    
    #stdout, stderr = muscle_cline()
    cmd = [muscle_path, "-align", input_file, "-output", output_file]
    subprocess.run(cmd)

    # Parse the output alignment file using AlignIO
    alignment = AlignIO.read(output_file, "fasta")

    # Print the alignment summary
    print(alignment)

    # Save the alignment to a file
    AlignIO.write(alignment, "output.fasta", "fasta")
    return alignment

def out_tree(tree_output_path:os.PathLike, tree):

    with open(tree_output_path, "w+", encoding="ascii") as out_file:
        print(tree, file=out_file)

def main():
    args = parse_args()
    if args[0] == None:
        args = get_inputs()
    print(args)
    muscle_path = get_muscle()
    fastapath = args[0]
    out_aln = os.path.join(args[1], 'out.aln')
    #alignpath = out_aln
    treepath = args[3]
    tree_out = os.path.join(args[1], 'out.tree')
    
    alignment = run_muscle(muscle_path=muscle_path, input_file=fastapath, output_file=out_aln)
    fasta_file = read_fastafile(fastapath)
    #alignment_file = read_alignmentfile(alignpath)
    tree_file = read_treefile(treepath)
    seqIds = get_seqIds(fasta_file)
    identpairs = get_identpair(seqIds)
    tree_file = get_phylo(alignment)
    treepath = tree_file
    newtree = replace_Ids(tree_file, identpairs)
    out_tree(tree_output_path=tree_out, tree=newtree)
    # export newtree to file
    

    


if __name__ == "__main__":

    main()