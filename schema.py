#input: fasta file with sequences of both parents, pdb file of one parent


#RASPP??
#n_cuts=1? proc vic, co s tim obecne

from data_structures import PdbProtChain
from data_structures import ProtChain
from data_structures import Atom
from data_structures import Residue
from data_structures import Point
from pdb_structure_parser import parse_pdb_file
import sys
import tempfile
import subprocess
from pathlib import Path
from math import sqrt

# function returns list of sequences contained in given fasta file
def get_chains_from_fasta_file(file_name: str)-> list[ProtChain]:
    with open(file_name,"r",encoding="utf-8") as file:
        content = file.read()

    return get_chains_from_fasta_string(content)
    
    

#function returns array of ProtChains from fasta string
def get_chains_from_fasta_string(fasta_string: str)-> list[ProtChain]:
    lines = fasta_string.splitlines()
    chains = []
    i=0

    while i <len(lines):
        #start of a new sequence
        if(lines[i].startswith(">")):
            chain_id = lines[i]
            i+=1
            sequence = ""

            while i < len(lines):
                line = lines[i].strip()

                if(line.startswith(">")): break

                if(line != ""): sequence = sequence + lines[i]

                i +=1

            chain = get_ProtChain_strucutre(chain_id,sequence)
            chains.append(chain)
            
        else:
            i+=1

    return chains

# function returns ProtChain structure with filled values
def get_ProtChain_strucutre(name_line: str, sequence_line: str)-> ProtChain:
    name = name_line.lstrip(">").strip()
    sequence = sequence_line.strip()
    return ProtChain(id=name, sequence = sequence)


# function takes parent chains and structure chain and makes 
def do_msa(parentChains: list[ProtChain], pdbChain: PdbProtChain) -> int | None:
    #make text which will be writte into mafft input file
    fasta_text = get_fasta_string(parentChains,pdbChain)
    if(fasta_text is None): return None

    mafft_result = run_mafft(fasta_text)
    if(mafft_result is None): return None

    #parse chains out of mafft fasta output
    chains = get_chains_from_fasta_string(mafft_result)

    #update structures with aligned sequences
    for i in range(0,len(chains)):
        if(chains[i].id == pdbChain.id):
            pdbChain.aligned_sequence = chains[i].aligned_sequence

        if(chains[i].id == parentChains[0].id):
            parentChains[0].aligned_sequence = chains[i].aligned_sequence

        if(chains[i].id == parentChains[1].id):
            parentChains[1].aligned_sequence = chains[i].aligned_sequence

    return 0

# function returns string reprezentation of given sequences in fasta format
def get_fasta_string(parent_chains: list[ProtChain], pdb_chain: PdbProtChain) -> str | None:
    #get all IDs
    ids = []
    ids.append(pdb_chain.id)
    for i in range(0,len(parent_chains)): ids.append(parent_chains[i].id)
    #get all protein chains
    chains = []
    chains.append(pdb_chain.sequence)
    for i in range(0,len(parent_chains)): chains.append(parent_chains[i].sequence)

    if(len(ids) != len(chains)):
        print("Error: Number of protein IDs is not equal to number of protein chains\n")
        return None
    
    fasta_text = ""
    #make fasta text from IDs and chains
    for i in range(0,len(ids)):
        id_line = ">" + ids[i] + "\n"
        sequence_line = chains[i] +  "\n"
        fasta_text = fasta_text + id_line + sequence_line

    return fasta_text

#function returns text result of mafft tool
def run_mafft(fasta_text: str):
    #make tmp dir 
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)

        input_file_path = temp_dir / "input.fasta"
        output_file_path = temp_dir / "output.fasta"
        #write to input file
        input_file_path.write_text(fasta_text,encoding="utf-8")

        #run mafft and save the result to out file
        with open(output_file_path,"w",encoding="utf-8") as out:
            subprocess.run(
                ["mafft","--auto",str(input_file_path)],
                stdout=out,
                check=True
            )
        #return output file content
        return output_file_path.read_text(encoding="utf-8")
    
    print("Error: mafft or tempfile. failed\n")
    return None

# return map of: residue in alignment -> pdb sequence index
def make_alignment_to_pdb_mapping(alignment_matrix: list[list[str]]) -> list[int | None] | None:
    map=[]
    aligned_sequence = alignment_matrix[0]
    index = 0

    for char in aligned_sequence:
        if(char == "-"):
            map.append(None)
        else:
            map.append(index)
            index+=1
    
    return map

def get_aligned_contact_matrix(pdb_chain: PdbProtChain, mapping: list[int]) -> list[list[float]] | None :
    dimension = len(pdb_chain.aligned_sequence)
    aligned_contact_matrix = [[0.0 for _ in range(dimension)] for _ in range(dimension)]
    aligned_sequence = pdb_chain.aligned_sequence

    #get not aligned contact matrix
    pdb_contact_matrix = get_pdb_contact_matrix(pdb_chain)
    if(pdb_contact_matrix is None): return None

    #transform not aligned contact matrix to aligned contact matrix
    #loop rows
    for i in range(len(aligned_sequence)):
        #loop columns
        for j in range(i+1, len(aligned_sequence)):
            if(aligned_sequence[i]=="-" or aligned_sequence[j]=="-"): continue
            else:
                #get row value in not aligned matrix
                m = mapping[i]
                #get column value in not aligned matrix
                n = mapping[j]
                #assign
                aligned_contact_matrix[i][j] = pdb_contact_matrix[m][n]
                aligned_contact_matrix[j][i] = pdb_contact_matrix[m][n]

    return aligned_contact_matrix

def get_pdb_contact_matrix(pdb_chain: PdbProtChain)-> list[list[float]] | None:
    #init matrix with zeros
    dimension = len(pdb_chain.residues)
    contact_matrix = [[0.0 for _ in range(dimension)] for _ in range(dimension)]

    #loop rows
    for i in range(dimension):
        #loop columns
        for j in range(i+1,dimension):
            value = get_contact_value(
                pdb_chain.residues[i],
                pdb_chain.residues[j]
                )
            
            if(value is None): 
                return None
            
            contact_matrix[i][j] = value
            contact_matrix[j][i] = value
        
        
    return contact_matrix

#return float from <0,1> depending on number of contacts
#1 - all atoms are in contact
#0 - no atoms are in contact
def get_contact_value(res1: Residue, res2: Residue) -> float | None:
    atoms1 = remove_hydrogen_from_list(res1.atoms)
    atoms2 = remove_hydrogen_from_list(res2.atoms)

    normalization_denominator = len(atoms1) * len(atoms2)
    
    if(normalization_denominator == 0):
        print("Error: Residue can not have no atoms\n")
        return None

    contacts = 0

    for atom1 in atoms1:
        for atom2 in atoms2:
            distance = compute_distance_3D(atom1.coords,atom2.coords)
            if(distance <= 4.5): contacts +=1

    return contacts/normalization_denominator    

#functions removes hydrogens from list of atoms and returns a new list
def remove_hydrogen_from_list(atoms: list[Atom])-> list[Atom]:
    new_list = []

    for atom in atoms:
        if(atom.element!="H"):
            new_list.append(atom)

    return new_list

#function computes distance between two points in 3D space
def compute_distance_3D(point_1: Point, point_2: Point) -> float:
    x_distance_squared = (point_1.x - point_2.x)**2
    y_distance_squared = (point_1.y - point_2.y)**2
    z_distance_squared = (point_1.z - point_2.z)**2

    return sqrt(x_distance_squared+y_distance_squared+z_distance_squared)

#function sets irrelevant contacts to zero
def remove_irrelevant_contacts(aligned_contact_matrix: list[list[float]], par1: ProtChain, par2: ProtChain) -> None:
    par1_aligned = list(par1.aligned_sequence)
    par2_aligned = list(par2.aligned_sequence)
    
    for i in range(len(aligned_contact_matrix)):
        for j in range(i+1,len(aligned_contact_matrix[i])):
            should_be_zero = False

            if(par1_aligned[i] == par2_aligned[i]): should_be_zero = True

            if(par1_aligned[j] == par2_aligned[j]): should_be_zero = True

            if(par1_aligned[i] == "-" or par2_aligned[i] == "-"): should_be_zero = True

            if(par1_aligned[j] == "-" or par2_aligned[j] == "-"): should_be_zero = True

            if(should_be_zero):
                aligned_contact_matrix[i][j] = 0.0
                aligned_contact_matrix[j][i] = 0.0

#function transforms matrix to array with non zero contact values only
# template of item: [i index, j index, weigth]
def transform_contact_matrix(matrix: list[list[float]]) -> list[list[int|float]]:
    new_arr = []

    for i in range(len(matrix)):
        for j in range(i+1,len(matrix[i])):
            if(matrix[i][j] != 0):
                new_arr.append([i,j,matrix[i][j]])

    return new_arr

#function returns array with every index evaluated
def get_schema_profile(contacts: list[list[int|float]], aligned_seq_length: int) -> list[float]:
    cut_index_scores = []
    
    for i in range(0,aligned_seq_length-1):
        score = 0.0

        for j in range(len(contacts)):
            contact = contacts[j]
            if(is_between(contact[0],contact[1],i)): score += contact[2]

        cut_index_scores.append(score)

    return cut_index_scores

def get_best_split_index(indides_scores: list[float], map: list[int|None])-> int:
    best_score_index = None

    for i in range(len(indides_scores)):
        if(map[i] is None): continue

        if(best_score_index is None): 
            best_score_index = i
            continue

        if(indides_scores[i] < indides_scores[best_score_index]):
            best_score_index = i

    return map[best_score_index]


def is_between(num1: int, num2: int, value:int) -> bool:
    minimum = min(num1,num2)
    maximum = max(num1,num2)

    if(minimum<=value and value < maximum): return True
    else: return False

def get_split_indices(pdb_chain: PdbProtChain, parent_chains: list[ProtChain]) -> list[int] | None:
    return_value = do_msa(parent_chains, pdb_chain)
    if(return_value is None):  return None

    alignment_matrix = [
        list(pdb_chain.aligned_sequence),
        list(parent_chains[0].aligned_sequence),
        list(parent_chains[1].aligned_sequence),
    ]

    alignment_to_pdb_map = make_alignment_to_pdb_mapping(alignment_matrix)

    aligned_contact_matrix = get_aligned_contact_matrix(pdb_chain,alignment_to_pdb_map)

    remove_irrelevant_contacts(aligned_contact_matrix,parent_chains[0],parent_chains[1])

    transformed_contact_matrix = transform_contact_matrix(aligned_contact_matrix)

    indices_split_scores = get_schema_profile(transformed_contact_matrix,len(pdb_chain.aligned_sequence))

    #remove list after implementation of multi-cut functionality
    return list(get_best_split_index(indices_split_scores,alignment_to_pdb_map))

if __name__ == "__main__":
    sys.exit(get_split_indices())