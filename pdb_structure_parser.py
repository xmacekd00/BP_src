from data_structures import Atom
from data_structures import Residue
from data_structures import Point
from data_structures import PdbProtChain
import sys

res_name_to_char = {
    "ALA": "A",
    "ARG": "R",
    "ASP": "D",
    "ASN": "N",
    "CYS": "C",
    "GLU": "E",
    "GLN": "Q",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

def parse_pdb_file(folder_base_path: str) -> PdbProtChain | None:
    file_name = folder_base_path + "asr/structure.pdb"
    #get all lines from pdb file
    with open(file_name,"r",encoding="utf-8") as file:
        lines = file.readlines()

    chain = PdbProtChain(id="",sequence="",residues=[])
    line_index = get_first_atom_line_index(lines)
    if(line_index==-1): 
        print("Error: Pdb file does not contain any ATOM records\n")  
        return None
    #index of residue starting from 0
    res_seq_index = 0

    #loop through residues
    while(line_index<len(lines) and lines[line_index].startswith("ATOM")):
        res = Residue()
        line = lines[line_index]
        chain.id = get_chain_id(line)
        #init residue
        res.seq_index= res_seq_index
        res.pdb_index = get_residue_pdb_index(line)
        res.pdb_insertion_code=get_insertion_code(line)
        res.res_name=get_residue_name(line)
        res.res_char=res_name_to_char[res.res_name]
        #get list of atoms contained in residue
        atoms_index_tuple = get_res_atoms(line_index,lines)
        res.atoms = atoms_index_tuple[0]

        res_seq_index += 1
        line_index = atoms_index_tuple[1]

        #add residue to chain structure
        chain.residues.append(res)
        chain.sequence = chain.sequence+res.res_char

    return chain


    



def get_res_atoms(index: int, lines: list[str]) -> tuple[list[Atom],int]:
    atoms=[]

    first_line = lines[index]
    res_id = get_residue_pdb_index(first_line)
    chain_id = get_chain_id(first_line)
    insertion_code = get_insertion_code(first_line)

    while(index < len(lines) and get_record_type(lines[index]) == "ATOM"):
        line=lines[index]

        if(get_residue_pdb_index(line) != res_id
           or get_chain_id(line) != chain_id
           or get_insertion_code(line) != insertion_code): break
              
        #get atom info        
        atom = Atom()
        atom.pdb_index= get_atom_index(line)
        atom.name=get_atom_name(line)
        atom.element=get_element_name(line)
        atom.coords = Point(
              get_x_coord(line),
              get_y_coord(line),
              get_z_coord(line)
        )
        #append to list
        atoms.append(atom)

        index+=1

    return atoms,index

#skip pdb header lines and return index of first ATOM line
def get_first_atom_line_index(lines: list[str])-> int:
    for i in range(0,len(lines)):
        if(lines[i].startswith("ATOM")): return i
    
    return -1


def get_record_type(line: str) -> str:
        return line[0:6].strip()

def get_atom_index(line: str) -> int:
        return int(line[6:11].strip())

def get_atom_name(line: str) -> str:
        return line[12:16].strip()

def get_residue_name(line: str) -> str:
        return line[17:20].strip()

def get_chain_id(line: str) -> str:
        return line[21].strip()

def get_residue_pdb_index(line: str) -> int:
        return int(line[22:26].strip())

def get_insertion_code(line: str) -> str:
        return line[26].strip()          

def get_x_coord(line: str) -> float:
        return float(line[30:38].strip())

def get_y_coord(line: str) -> float:
        return float(line[38:46].strip())

def get_z_coord(line: str) -> float:
        return float(line[46:54].strip())

def get_element_name(line: str) -> str:
      return line[76:78].strip()