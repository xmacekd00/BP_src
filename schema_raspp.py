import sys
from decimal import Decimal

import schemarecomb as sr
from schemarecomb.optimizers import RASPP


def get_split_indices(parents: sr.ParentSequences, number_of_cuts: int)-> list[int]:
    raspp = RASPP(parents,number_of_cuts)

    alignment_length = len(parents.alignment)

    #potential space for optimalization
    min_part_len = 20
    max_part_len = alignment_length - number_of_cuts * min_part_len

    starting_point_node = raspp.columns[0][0]

    # node -> (lowest energy path cost, path to node)    
    current_states = {starting_point_node: (Decimal(0), [])}

    #loop cuts
    for n in range(number_of_cuts):

        next_states = {}

        for node, (energy, path) in current_states.items():

            #loop through nodes which are the potential position for next cut
            for edge in node.out_edges:
                next_node = edge.out_node

                part_length = next_node.position - node.position

                #part length too small/big
                if(not (min_part_len <= part_length <= max_part_len)): continue
                
                new_energy = energy + edge.e_diff
                extended_path = path + [next_node.position]

                old_state = next_states.get(next_node)

                #decide if node has to be updated with new potential shotest path
                if (old_state is None or new_energy < old_state[0]):
                    next_states[next_node] = (new_energy,extended_path)

        current_states = next_states

    
    candidates = []

    #finalize composition of candidates
    for node, (energy, path) in current_states.items():
        last_part_length = alignment_length - node.position

        if (min_part_len <= last_part_length <= max_part_len): candidates.append((energy, path))

    if(not candidates): raise ValueError("Error: RASPP failed")

    return get_shortest_path(candidates)

def get_shortest_path(candidates: list[tuple])->list[int]:
    lowest_energy = None
    shortest_path = None

    for candidate in candidates:
        energy = candidate[0]
        path = candidate[1]

        if(lowest_energy is None or energy < lowest_energy):
            lowest_energy = energy
            shortest_path = path

    return shortest_path

def main():
    parents_fasta_path = sys.argv[1]
    pdb_path = sys.argv[2]
    number_of_cuts = int(sys.argv[3])

    pdb = sr.PDBStructure.from_pdb_file(pdb_path,chain="A")

    parents = sr.ParentSequences.from_fasta(parents_fasta_path,pdb_structure=pdb,prealigned=True)

    split_indices = get_split_indices(parents,number_of_cuts)

    print(",".join(str(index) for index in split_indices))

    return 0

if __name__ == "__main__":
    sys.exit(main())