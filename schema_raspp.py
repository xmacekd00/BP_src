import sys
import schemarecomb as sr

def main():
    fasta_path = sys.argv[1]
    pdb_path = sys.argv[2]
    number_of_cuts = int(sys.argv[3])

    pdb = sr.PDBStructure.from_pdb_file(pdb_path,chain="A")

    parents = sr.ParentSequences.from_fasta(fasta_path,pdb_structure=pdb_path,prealigned=True)

    libraries = sr.generate_libraries(parents,number_of_cuts)

    if(not libraries):
        print("Error: no schema libraries generated\n")
        return 1
    
    best_lib = max(libraries,key= lambda lib: lib.mutation_rate - lib.energy)

    alignment_length = len(parents.alignment)

    split_indices = [
        breakpoint.position for breakpoint in best_lib.breakpoints if 0 < breakpoint.position < alignment_length
    ]

    print(",".join(str(i) for i in split_indices))

    return 0

if __name__== "__main__":
    sys.exit(main())