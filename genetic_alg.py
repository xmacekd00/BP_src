#celkove splitting...asi rozjet pres venv nejaky stary RASP/SCHEMA

from data_structures import ProtChain
from data_structures import PdbProtChain
from pathlib import Path
from schema import get_split_indices
from schema import get_chains_from_fasta_file
from pdb_structure_parser import parse_pdb_file
from Bio import Phylo
import sys
import random
import subprocess
import shutil
import tempfile
import os
import csv
import itertools

split_indices = []
population_size = 0
sequences_to_next_generation_count = 0
AGGREPROT_BIN = "/home/david-macek/Documents/VUT_FIT/BP/solution/.venv-aggreprot/bin/aggreprot-predictor"
AMINO_ACID_ORDER = "ARNDCQEGHILKMFPSTWYV"


BASE_DIR = Path(__file__).resolve().parent.parent
ESM_PYTHON = BASE_DIR / ".venv-esm/bin/python"
ESM_SCRIPT = Path(__file__).resolve().parent / "esm_score.py"

SCHEMA_PYTHON = BASE_DIR / ".venv-schema/bin/python"
SCHEMA_SCRIPT = Path(__file__).resolve().parent / "schema_raspp.py"

#functions returns a list of names of all direct ancestor nodes of query
def get_query_ancestors(tree_path: str)-> list[str] | None:
    tree = Phylo.read(tree_path,"newick")
    query = tree.find_any(name="query")

    if(query is None):
        print(f"Error: no leaf node named query in {tree_path}\n")
        return None
    
    path = tree.get_path(query)

    #add root and remove qeury node
    ancestor_nodes = [tree.root] + path[:-1]

    return [str(int(node.confidence)) for node in ancestor_nodes]

#function reads given full_msa input file and returns query ProtChain representation
def get_query_representation(base_path: str)-> ProtChain:
    full_msa_data = load_sequences_from_full_msa_file(base_path)

    #get query msa
    query_msa = ""
    for i in range(len(full_msa_data)):
        if(full_msa_data[i]=="query"):
            query_msa=full_msa_data[i+1]

    query_seq = query_msa.replace("-","")

    query: ProtChain = ProtChain(
        id="query",
        sequence=query_seq,
        residues=list(query_seq),
        aligned_sequence=query_msa,
        score=0.0
        )
    
    return query

# function returns list of protein chains parsed out of .dat files
def init_population(asr_folder_base_path: str) -> list[ProtChain | PdbProtChain]:
    #get numbers of ancestor nodes
    ancestor_nodes = get_query_ancestors(asr_folder_base_path + "asr/ancestral_tree.tre")
    if ancestor_nodes is None: return None

    #get query representation
    query: ProtChain = get_chains_from_fasta_file(asr_folder_base_path + "query.fasta")[0]

    dat_folder_path = asr_folder_base_path + "asr/lazarus_tree_nodes/tree1/"

    #make list of .dat file names from node numbers
    ancestor_dat_files_names = [dat_folder_path + "node" + number + ".dat" for number in ancestor_nodes]

    chains = []
    chains.append(query)
    index = 0

    #extract from all .dat files protein sequences and append them to list
    for file in ancestor_dat_files_names:
        aligned_sequence = extract_sequence_from_dat_file(file)
        sequence = aligned_sequence.replace("-","")

        chain = ProtChain(
            id="node"+str(index),
            sequence=sequence,
            residues=list(sequence),
            aligned_sequence=aligned_sequence)
        index+=1

        chains.append(chain)
        
    return chains

#function returns a protein sequence consisted of most probable residue on every position
def extract_sequence_from_dat_file(path: Path) -> str:
    with open(path, "r") as f:
        lines = f.readlines()

    sequence = ""
    
    for line in lines:
        parts = line.split()
        #invalid row
        if not parts or len(parts)<2: continue
        #most probable residue
        residue = parts[1]

        sequence += residue
        
    return sequence
    
#function computes fitness score of each individual in population and sets its .score
def eval_population(population: list[ProtChain])-> int | None:
    
    for chain in population:
        #already scored
        if(chain.score != 0.0): continue

        score = get_esm_score(chain.sequence)
        chain.score = score
    return 0
        
#function runns esm scoring script and returns score for given sequence
def get_esm_score(seq: str) -> float:
    result = subprocess.run(
        [
            str(ESM_PYTHON),
            str(ESM_SCRIPT),
            seq
        ],
        capture_output=True,
        text=True,
        check=True
    )

    return float(result.stdout.strip())

#function returns list of pairs of sequences ready for crossover
#random weighted choice based on a fitnes score
def get_crossover_pairs(population: list[ProtChain])-> list[list[ProtChain]]:
    
    chains = []
    weights = []
    #put weigths and sequences into separate lists
    for i in range(len(population)):
        weights.append(population[i].score)
        chains.append(population[i])

    pairs = []
    pairs_to_generate = population_size - sequences_to_next_generation_count
    #generate pairs
    for i in range(pairs_to_generate):
        parent1 = random.choices(chains, weights=weights, k=1)[0]
        parent2 = random.choices(chains, weights=weights, k=1)[0]

        #prevent choice of two same sequences
        while parent2 == parent1:
            parent2 = random.choices(chains, weights=weights, k=1)[0]

        pairs.append([parent1, parent2])


    return pairs

#function makes a new population for the next iteration of genetic algorithm
def do_crossover(
        population: list[ProtChain],
        generation_number: int
        )->list[ProtChain]:
    
    crossover_pairs = get_crossover_pairs(population)
    children_index=0

    new_generation = []

    for i in range(len(crossover_pairs)):
        parent1 = crossover_pairs[i][0]
        parent2 = crossover_pairs[i][1]
        
        #split sequences
        par1_splitted = get_splitted_sequences(parent1,split_indices)
        par2_splitted = get_splitted_sequences(parent2,split_indices)
        
        #calculate which parts are gonna be used from par1 and which from par2
        new_aligned_seq=""
        #iterate through parts and build the sequence of par1 and par2
        for j in range(len(par1_splitted)):
            par1_score = sum(par1_splitted[j].aggreprot_scores)
            par2_score = sum(par2_splitted[j].aggreprot_scores)
            #score transformed to range <0,1>
            weight1 = 1 / (par1_score + 1e-6)
            weight2 = 1 / (par2_score + 1e-6)
            #random weighted choice of part from par1 or par2
            choice = random.choices(
                [par1_splitted[j],par2_splitted[j]],
                weights=[weight1,weight2],
                k=1)[0]
            
            new_aligned_seq+=choice.aligned_sequence

        new_seq = new_aligned_seq.replace("-","")

        new_chain = ProtChain(
            id=f"gen{generation_number}_child{children_index}",
            sequence=new_seq,
            residues=list(new_seq),
            aligned_sequence=new_aligned_seq
            )
        new_generation.append(new_chain)
        children_index +=1
    return new_generation


#function splits sequence by given indices and returns list of ProtChain structures
#One ProtChain strcture is one part of parent sequence
def get_splitted_sequences(chain: ProtChain, split_indices: list[int])-> list[ProtChain]:
    parts: list[ProtChain] = []
    indices = [0] + split_indices + [len(chain.aligned_sequence)]


    for i in range(len(indices)-1):
        #start and end index of current part in aligned indexing
        start = indices[i]
        end = indices[i+1]
        
        #get aligned and raw part
        aligned_part_sequence = chain.aligned_sequence[start:end]
        part_sequence = aligned_part_sequence.replace("-","")

        #transform aligned sequence start index to raw sequence start index
        raw_start = sum(
            residue != "-"
            for residue in chain.aligned_sequence[:start])
        
        raw_end = sum(
            residue != "-"
            for residue in chain.aligned_sequence[:end]
        )

        part_score = chain.aggreprot_scores[raw_start:raw_end]

        part = ProtChain(
            id=f"{chain.id}_part{i}",
            sequence=part_sequence,
            residues=list(part_sequence),
            aggreprot_scores=part_score,
            aligned_sequence=aligned_part_sequence)
        parts.append(part)

    return parts

#function appends 10 best proteins to new population based on their score
def add_best_performers_to_new_generation(population: list[ProtChain], new_population: list[ProtChain]):
    chains = []
    scores = []
    #put scores and sequences into separate lists
    for i in range(len(population)):
        scores.append(population[i].score)
        chains.append(population[i])


    for i in range(sequences_to_next_generation_count):
        #find index of protein with best score
        best_performer_index = scores.index(max(scores))
        new_population.append(chains[best_performer_index])
        #remove already appended best protein 
        scores.pop(best_performer_index)
        chains.pop(best_performer_index)

#function runs aggreprot tool on each chain and assigns result values to each chain.aggreprot_score
def load_aggreprot_scores(population: list[ProtChain]):
    tmp_dir = tempfile.mkdtemp(dir=".")

    try:
        #create required temporary files and folders
        tmp_path = Path(tmp_dir)
        fasta_path = tmp_path / "population.fasta"
        aggreprot_out_dir = tmp_path / "aggreprot_out"

        #write all sequences to fasta file
        with open(fasta_path,"w") as f:
            for i in range(len(population)):
                f.write(f">{population[i].id}\n")
                f.write(f"{population[i].sequence}\n")

        run_aggreprot(aggreprot_out_dir,fasta_path)

        aggreprot_result_path = aggreprot_out_dir / "batch-profile-final-mean.txt"

        with open(aggreprot_result_path,"r") as f:
            lines = f.readlines()

        process_aggreprot_output(lines, population)

    finally:
        #remove the temporary folder with all content
        if(os.path.exists(tmp_dir)):
            shutil.rmtree(tmp_dir)
    
    
#function runs aggreprot tool with given arguments 
def run_aggreprot(out_dir:Path, fasta_path: Path):

    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = "-1"
    env["TF_CPP_MIN_LOG_LEVEL"] = "3"

    cmd = [
        AGGREPROT_BIN,
        "predict-sequential-batch",
        "--out-dir",
        str(out_dir),
        str(fasta_path)
    ]

    subprocess.run(cmd, check=True, env=env)

#function parses output of aggreprot tool and assigns results to chains
def process_aggreprot_output(lines: list[str],population: list[ProtChain]):
    chain_id_map = {}
    for chain in population:
        chain_id_map[chain.id]=chain

    for i in range(len(lines)):
        parts = lines[i].split()

        if(len(parts)<3): continue

        line_id = parts[0]
        scores = [float(x) for x in parts[2:]]

        if(line_id in chain_id_map):
            chain_id_map[line_id].aggreprot_scores = scores

def load_sequences_from_full_msa_file(base_path: str)->list[str]:
    path = base_path + "asr/full_msa.fasta"
    #read full_msa.fasta
    with open(path,"r") as f:
        lines  = f.readlines()

    #remove sequence id lines and blank lines
    sequences = []
    current_seq = ""

    #append first sequences id
    sequences.append(lines[0][1:].strip())

    for i in range(1,len(lines)):
        #new sequence begins, so append the previous one
        if(lines[i].startswith(">")):
            sequences.append(current_seq.strip())
            current_seq = ""
            #append next sequence id
            sequences.append(lines[i][1:].strip())
        else:
            current_seq += lines[i].strip()
    #append last sequnce
    sequences.append(current_seq)

    return sequences

#function computes conservation rate of each column from full_msa
def compute_column_conservation(folder_base_path: str)-> list[float] | None:
    full_msa_data = load_sequences_from_full_msa_file(folder_base_path)

    sequences = [full_msa_data[i] for i in range(len(full_msa_data)) if i%2==1]
    

    if(len(sequences)==0):
        print("Error: full_msa.fasta is empty or bug in compute_column_conservation\n")
        return None
    
    col_conservation_score = []

    for i in range(len(sequences[0])):
        residue_count_dict = {}

        for j in range(len(sequences)):
            residue = sequences[j][i]
            
            #skip gaps
            if(residue=="-"):continue

            if(residue not in residue_count_dict):
                residue_count_dict[residue] = 1
            else:
                residue_count_dict[residue] += 1

        #all items in this column are gaps
        if(not residue_count_dict):
            col_conservation_score.append(0.0)
            continue
        
        #compute which residues appears the most
        max_residue = max(residue_count_dict,key=residue_count_dict.get)
        max_res_count = residue_count_dict[max_residue]

        total_res_count = sum(residue_count_dict.values())

        col_conservation_score.append(max_res_count/total_res_count)

    return col_conservation_score

#function applies mutation on each individual of a population
def do_mutation(
        population: list[ProtChain],
        conservation_scores: list[float],
        posterior_prob: list[list[float]],
        gap_prob: list[float]):
    
    #square the value to lower the probability of stable residues and to rise the probability of unstable residues
    mutation_weights = [(1.0 - score) ** 2 for score in conservation_scores]


    for chain in population:
        #get list of indices, where could be mutation applied
        possible_mutation_indices = [
            i
            for i, residue in enumerate(chain.aligned_sequence)
            if residue != "-"
            and gap_prob[i] < 0.5
            and mutation_weights[i] > 0
        ]

        #no residue can not be mutated
        if(not possible_mutation_indices): continue
        
        #choose only those weights that are 
        possible_mut_weights = [mutation_weights[i] for i in possible_mutation_indices]

        mutation_index = random.choices(
            possible_mutation_indices,
            weights=possible_mut_weights,
            k=1
            )[0]
        
        old_residue = chain.aligned_sequence[mutation_index]

        possible_aa = []
        aa_weights = []

        for index, residue in enumerate(AMINO_ACID_ORDER):
            #t choosing the same residue to increase the diversity of population
            if(residue == old_residue): continue

            possible_aa.append(residue)
            aa_weights.append(posterior_prob[mutation_index][index])

        #no mutation
        if (sum(aa_weights) == 0): continue

        new_residue =  random.choices(
            possible_aa,
            weights=aa_weights,
            k=1
        )[0]
        
        #change old residue to new one
        aligned_list = list(chain.aligned_sequence)
        aligned_list[mutation_index] = new_residue

        chain.aligned_sequence = "".join(aligned_list)
        chain.sequence = chain.aligned_sequence.replace("-","")
        chain.residues = list(chain.sequence)
        chain.score = 0.0
        chain.aggreprot_scores = []

#function makes from csv file posterior and gap probabilities lists
def load_posterior_probabilities(folder_base_path: str
                                 ) -> tuple[list[list[float]],list[float]]:
    csv_path = folder_base_path + "asr/ancestral_profile.csv"
    rows = []

    with open(csv_path,"r") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    max_index = max(int(row["position"]) for row in rows)

    #init helper lists with zeros
    aa_sums = [[0.0] * len(AMINO_ACID_ORDER) for _ in range(max_index)]
    gap_sums = [0.0] * (max_index)
    counts = [0] * (max_index)

    for row in rows:
        index = int(row["position"])-1

        #parse gap probability
        gap_prob = row["-"]
        if(gap_prob != "-"):
            gap_sums[index] += float(gap_prob)

        #parse amino acids probabilities
        for acid_index, acid in enumerate(AMINO_ACID_ORDER):
            prob = row[acid]

            if(prob != "-"):
                aa_sums[index][acid_index] += float(prob)

        counts[index] +=1

    #init output lists with zeros
    posterior_prob = [[0.0] for _ in range(max_index)]
    gap_prob = [0.0] * max_index

    for index in range(max_index):
        #compute the divisor for normalization
        total_weigth = sum(aa_sums[index])

        #normalize the sums of probabilities
        if(total_weigth>0):
            posterior_prob[index] = [
                value / total_weigth
                for value in aa_sums[index]
            ]
        
        #average of probabilities
        if(counts[index] > 0):
            gap_prob[index] = gap_sums[index] / counts[index]

    return posterior_prob, gap_prob

#function takes each sequence divided into fragments and makes a non-fragmented seuquence 
def assamble_sequences(combinations, list_of_parts: list[list[str]])->list[str]:
    expanded_seqs = []
    #assamble sequences from parts
    for combination in combinations:
        new_seq = ""

        for fragmen_index, parent_index in enumerate(combination):
            new_seq += list_of_parts[parent_index][fragmen_index]
        
        expanded_seqs.append(new_seq)

    return expanded_seqs

#prepare a list of splitted sequences at given indices
def get_fragmented_population(population: list[ProtChain])-> list[list[str]]:
    list_of_parts: list[list[str]] = []
    
    #make a list of fragments from each sequence in population 
    for i in range(len(population)):
        seq = population[i].aligned_sequence
        boundaries = [0] + split_indices + [len(seq)]
        #append list of fragmented sequence
        list_of_parts.append([
            seq[boundaries[i]:boundaries[i+1]] 
            for i in range(len(boundaries)-1)
            ])
        
    return list_of_parts


#function randomly expands population to targeted size
def expand_population(
        population: list[ProtChain],
        iteration: int,
        target_size: int =100
        )->list[ProtChain]:
    
    population_len = len(population)
    fragment_count = len(split_indices)+1
    sequences_to_create = target_size - population_len
    total_combinations = population_len ** fragment_count

    #prepare a list of splitted sequences at given indices
    list_of_parts = get_fragmented_population(population)
    
    #take all combinations if number of combinations is lesser than targeted size
    if(total_combinations <= target_size):
        combinations = itertools.product(
            range(population_len),
            repeat=fragment_count
        )

    else:
        #unique values only
        combinations = set()

        while(len(combinations)<target_size):
            combination = tuple(
                random.randrange(population_len)
                for i in range(fragment_count)
            )

            combinations.add(combination)
        
        combinations = list(combinations)


    new_combinations = []
    #remove existing combinations in base population
    for combination in combinations:
        if(len(set(combination))>1):
            new_combinations.append(combination)
    
    expanded_seqs: list[str] = assamble_sequences(new_combinations,list_of_parts)

    #make ProtChains from sequences
    expanded_protchains: list[ProtChain] = [
        ProtChain(
            id="expanded_"+ str(iteration) +"_"+str(i),
            sequence=expanded_seqs[i].replace("-",""),
            residues=list(expanded_seqs[i].replace("-","")),
            aligned_sequence=expanded_seqs[i],
            score=0.0,
            )
        for i in range(len(expanded_seqs))]

    #append population base to expanded population
    expanded_population = []
    expanded_population.extend(population)
    #append new expanded chains
    index = 0
    while(
        len(expanded_population) < target_size 
        and index < len(expanded_protchains)
        ):
        expanded_population.append(expanded_protchains[index])
        index+=1

    return expanded_population

#function makes temporary fasta file with all population chains,runs raspp/schema and returns split indices
def get_schema_split_inidices(population: list[ProtChain],
                              pdb_path: str,
                              num_of_cuts: int
                              )->list[int] | None:
        #make tmp fasta file for schema        
        with tempfile.TemporaryDirectory(dir=".") as tmp_dir:
            fasta_path = Path(tmp_dir) / "schema_parents.fasta"

        #write all chains
        with open(fasta_path, "w") as f:
            for chain in population:
                f.write(f">{chain.id}\n")
                f.write(f"{chain.aligned_sequence}\n")
        #run raspp/schema
        result = subprocess.run(
            [
                str(SCHEMA_PYTHON),
                str(SCHEMA_SCRIPT),
                str(fasta_path),
                pdb_path,
                str(num_of_cuts)
            ],
            capture_output=True,
            text=True
        )

        if(result.returncode != 0):
            print(f"Error: raspp/schema - {result.stderr}\n")
            return None

        output = result.stdout.strip()

        if(not output):
            print("Error: schema did not work")
            return None

        split_indices = [
            int(value)
            for value in output.split(",")
        ]

        return split_indices

#argv[1] - path to ASR folder
def main():
    global split_indices
    global population_size
    global sequences_to_next_generation_count

    #after test should be = argv[0]
    asr_folder_base_path = "/home/david-macek/Documents/VUT_FIT/BP/materialy/ASR_data/p1j0p2/"

    #population init
    population = init_population(asr_folder_base_path)
    if(len(population)<2):
        print("Error: population too small - less than 2 proteins\n")
        return 1
    
    
    #top performing 10% of current generation will be promoted to next generation without a change
    population_size = len(population)
    sequences_to_next_generation_count = population_size // 10

    #pdb_chain = parse_pdb_file(sys.argv[1])
    pdb_chain = parse_pdb_file(asr_folder_base_path)
    if(pdb_chain is None): return None

    #get the indices of cuts
    split_indices = get_split_indices(
        pdb_chain,
        parent_chains = [population[0],population[1]]
        )

    #compute conservation scores for every index
    conservation_scores = compute_column_conservation(asr_folder_base_path)
    if(conservation_scores is None): return None
        
    #get posterior probabilities and gap probabilities lists
    posterior_prob, gap_prob = load_posterior_probabilities(asr_folder_base_path)

    generation_index = 0
    while generation_index != 100:
        #expand population
        expanded_population = expand_population(population,generation_index)

        #run aggreprot tool on population
        load_aggreprot_scores(expanded_population)

        #evaluation of popuplation stage
        eval_population(expanded_population)

        #get top 10% performing elite
        elite = []
        add_best_performers_to_new_generation(expanded_population,elite)

        #crossover
        new_generation = do_crossover(expanded_population,generation_index+1)
        
        #mutation stage - mutate only children
        do_mutation(new_generation,conservation_scores,posterior_prob,gap_prob)

        #add elite to new generation
        new_generation += elite

        #next iteration(generation)
        population = new_generation
        generation_index += 1
        

if __name__ == "__main__":
    sys.exit(main())

