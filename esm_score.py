import torch
import esm
import sys

#for better results could be used differnt models like t33_650M
#this model is used for testing and debugging, because it is significantly faster
model, alphabet = esm.pretrained.esm2_t30_150M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()

def get_plausibility_score(sequence):
    data = [("protein", sequence)]
    _, _, tokens = batch_converter(data)
    
    with torch.no_grad():
        results = model(tokens)
        logits = results["logits"]
    
    toks = tokens[0, 1:-1]
    log_probs = torch.log_softmax(logits[0, 1:-1], dim=-1)
    target_log_probs = log_probs[torch.arange(len(toks)), toks]
    score = torch.mean(target_log_probs).item()
    return score

def evaluate_individual(sequence):
    predicted_tm = 1
    
    plausibility = get_plausibility_score(sequence)
    
    fitness = predicted_tm + (10.0 * plausibility)
    
    return fitness

if __name__ == "__main__":
    sequences_path = sys.argv[1]

    with open(sequences_path, "r") as f:
        sequences = [line.strip() for line in f if line.strip()
        ]

    for sequence in sequences:
        score = evaluate_individual(sequence)
        #"score:" is there to mark lines with a score
        #on stdout could be written different kind of stuff like downloading the model
        print(f"score:{score}")

