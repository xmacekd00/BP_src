import torch
import esm
import sys

model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
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
    
    return torch.mean(target_log_probs).item()

def evaluate_individual(sequence):
    predicted_tm = 1
    
    plausibility = get_plausibility_score(sequence)
    
    fitness = predicted_tm + (10.0 * plausibility)
    
    return fitness

if __name__ == "__main__":
    sequence = sys.argv[1]
    score = evaluate_individual(sequence)
    print(score)

