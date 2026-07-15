import sys
import torch
import esm


def score_sequence(sequence: str) -> float:
    model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    model.eval()

    batch_converter = alphabet.get_batch_converter()

    _, _, tokens = batch_converter([
        ("protein", sequence)
    ])

    position_scores = []

    with torch.no_grad():
        for sequence_index in range(len(sequence)):
            token_index = sequence_index + 1  # první token je BOS

            masked_tokens = tokens.clone()
            original_token = tokens[0, token_index].item()

            masked_tokens[0, token_index] = alphabet.mask_idx

            output = model(masked_tokens)
            logits = output["logits"][0, token_index]

            log_probabilities = torch.log_softmax(logits, dim=-1)

            position_scores.append(
                log_probabilities[original_token].item()
            )

    return sum(position_scores) / len(position_scores)


if __name__ == "__main__":
    sequence = sys.argv[1].strip().upper()
    score = score_sequence(sequence)

    print(score)