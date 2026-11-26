import random

#create constants for the GC content range and sequence length
GC_CONTENT_MIN = 0.3
GC_CONTENT_MAX = 0.7
SEQUENCE_LENGTH = 1000

def generate_dna_sequence(length, gc_content):
    """Generate a random DNA sequence with specified length and GC content."""
    sequence = []
    for _ in range(length):
        if random.random() < gc_content:
            sequence.append(random.choice("GC"))
        else:
            sequence.append(random.choice("AT"))
    return "".join(sequence)

def wrap_sequence(sequence, width=80):
    """Wrap a sequence at specified width."""
    return "\n".join(sequence[i:i+width] for i in range(0, len(sequence), width))

# Generate the FASTA file
# create the fasta file in the revised_analyzer directory
with open("../revised_analyzer/sample.fasta", "w") as f:
    for i in range(10):
        seq_id = f"seq_{i+1}"
        description = f"Sequence {i+1} generated for testing purposes"
        gc_content = random.uniform(GC_CONTENT_MIN, GC_CONTENT_MAX)
        sequence = generate_dna_sequence(SEQUENCE_LENGTH, gc_content)
        wrapped_sequence = wrap_sequence(sequence)

        f.write(f">{seq_id} | {description}\n")
        f.write(f"{wrapped_sequence}\n")
