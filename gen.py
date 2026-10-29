import random

def generate_fasta(filename, size_in_mb):
    bases = ['A', 'C', 'G', 'T', 'N']
    target_size = size_in_mb * 1024 * 1024
    
    print(f"Generating file: {filename} (Size: {size_in_mb} MB)...")
    with open(filename, 'w') as f:
        f.write(f">Test_Sequence_{filename}\n")
        written = 0
        while written < target_size:
            line = "".join(random.choices(bases, k=60))
            f.write(line + "\n")
            written += 60
    print(f"Finished {filename} successfully! \n")

files_to_generate = {
    "small.fasta": 2,
    "medium.fasta": 20,
    "large.fasta": 60,
    "xlarge.fasta": 120
}

for filename, size in files_to_generate.items():
    generate_fasta(filename, size)

print("All requested files have been generated and are ready for testing!")