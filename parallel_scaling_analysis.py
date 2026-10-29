import subprocess

files = ["small.fasta", "medium.fasta", "large.fasta", "xlarge.fasta"]
pattern = "ATCG"
num_processes = 4

print("Starting automated benchmark...\n")
print(f"{'File Name':<15} | {'Parallel Time (s)':<18} | {'Serial Time (s)':<16} | {'Speedup':<10}")
print("-" * 65)

for file in files:
    input_data = f"{file}\n{pattern}\n"
    
    cmd = f"mpiexec -n {num_processes} .\dna_search_mpi.exe"
    result = subprocess.run(cmd, input=input_data, text=True, capture_output=True, shell=True)
    
   
    output = result.stdout
    par_time, ser_time, speedup = "N/A", "N/A", "N/A"
    
    for line in output.split('\n'):
        if "Parallel Time (slowest process):" in line:
            par_time = line.split(":")[-1].strip().split()[0]
        elif "Serial Time:" in line:
            ser_time = line.split(":")[-1].strip().split()[0]
        elif "Speedup:" in line:
            speedup = line.split(":")[-1].strip()
            
    print(f"{file:<15} | {par_time:<18} | {ser_time:<16} | {speedup:<10}")

print("\nBenchmark completed successfully!")