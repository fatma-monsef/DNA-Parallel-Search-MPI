#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <mpi.h>

#define MAX_PATTERN_LEN 100

/* ========================================================================= */
/* TEAM MEMBER 1: Data Cleaning & Loading Module                             */
/* Responsible for reading the FASTA file, stripping headers/whitespace,     */
/* converting bases to uppercase, and building the raw DNA sequence string.  */
/* ========================================================================= */
char* read_and_clean_fasta(const char* filename, long* out_dna_len) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error: Could not open file %s\n", filename);
        return NULL;
    }

    long capacity = 1000;
    long dna_len = 0;
    char *dna = (char*)malloc(capacity * sizeof(char));
    if (!dna) {
        fclose(file);
        return NULL;
    }

    int ch;
    int in_header = 0;

    while ((ch = fgetc(file)) != EOF) {
        if (ch == '>') {
            in_header = 1;
            continue;
        }
        if (ch == '\n') {
            in_header = 0;
            continue;
        }
        if (in_header) continue;
        if (isspace(ch)) continue;

        char base = toupper((char)ch);
        if (base == 'A' || base == 'C' || base == 'G' || base == 'T' || base == 'N') {
            if (dna_len >= capacity - 1) {
                capacity *= 2;
                char *temp = (char*)realloc(dna, capacity * sizeof(char));
                if (!temp) {
                    free(dna);
                    fclose(file);
                    return NULL;
                }
                dna = temp;
            }
            dna[dna_len++] = base;
        }
    }

    dna[dna_len] = '\0';
    *out_dna_len = dna_len;
    fclose(file);
    return dna;
}

int main(int argc, char** argv) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    long dna_len = 0;
    char *dna = NULL;
    char pattern[MAX_PATTERN_LEN];
    long pattern_len = 0;

    /* ========================================================================= */
    /* TEAM MEMBER 2: Master Root & Broadcast Module                             */
    /* Process 0 handles user input, reads the pattern, and broadcasts both the  */
    /* DNA sequence properties and pattern data to all worker processes.         */
    /* ========================================================================= */
    if (rank == 0) {
        const char *filename = "dna.fasta";
        printf("Process 0: Reading DNA file...\n");
        dna = read_and_clean_fasta(filename, &dna_len);

        if (!dna) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        printf("Process 0: DNA loaded successfully. Total length = %ld\n", dna_len);

        printf("Enter the pattern to search for: ");
        scanf("%s", pattern);
        pattern_len = strlen(pattern);

        for (int i = 0; i < pattern_len; i++) {
            pattern[i] = toupper(pattern[i]);
        }
    }

    MPI_Bcast(&dna_len, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&pattern_len, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(pattern, MAX_PATTERN_LEN, MPI_CHAR, 0, MPI_COMM_WORLD);

    /* ========================================================================= */
    /* TEAM MEMBER 3: Workload Distribution & Overlap Management Module          */
    /* Calculates balanced chunks per process, manages boundary overlaps (L-1),  */
    /* and sets up MPI_Scatterv parameters for data distribution.                */
    /* ========================================================================= */
    long base = dna_len / size;
    long remainder = dna_len % size;
    long overlap = pattern_len - 1;

    int *sendcounts = (int*)malloc(size * sizeof(int));
    int *displs = (int*)malloc(size * sizeof(int));

    long current_disp = 0;
    for (int i = 0; i < size; i++) {
        long core_size = base + (i < remainder ? 1 : 0);
        long current_overlap = (i == size - 1) ? 0 : overlap;
        sendcounts[i] = core_size + current_overlap;
        displs[i] = current_disp;
        current_disp += core_size;
    }

    long my_core_size = base + (rank < remainder ? 1 : 0);
    long my_overlap = (rank == size - 1) ? 0 : overlap;
    long my_recv_size = my_core_size + my_overlap;

    char *local_dna = (char*)malloc((my_recv_size + 1) * sizeof(char));
    if (!local_dna) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Scatterv(dna, sendcounts, displs, MPI_CHAR,
                 local_dna, my_recv_size, MPI_CHAR, 0, MPI_COMM_WORLD);
    local_dna[my_recv_size] = '\0';

    /* Start timing *after* setup and input completion for accurate benchmarking */
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    /* ========================================================================= */
    /* TEAM MEMBER 4: Parallel Search Algorithm Module                           */
    /* Each process searches its assigned core window for exact pattern matches  */
    /* while utilizing the overlap region to prevent boundary missing matches.   */
    /* ========================================================================= */
    long local_count = 0;
    for (long i = 0; i < my_core_size; i++) {
        int match = 1;
        for (long j = 0; j < pattern_len; j++) {
            if (local_dna[i + j] != pattern[j]) {
                match = 0;
                break;
            }
        }
        if (match) {
            local_count++;
        }
    }

    long total_matches = 0;
    MPI_Reduce(&local_count, &total_matches, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    double end_time = MPI_Wtime();

    /* ========================================================================= */
    /* TEAM MEMBER 5: Results Aggregation & Performance Testing Module           */
    /* Root process gathers final match counts, computes total execution time,   */
    /* and prints comprehensive performance reports and benchmarking metrics.    */
    /* ========================================================================= */
    if (rank == 0) {
        printf("\n========================================\n");
        printf("         MPI Parallel Execution Results \n");
        printf("========================================\n");
        printf("Total Processes: %d\n", size);
        printf("Total DNA Length: %ld\n", dna_len);
        printf("Pattern: %s\n", pattern);
        printf("Total Matches Found: %ld\n", total_matches);
        printf("Total Execution Time: %f seconds\n", end_time - start_time);
        printf("========================================\n");

        free(dna);
    }

    free(sendcounts);
    free(displs);
    free(local_dna);

    MPI_Finalize();
    return 0;
}
