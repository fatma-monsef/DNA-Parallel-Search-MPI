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

        char base = (char)toupper(ch);
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

/* Validate the user-provided DNA pattern. */
int is_valid_pattern(const char* pattern, long pattern_len) {
    long i;
    if (pattern_len <= 0) return 0;
    for (i = 0; i < pattern_len; i++) {
        char base = (char)toupper((unsigned char)pattern[i]);
        if (base != 'A' && base != 'C' && base != 'G' &&
            base != 'T' && base != 'N') {
            return 0;
        }
    }
    return 1;
}

/* Serial reference search for correctness and speedup. */
long serial_count_matches(const char* dna, long dna_len,
                          const char* pattern, long pattern_len) {
    long count = 0;
    long i;
    if (pattern_len > dna_len) return 0;
    for (i = 0; i <= dna_len - pattern_len; i++) {
        long j;
        int match = 1;
        for (j = 0; j < pattern_len; j++) {
            if (dna[i + j] != pattern[j]) {
                match = 0;
                break;
            }
        }
        if (match) count++;
    }
    return count;
}

/* Read a full line from the user (supports spaces) and strip \n / \r.
   Also strips surrounding quotes (useful when a path is copied from Windows). */
void read_line(char *buffer, int size) {
    if (fgets(buffer, size, stdin)) {
        size_t len = strlen(buffer);
        while (len > 0 && (buffer[len - 1] == '\n' || buffer[len - 1] == '\r')) {
            buffer[--len] = '\0';
        }
        if (len >= 2 && buffer[0] == '"' && buffer[len - 1] == '"') {
            memmove(buffer, buffer + 1, len - 2);
            buffer[len - 2] = '\0';
        }
    } else {
        buffer[0] = '\0';
    }
}

int main(int argc, char** argv) {
    int rank, size;

    MPI_Init(&argc, &argv);
    /* Show input prompts immediately in the MPI console. */
    setvbuf(stdout, NULL, _IONBF, 0);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    long dna_len = 0;
    char *dna = NULL;
    char pattern[MAX_PATTERN_LEN];
    long pattern_len = 0;

    /* Shared validation and performance variables. */
    int valid_input = 1;
    long local_search_count = 0;
    double local_elapsed = 0.0;
    double max_parallel_time = 0.0;
    double serial_time = 0.0;
    long serial_matches = 0;

    memset(pattern, 0, MAX_PATTERN_LEN);

    /* ========================================================================= */
    /* TEAM MEMBER 2: Master Root & Broadcast Module                             */
    /* Process 0 handles user input, reads the pattern, and broadcasts both the  */
    /* DNA sequence properties and pattern data to all worker processes.         */
    /* ========================================================================= */
    if (rank == 0) {
        char filename[512];
        char pattern_input[512];

        /* ===== 1) Ask the user for the FASTA file path (with retry) ===== */
        while (1) {
            printf("Please enter the path of the FASTA file: ");
            fflush(stdout);
            read_line(filename, sizeof(filename));

            if (strlen(filename) == 0) {
                printf("Error: File path cannot be empty. Try again.\n");
                continue;
            }

            printf("Process 0: Reading DNA file...\n");
            dna = read_and_clean_fasta(filename, &dna_len);

            if (dna == NULL) {
                printf("Please check the path and try again.\n\n");
                continue;
            }
            if (dna_len == 0) {
                printf("Error: No valid DNA bases found in the file. Try another file.\n\n");
                free(dna);
                dna = NULL;
                continue;
            }
            break;
        }
        printf("Process 0: DNA loaded successfully. Total length = %ld\n", dna_len);

        /* ===== 2) Ask the user for the pattern (with retry) ===== */
        while (1) {
            printf("Please enter the DNA pattern to search for (A, C, G, T, N): ");
            fflush(stdout);
            read_line(pattern_input, sizeof(pattern_input));
            pattern_len = (long)strlen(pattern_input);

            if (pattern_len == 0) {
                printf("Error: Pattern cannot be empty. Try again.\n");
                continue;
            }
            if (pattern_len >= MAX_PATTERN_LEN) {
                printf("Error: Pattern is too long (max %d characters). Try again.\n",
                       MAX_PATTERN_LEN - 1);
                continue;
            }

            for (long i = 0; i < pattern_len; i++) {
                pattern_input[i] = (char)toupper((unsigned char)pattern_input[i]);
            }

            if (!is_valid_pattern(pattern_input, pattern_len)) {
                printf("Error: Pattern must contain only A, C, G, T, or N. Try again.\n");
                continue;
            }
            if (pattern_len > dna_len) {
                printf("Error: Pattern is longer than the DNA sequence (%ld). Try again.\n",
                       dna_len);
                continue;
            }
            break;
        }

        strcpy(pattern, pattern_input);
        printf("Searching for pattern \"%s\" ...\n", pattern);
    }

    /* Tell every process whether the input is valid. */
    MPI_Bcast(&valid_input, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!valid_input) {
        free(dna);
        MPI_Finalize();
        return 1;
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
        sendcounts[i] = (int)(core_size + current_overlap);
        displs[i] = (int)current_disp;
        current_disp += core_size;
    }

    long my_core_start = rank * base + (rank < remainder ? rank : remainder);
    long my_core_size = base + (rank < remainder ? 1 : 0);
    long my_overlap = (rank == size - 1) ? 0 : overlap;
    long my_recv_size = my_core_size + my_overlap;

    char *local_dna = (char*)malloc((my_recv_size + 1) * sizeof(char));
    if (!local_dna) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Scatterv(dna, sendcounts, displs, MPI_CHAR,
                 local_dna, (int)my_recv_size, MPI_CHAR, 0, MPI_COMM_WORLD);
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

    /* Search only valid global start positions. This prevents reading
       beyond local_dna near the end of the DNA sequence. */
    if (my_core_start < dna_len - pattern_len + 1) {
        long remaining_valid_starts = dna_len - pattern_len + 1 - my_core_start;
        local_search_count = (remaining_valid_starts < my_core_size)
                           ? remaining_valid_starts : my_core_size;
    }

    for (long i = 0; i < local_search_count; i++) {
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

    /* Use the slowest process time as the parallel time. */
    local_elapsed = end_time - start_time;
    MPI_Reduce(&local_elapsed, &max_parallel_time, 1, MPI_DOUBLE,
               MPI_MAX, 0, MPI_COMM_WORLD);

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

        /* Compare with a serial reference and print speedup/efficiency. */
        {
            double serial_start = MPI_Wtime();
            serial_matches = serial_count_matches(dna, dna_len, pattern, pattern_len);
            serial_time = MPI_Wtime() - serial_start;

            printf("Serial Matches: %ld\n", serial_matches);
            printf("Results Match: %s\n",
                   (serial_matches == total_matches) ? "YES" : "NO");
            printf("Parallel Time (slowest process): %.9f seconds\n",
                   max_parallel_time);
            printf("Serial Time: %.9f seconds\n", serial_time);
            if (max_parallel_time > 0.0) {
                printf("Speedup: %.4f\n", serial_time / max_parallel_time);
                printf("Efficiency: %.4f\n",
                       serial_time / max_parallel_time / size);
            }
        }

        free(dna);
    }

    free(sendcounts);
    free(displs);
    free(local_dna);

    MPI_Finalize();
    return 0;
}
