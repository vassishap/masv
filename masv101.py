# ######################################################################################
#
# masv101.py
#
# Version: MASV v.1.0.1
#
# Description:
#   This script performs denoising of sequencing data from a FASTA file.
#   It identifies unique sequences, dereplicates them, and then classifies them
#   as either valid sequence variants or noise based on k-mer profile comparisons
#   and abundance ratios.
#
# Usage:
#   python masv101.py -i <input.fasta> -f <degree_of_freedom> -a <abundance_ratio> -s <save_spurious_as_variant>
#
# Arguments:
#   -i, --input_fasta: Path to the input FASTA file.
#   -f, --freedom: An integer representing the degree of freedom for noise filtering.
#   -a, --abundance_ratio: A float representing the minimum abundance ratio (parent/noise) required for filtering.
#   -s, --save_spurious: A boolean flag (True/False) to save sequences classified as 'SPURIOUS VARIANT' to the variants output file.
#
# ######################################################################################

import argparse, sys, time
from collections import Counter

# --- Argument Parsing ---
# Sets up the command-line interface for the user to provide input files and parameters.
parser = argparse.ArgumentParser(description="MASV v.1.0.1: A high-resolution denoising script for amplicon data.")
parser.add_argument('-i', dest='input_fasta', type=str, required=True,
                    help='Path to the input FASTA file.')
parser.add_argument('-f', dest='freedom', type=int, required=False, default=1,
                    help="Degree of freedom (integer). A multiplier for the denoising thresholds. "
                         "'-f 1' is the most stringent (default). '-f 2' is more lenient, etc.")
parser.add_argument('-a', dest='abundance_ratio', type=float, required=False, default=1.45,
                    help="Minimum abundance ratio (Parent Sequence Size / Noise Sequence Size) "
                         "required to classify a sequence as noise. Default is 1.45.")
parser.add_argument('-s', dest='save_spurious', type=lambda x: (x.lower() == 'true'), required=False, default=False,
                    help="Boolean flag (True/False). If True, singleton sequences classified as 'SPURIOUS VARIANT' "
                         "are saved to the 'variants.fa' file. If False, they are saved to 'noise.fa'. Default is False.")
args = parser.parse_args()


class Sequence:
    """
    Represents a biological sequence and its associated k-mer profile.

    This class stores the sequence information, calculates its k-mer profile,
    and provides a method to compare its profile with another sequence.
    The sequence is "doubled" (e.g., 'A' becomes 'AA') before k-mer generation
    to handle short sequences and minimize edge effects.
    """

    def __init__(self, title, sequence, k=2):
        """
        Initializes a Sequence object.

        Args:
            title (str): The FASTA header/title for the sequence.
                         The size (abundance) is parsed from this title.
            sequence (str): The DNA sequence string.
            k (int): The size of the k-mers to generate.
        """
        self.title = title  # Sequence identifier from FASTA header
        self.size = int(title.rsplit("=", 1)[-1][:-1]) # Abundance, parsed from title
        self.sequence = sequence.upper()  # The DNA sequence, converted to uppercase
        self.length = len(sequence)  # Length of the original sequence
        self.k = k  # k-mer size
        # Create a "doubled" sequence to enhance k-mer profile differences
        self.doubled_sequence = self.sequence.replace('A', 'AA').replace('C', 'CC').replace('T', 'TT').replace('G', 'GG')
        # Generate the k-mer profile from the doubled sequence
        self.kmers = self.generate_kmers(self.doubled_sequence, k)

    def generate_kmers(self, sequence, k):
        """
        Generates a k-mer profile from a sequence.

        Args:
            sequence (str): The sequence to process.
            k (int): The size of the k-mers.

        Returns:
            dict: A dictionary mapping each k-mer to its frequency count.
                  Only k-mers with valid DNA characters (A, C, T, G) are included.
        """
        kmers = {}
        valid_chars = {"A", "C", "T", "G"}
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            # Ensure all characters in the k-mer are valid DNA bases
            if all(c in valid_chars for c in kmer):
                kmers[kmer] = kmers.get(kmer, 0) + 1
        return kmers

    def compare_kmer_profiles(self, other):
        """
        Compares the k-mer profile of this sequence with another sequence.

        This method calculates the absolute difference in counts for each k-mer
        between the two profiles. It distinguishes between "perfect" k-mers
        (homopolymers like 'AA', 'CC') and "imperfect" k-mers.

        Args:
            other (Sequence): The other Sequence object to compare against.

        Returns:
            list: A list containing three integers:
                  [total_mismatch_count, perfect_kmer_mismatch, imperfect_kmer_mismatch]
        """
        mismatch_count_all = 0
        mismatch_count_actg = 0

        # Define perfect k-mers (homopolymers of length 2)
        perfect_kmers = {'AA', 'CC', 'TT', 'GG'}

        # Get the union of all k-mers present in either sequence
        keys_union = set(self.kmers) | set(other.kmers)

        for kmer in keys_union:
            # Calculate the absolute difference in frequency for the current k-mer
            diff = abs(self.kmers.get(kmer, 0) - other.kmers.get(kmer, 0))
            mismatch_count_all += diff
            # Check if the k-mer is one of the "perfect" k-mers
            if kmer in perfect_kmers:
                mismatch_count_actg += diff

        # Calculate the mismatch count for "imperfect" k-mers
        mismatch_count_other = mismatch_count_all - mismatch_count_actg

        return [mismatch_count_all, mismatch_count_actg, mismatch_count_other]


def dereplicate_fasta(fasta_file):
    """
    Reads a FASTA file, finds unique sequences, and counts their occurrences.

    Args:
        fasta_file (str): The path to the input FASTA file.

    Returns:
        dict: A dictionary mapping a new unique identifier (e.g., "Uniq0001;size=100;")
              to each unique sequence.
    """
    global total_input_seqs

    seq_counter = Counter()
    seq_chunks = []
    # Read the FASTA file line by line
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # When a new header is found, process the previous sequence
                if seq_chunks:
                    seq = ''.join(seq_chunks)
                    seq_counter[seq] += 1
                    seq_chunks = []
            else:
                # Append sequence lines to the chunk
                seq_chunks.append(line)
        # Process the last sequence in the file
        if seq_chunks:
            seq = ''.join(seq_chunks)
            seq_counter[seq] += 1

    # Sort sequences by abundance in descending order
    sorted_seqs = seq_counter.most_common()
    # Create new headers for unique sequences including their size (count)
    uniq_sequences = {f"Uniq{i + 1:04d};size={count};": seq
                      for i, (seq, count) in enumerate(sorted_seqs)}

    total_input_seqs = sum(seq_counter.values())
    print(f"Total input sequences: {total_input_seqs}")
    print(f"Dereplicated sequences: {len(uniq_sequences)}")

    return uniq_sequences


def sort_variants(seq1, seq2):
    """
    Compares two sequences to determine if one is a noisy variant of the other.

    A sequence (seq2) is classified as noise if it has a similar length,
    is significantly less abundant than the other sequence (seq1), and their
    k-mer profiles are very similar based on predefined thresholds.

    Args:
        seq1 (Sequence): The more abundant sequence (potential parent).
        seq2 (Sequence): The less abundant sequence (potential noisy variant).
    """
    global variants, noise, filtered, sequences, count, fx, ax

    # Condition 1: Length difference must be at most 1 base.
    if abs(seq1.length - seq2.length) <= 1:
        # Condition 2: Abundance ratio must be significant (parent must be >= ax more abundant).
        # Use the global abundance ratio 'ax' (from -a argument)
        if seq1.size / seq2.size >= ax:
            mismatch_count = seq1.compare_kmer_profiles(seq2)
            # Precompute mismatch values for clarity
            p = mismatch_count[1]   # Perfect k-mer mismatch
            im = mismatch_count[2]  # Imperfect k-mer mismatch
            imp = im - p            # Difference between imperfect and perfect mismatches
            len_diff = abs(seq1.length - seq2.length)
            # Condition 3: K-mer profile differences must be within thresholds defined by freedom factor `fx`.
            if p <= 4 * fx and im <= 4 * fx and imp <= 2 * fx:
                # If all conditions are met, classify seq2 as noise
                noise[seq2.title] = seq2.sequence
                filtered[seq2.title] = [seq2.title, seq1.title, "NOISY VARIANT", p, im, len_diff]
                count = count + 1


def write_results(variants, noise, filtered):
    """
    Writes the classification results to output files.

    - asv_tab.txt: A table summarizing the classification of each sequence.
    - variants.fa: A FASTA file of sequences classified as valid variants.
    - noise.fa: A FASTA file of sequences classified as noise.

    Args:
        variants (dict): Dictionary of valid sequence variants.
        noise (dict): Dictionary of noisy sequences.
        filtered (dict): Dictionary containing classification details for the table.
    """
    global weights
    # Write the summary table (ASV table)
    asv_table = 'asv_tab.txt'
    with open(asv_table, 'w') as at:
        at.write('title\tclosest neighbor(s)\tdescription\tperfect k-mer\timperfect k-mer\tlength difference\n')
        for k, v in filtered.items():
            at.write(
                str(v[0]) + '\t' + str(v[1]) + '\t' + str(v[2]) + '\t' + str(v[3]) + '\t' + str(v[4]) + '\t' + str(v[5]) + '\n')

    # Write the FASTA file for valid variants
    out_fasta_variants = 'variants.fa'
    with open(out_fasta_variants, 'w') as of:
        for i in range(len(variants)):
            title = str(list(variants.keys())[i])
            of.write('>' + title + 'neighbors=' + str(weights.get(title, 0)) + ';' + '\n')
            of.write(str(list(variants.values())[i]) + '\n')

    # Write the FASTA file for noisy sequences
    out_fasta_noise = 'noise.fa'
    with open(out_fasta_noise, 'w') as of:
        for i in range(len(noise)):
            of.write('>' + str(list(noise.keys())[i]) + '\n')
            of.write(str(list(noise.values())[i]) + '\n')


def progress_bar(i):
    """
    Displays a simple progress bar in the console.

    Args:
        i (int): The current iteration index.
    """
    progress = (i + 1) / len(sequences)
    bar_length = 40  # Visual length of the progress bar
    filled_length = int(bar_length * progress)
    bar = '#' * filled_length + '-' * (bar_length - filled_length)
    percent = progress * 100
    # Write the progress bar to the same line in the console
    sys.stdout.write(f"\r[{bar}] {percent:.1f}%")
    sys.stdout.flush()


# --- Main Execution ---
if __name__ == "__main__":
    print('\nscript MASV v.1.0.1')
    start_time = time.time()

    # Get arguments from command line
    input_fasta = args.input_fasta
    fx = int(args.freedom)
    # Get the new abundance ratio argument
    ax = args.abundance_ratio
    # Get the new spurious variant saving flag
    save_spurious = args.save_spurious

    # Step 1: Dereplicate the input FASTA file
    uniq_fasta = dereplicate_fasta(input_fasta)

    # Step 2: Create Sequence objects for each unique sequence
    sequences = [Sequence(title, seq) for title, seq in uniq_fasta.items()]

    # Initialize dictionaries to store results
    filtered = {} # Stores detailed classification info for the output table
    variants = {} # Stores sequences identified as true variants
    noise = {}    # Stores sequences identified as noise
    weights = {}  # Stores the count of noisy neighbors for each variant

    print('Identifying sequence variants:')
    # Step 3: Compare all unique sequences against each other
    for i, seq1 in enumerate(sequences):
        count = 0
        # Check if seq1 has not already been classified as noise
        if seq1.title not in filtered.keys():
            # Compare seq1 with all subsequent sequences in the list
            for j, seq2 in enumerate(sequences):
                # Ensure we don't compare a sequence to itself and that seq2 is not already classified
                if i < j and seq2.title not in filtered.keys():
                    sort_variants(seq1, seq2)

            # Store the number of noisy sequences found for seq1
            weights[seq1.title] = count

            # Classify seq1 as a variant if its abundance is 2 or more
            if seq1.size >= 2:
                variants[seq1.title] = seq1.sequence
                filtered[seq1.title] = [seq1.title, '*', "VARIANT", '*', '*', '*']
            # Sequences with size < 2 are not considered variants at this stage
            else:
                pass
        # Update the progress bar after processing each sequence
        progress_bar(i)

    # Step 4: Final classification for any remaining sequences
    # Any sequence not classified as a variant or noise is considered a "spurious variant"
    for seq1 in sequences:
        if seq1.title not in filtered.keys():
            # Use the 'save_spurious' flag (-s argument) to decide where to save the sequence
            if save_spurious:
                # Save to variants output
                variants[seq1.title] = seq1.sequence
            else:
                # Save to noise output (default behavior)
                noise[seq1.title] = seq1.sequence
            # The classification in the table remains the same
            filtered[seq1.title] = [seq1.title, '*', "SPURIOUS VARIANT", '*', '*', '*']

    # Step 5: Write all results to output files
    write_results(variants, noise, filtered)

    # --- Script Completion ---
    end_time = time.time()
    total_time = end_time - start_time
    print(f"\nDetected sequence variants: {len(variants)}")
    print(f"Script execution time: {total_time:.2f} seconds ({total_time / 60:.2f} minutes)\n")

