# MASV
## A high-resolution and transparent Python script for denoising fungal amplicon sequence variants
MASV is a lightweight, transparent Python script for high-resolution denoising of amplicon sequencing data. It is designed to efficiently separate true, low-abundance sequence variants from sequencing noise. 
> It is particularly well-suited for high-resolution analyses, such as deep sequencing of single organisms from cultures or herbarium specimens, where the goal is to resolve the full spectrum of intra-species and allelic variation.
## Installation
No installation is required. MASV is a single standalone script.
Simply download the masv100.py file.
The script is written for Python 3 and only uses standard libraries.
## Usage
Run the script from your terminal. The two required arguments are the input FASTA file (-i) and the degree of freedom (-f).
```bash
python MASV_final.py -i <input.fasta> -f <degree_of_freedom>
```
Example:
```bash
python MASV_final.py -i sample_data.fasta -f 1
```
**Arguments**  
1. -i, --input_fasta (Required): Path to the input FASTA file.  
2. -f, --freedom (Required): An integer representing the "degree of freedom" for noise filtering. This value acts as a multiplier on the core denoising thresholds.  
3. -f 1: (Recommended) The most stringent (default) setting. It uses the empirically-derived base thresholds (4, 4, 2) designed to capture noise originating from a single nucleotide error.  
 -f 2: A more lenient setting. It doubles the thresholds (to 8, 8, 4), allowing variants with approximately two nucleotide differences to be classified as noise.  

**Output Files**  
MASV generates three files in the directory where it is run:  
1. variants.fa  
A FASTA file containing all sequences classified as "VARIANT".
> The header of each sequence includes a neighbors=... tag, which counts how many "NOISY VARIANT" sequences were associated with it.
2. noise.fa  
A FASTA file containing all sequences classified as "NOISY VARIANT" or "SPURIOUS VARIANT".
3. asv_tab.txt  
The "transparency report." This is a tab-delimited log file detailing the classification for every unique sequence.
> It lists the sequence's final description ("VARIANT", "NOISY VARIANT", etc.), its closest neighbor (if any), and the perfect k-mer, imperfect k-mer, and length difference metrics used for the decision.
