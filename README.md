# MASV (current version 1.0.1)
## A high-resolution and transparent Python script for denoising amplicon sequence variants
MASV is a lightweight, transparent Python script for high-resolution denoising of amplicon sequencing data. It is designed to efficiently separate true, low-abundance sequence variants from sequencing noise. 
>  This script was originally designed for analsys of the data coming from Illumina based sequencing of fungal ITS region. However, the core algorithm is not marker-specific and can be used for analyzing amplicon-based HTS data from any molecular marker or organism.

It is particularly well-suited for high-resolution analyses, such as deep sequencing of single organisms from cultures or herbarium specimens, where the goal is to resolve the full spectrum of intra-species and allelic variation.

## Installation
The script is written for Python 3 and only uses standard libraries.  
No installation is required. MASV is a single standalone script.  
Simply download the masv101.py file.  
## General workflow  
![My Project Logo](./wf_masv.png)  
> Full information on MASV's algorithm workflow can be found in this publication: link
## Usage
Run the script from your terminal. The only required arguments is the input FASTA file (-i).
```bash
python masv101.py -i <input.fasta>
```
**Arguments**  
-i, --input_fasta (Required): Path to the input FASTA file.

-f, --freedom (Optional): An integer representing the "degree of freedom" for noise filtering. This value acts as a multiplier on the core denoising thresholds.

>-f 1: (Default) The most stringent setting. It uses the empirically-derived base thresholds (see note below) designed to capture noise originating from a single nucleotide error.
>
>-f 2: A more lenient setting. It doubles the thresholds, allowing variants with approximately two nucleotide differences to be classified as noise.

-a, --abundance_ratio (Optional): A float representing the minimum abundance ratio (Parent Size / Noise Size) required to classify a sequence as noise. The default is 1.45.

-s, --save_spurious (Optional): A boolean flag (True/False) that controls the output of sequences classified as "SPURIOUS VARIANT" (singletons that did not match a parent). The default is False.

Example of the default settings:
```bash
python masv101.py -i sample_data.fasta -f 1 -a 1.45 -s False
```

**Singletons Preservation (-s flag)**  
The -s True option can be used to preserve unclassified singleton sequences ("SPURIOUS VARIANT") by writing them to the variants.fa output file.

Example to preserve singletons:
```bash
python masv101.py -i sample_data.fasta -f 1 -a 1.45 -s True
```

**Output Files**  
MASV generates three files in the directory where it is run:  
1. variants.fa  
A FASTA file containing all sequences classified as "VARIANT".
> The header of each sequence includes a neighbors=... tag, which counts how many "NOISY VARIANT" sequences were associated with it. Sequences classified as "SPURIOUS VARIANT" are also added here if the -s True flag is used.
2. noise.fa  
A FASTA file containing all sequences classified as "NOISY VARIANT".
> Sequences classified as "SPURIOUS VARIANT" are added here if the -s False flag is used (default behavior).
3. asv_tab.txt  
The "transparency report." This is a tab-delimited log file detailing the classification for every unique sequence.
> It lists the sequence's final description ("VARIANT", "NOISY VARIANT", "SPURIOUS VARIANT"), its closest neighbor (if any), and the perfect k-mer, imperfect k-mer, and length difference metrics used for the decision.  

## A note on thresholds (-f and -a parameters)

The core thresholds (4, 4, 2) used by MASV are not arbitrary. They were empirically derived from an *in silico* simulation (see details in the publication).  
**Method**: 1,500 "child" sequences were generated from a "parent" *Russula* sp. ITS2 sequence.  Each child contained exactly one random SNP (n=750) or InDel (n=750).  
**Analysis**: Each of the 1,500 artificial "child" sequences was compared to the original, unaltered parent sequence.  
During each comparison, we calculated the three core metrics used by MASV to quantify sequence dissimilarity.  
*	Perfect k-mer difference (d[p]): The absolute difference in counts of 'perfect' k-mers (i.e., 'AA', 'CC', 'GG', 'TT').
*	Imperfect k-mer difference (d[im]): The absolute difference in counts of all other 'imperfect' k-mers (e.g., 'AC', 'AT', 'CA', etc.).
*	Net difference (d[net]): The value of the imperfect k-mer difference minus the perfect k-mer difference (d[im]−d[p]).   

**Result**: The maximum observed distortion from a single nucleotide error defined the base thresholds.
* Max d[p]: 4 
* Max d[im]: 4 
* Max d[net]: 2

Therefore, running MASV with -f 1 sets the filter to its most stringent, calibrated level, designed specifically to remove noise that is consistent with a single sequencing error.  

**IMPORTANT**

Long-read sequencing technologies (e.g., PacBio, Oxford Nanopore) possess substantially different error profiles and rates. Direct application of MASV to data from these or other non-Illumina platforms will likely result in suboptimal performance.

Recommendation: To maintain optimal resolution and accuracy when applying MASV to non-Illumina data, the platform-specific parameters (e.g., abundance-ratio, k-mer thresholds) should be re-evaluated and adjusted.

## Citation  
If you use MASV in your research, please cite:  
> MASV: A high-resolution and transparent Python script for denoising fungal amplicon sequence variants Vasilii Shapkin, Miroslav Kolařík, Petr Kohout, Tomáš Větrovský [Journal, Year]

