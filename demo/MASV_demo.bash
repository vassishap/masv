#!/bin/bash
##script MASV demo run: denoising of unique amplicon reads (ITS and 16S) and de novo chimera check
## MASV 2.0.0 (github.com/vassishap/masv, commit f0f8673)
## last modification 23 October 2026

{
######### Parameters #####################

ax_default=2				# default 2 (MASV -a)
ax_high_ITS=6				# higher ax used for ITS in the benchmark
ax_high_16S=6.3				# higher ax used for 16S in the benchmark
spurious=False				# default False (MASV -s)
threads=8					# default 2 (MASV -t)

######### Parameters - END ###############

# Build MASV once (needs g++ with C++17)
# git clone https://github.com/vassishap/masv.git
# cd masv && git checkout f0f8673
# g++ -O3 -std=c++17 -pthread masv.cpp -o masv
# cp masv /usr/local/bin/    # or adjust the path below

mkdir -p results

#Introduce the markers to analyse - input files must be named <marker>_uniques.fa
#Headers must look like >Uniq1;size=N;  (e.g. from usearch -fastx_uniques -sizeout -relabel Uniq)
#Only A/C/G/T bases are accepted
for A in ITS 16S
do

printf '\n\n//////////////////////////////// \n   Work on : '$A'\n'

if [ "$A" = "ITS" ]; then ax_high=$ax_high_ITS; else ax_high=$ax_high_16S; fi

for AX in $ax_default $ax_high
do

printf '\n******************************** \n   '$A' with ax = '$AX'\n\n'

# MASV writes asv_tab.txt, variants.fa and noise.fa to the current folder,
# so each run goes into its own folder
mkdir -p results/$A\_ax$AX
cd results/$A\_ax$AX

# Denoise: predict biological sequences (variants) from unique reads
masv -i ../../$A\_uniques.fa -a $AX -s $spurious -t $threads

# Fix headers for usearch (if input headers had no trailing ";", MASV writes ...;size=Nnoise=M;)
sed -E 's/;size=([0-9]+)noise=[0-9]+;?/;size=\1;/' variants.fa > $A\_ax$AX\_variants_sized.fa

# De novo chimera check
usearch -uchime3_denovo $A\_ax$AX\_variants_sized.fa -uchimeout $A\_ax$AX\_uchime.txt -chimeras $A\_ax$AX\_chimeras.fa -nonchimeras $A\_ax$AX\_nonchimeras.fa

# Optional: compare outputs with reference sequences (<marker>_refs.fa)
# usearch -usearch_global $A\_ax$AX\_nonchimeras.fa -db ../../$A\_refs.fa -id 0.97 -strand both -userout $A\_ax$AX\_refhits.txt -userfields query+target+id

cd ../..

done
done
# end of for...

} 2>&1 | tee `echo $0 | cut -d"." --complement  -f2-`_`date +%y%m%d_%H%M`.log
