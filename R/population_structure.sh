# ======================================= ADMIXTURE ======================================

#	1.	For each of the separate assemblies, we'll remove any individuals that have fewer
#		than 10K SNPs (called "relaxed").

# 		For ATL_MXPL, stringent (n=189)
vcftools --vcf ATL_MXPL.vcf --keep atlmxpl_indskeep.txt --recode --out ATL_MXPL_relaxed
# 		For CENTAM, 18 inds removed because of high MD, leaving relaxed (n=140)
vcftools --vcf new_CENTAM.vcf --keep centam_indskeep.txt --recode --out new_CENTAM_relaxed
# 		For PACMX, 29 inds removed because of high MD, leaving relaxed (n=245)
vcftools --vcf new_PACMX.vcf --keep pacmx_indskeep.txt --recode --out new_PACMX_relaxed
# 		For forreri, no filtering needs to be done (all 104 inds retained)

#	2.	Remove any sites that are more than 75% missing:
# 		kept 707436 out of a possible 718861 Sites
vcftools --vcf ATL_MXPL_relaxed.recode.vcf --max-missing 0.25 --recode --out ATL_MXPL_relaxed_0.25miss
# 		kept 434394 out of a possible 441907 Sites
vcftools --vcf new_CENTAM_relaxed.recode.vcf --max-missing 0.25 --recode --out new_CENTAM_relaxed_0.25miss
# 		kept 767517 out of a possible 775665 Sites
vcftools --vcf new_PACMX_relaxed.recode.vcf --max-missing 0.25 --recode --out new_PACMX_relaxed_0.25miss
# 		kept 546944 out of a possible 564916 Sites
vcftools --vcf forreri.vcf --max-missing 0.25 --recode --out forreri_0.25miss

#	3.	Use Plink to generate a variant-based missing data report (.lmiss):
plink --vcf ATL_MXPL_relaxed_0.25miss.recode.vcf --missing --const-fid --allow-extra-chr --out ATL_MXPL_relaxed_0.25miss
plink --vcf new_CENTAM_relaxed_0.25miss.recode.vcf --missing --const-fid --allow-extra-chr --out new_CENTAM_relaxed_0.25miss
plink --vcf new_PACMX_relaxed_0.25miss.recode.vcf --missing --const-fid --allow-extra-chr --out new_PACMX_relaxed_0.25miss
vcftools --vcf forreri_0.25miss.recode.vcf --missing-site --out forreri_0.25miss

#	4.	Identify which SNPs should be retained after LD-pruning using the `LD-pruning.R` script.

#	5.	Now, remove selected SNPs from the matrix using vcftools.
# 		For ATL_MXPL, LD-pruning kept 50228 out of a possible 707436 Sites
vcftools --vcf ATL_MXPL_relaxed_0.25miss.recode.vcf --positions ATL_MXPL_relaxed_0.25miss_ldpruned.txt --recode --out ATL_MXPL_relaxed_0.25miss_ldp
# 		For CENTAM, LD-pruning kept 36126 out of a possible 434394 Sites
vcftools --vcf new_CENTAM_relaxed_0.25miss.recode.vcf --positions new_CENTAM_relaxed_0.25miss_ldpruned.txt --recode --out new_CENTAM_relaxed_0.25miss_ldp
# 		For PACMX, LD-pruning kept 45130 out of a possible 767517 Sites
vcftools --vcf new_PACMX_relaxed_0.25miss.recode.vcf --positions new_PACMX_relaxed_0.25miss_ldpruned.txt --recode --out new_PACMX_relaxed_0.25miss_ldp
# 		For forreri, kept 65907 out of a possible 546944 Sites
vcftools --vcf forreri_0.25miss.recode.vcf --positions forreri_0.25miss_ldpruned.txt --recode --out forreri_0.25miss_ldp

#	6.	Make a .bed file (and associated bim and fam files) for ADMIXTURE using Plink
plink --vcf ATL_MXPL_relaxed_0.25miss_ldp.recode.vcf --out ATL_MXPL_relaxed_0.25miss_ldp --allow-extra-chr --make-bed --const-fid
plink --vcf new_CENTAM_relaxed_0.25miss_ldp.recode.vcf --out new_CENTAM_relaxed_0.25miss_ldp --allow-extra-chr --make-bed --const-fid
plink --vcf new_PACMX_relaxed_0.25miss_ldp.recode.vcf --out new_PACMX_relaxed_0.25miss_ldp --allow-extra-chr --make-bed --const-fid
plink --vcf forreri_0.25miss_ldp.recode.vcf --out forreri_0.25miss_ldp --allow-extra-chr --make-bed --const-fid

#	7.	Change chromosome names (the first column of the bim file needs to be the chrom name,
#		which has a specified list of names of which none apply to our data). The following just
#		replaces the RAD tag numbers of the first column with 0s.

awk '{$1=0;print $0}' ATL_MXPL_relaxed_0.25miss_ldp.bim > ATL_MXPL_relaxed_0.25miss_ldp.bim.tmp
mv ATL_MXPL_relaxed_0.25miss_ldp.bim.tmp ATL_MXPL_relaxed_0.25miss_ldp.bim

awk '{$1=0;print $0}' new_CENTAM_relaxed_0.25miss_ldp.bim > new_CENTAM_relaxed_0.25miss_ldp.bim.tmp
mv new_CENTAM_relaxed_0.25miss_ldp.bim.tmp new_CENTAM_relaxed_0.25miss_ldp.bim

awk '{$1=0;print $0}' new_PACMX_relaxed_0.25miss_ldp.bim > new_PACMX_relaxed_0.25miss_ldp.bim.tmp
mv new_PACMX_relaxed_0.25miss_ldp.bim.tmp new_PACMX_relaxed_0.25miss_ldp.bim

awk '{$1=0;print $0}' forreri_0.25miss_ldp.bim > forreri_0.25miss_ldp.bim.tmp
mv forreri_0.25miss_ldp.bim.tmp forreri_0.25miss_ldp.bim

#	8.	Run ADMIXTURE for 5 repetitions (per K) with different starting seeds. Run the
#		following within a repetition directory as the code won't append rep number to the
#		P or Q files.

PREFIX="ATL_MXPL_relaxed_0.25miss_ldp"
for K in 3 4 5 6 7 8 9 10 11 12 13 14 15; do admixture --cv ../../$PREFIX.bed $K --seed=$RANDOM | tee $PREFIX.${K}.out; done

PREFIX="new_PACMX_relaxed_0.25miss_ldp"
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; do admixture --cv ../../$PREFIX.bed $K --seed=$RANDOM | tee $PREFIX.${K}.out; done

PREFIX="new_CENTAM_relaxed_0.25miss_ldp"
for K in 3 4 5 6 7 8 9 10 11 12 13 14 15; do admixture --cv ../../$PREFIX.bed $K --seed=$RANDOM | tee $PREFIX.${K}.out; done

PREFIX="forreri_0.25miss_ldp"
for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv ../../$PREFIX.bed $K --seed=$RANDOM | tee $PREFIX.${K}.out; done

#	9.	Rename filenames so they all have replicates appended:

for f in *; do mv -v "$f" "${f%.*}.rep1.${f##*.}"; done