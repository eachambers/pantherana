# ================== Pooled assembly: RAxML-ng with asc bias correction ===================

# 1.	We're going to run three separate trees using our three concatenated SNP matrices. There are
#		a few steps we need to take before doing so. Import the full (invariant+variant) file into PAUP*

# 2.	Now, ask PAUP* to calculate state frequencies (per sample and mean) using:
> stateFreq
#			This will output a table (which I've saved as "data/pooled_statefreqs.txt").

# 3.	Run `State_freqs.R` code to calculate state frequencies for each of the pooled assemblies.

# 4.	We now need to deal with any sites that may resolve to being invariant (due to ambiguous
#		bases). To do so, run `Remove_invariant_sites.R` on the three SNP datasets.

#			rana_n-1.snps_min10K.phy has 350,292 SNPs
#			Phylip files with potentially invariant sites removed saved as:
#			rana_n-1.snps_min10K_invarrem.phy # 245,016 sites remain

# 4.	Now, we're ready to run RAxML-ng. Let's use the STAM ascertainment bias correction.
#		The numbers specified in the ASC {}s are the number of *invariant* A/C/G/T sites.The 
#		line of code to run will be the following:

./raxml-ng --model GTR+G+ASC_STAM{631110/495387/508960/626585} --all --tree rand --bs-trees 200 --threads 24 --msa /scratch/03123/eac3496/rana_n-1.snps_min10K_invarrem.phy --prefix min10K
