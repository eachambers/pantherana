# ======================================== HHSD ==========================================

# We ran HHSD on three datasets that correspond to three groupings within our 
# phylogenetic tree that have unclear species boundaries. In all cases, individuals were
# assigned populations based on the majority cluster for their ancestry coefficient values.
# Guide trees were based on our phylogeny. The three groups are as follows:
#			(1) The foothills group (magnaocularis, yavapaiensis, omiltemana OA, omiltemana
#				GE, Atenquique short, and Atenquique long) n=79, file naming is "foothills"
#			(2) The northern Atlantic coast and Central Mexican Plateau (berlandieri,
#				neovolcanica, tlaloci, spectabilis, and chichicuahutla) n=116, file naming is "mxpl"
#			(3) Pacific coastal lowlands, i.e., forreri complex (forreri+cora+adleri, floresi, 
#				hillisi, miadis, Arcelia) n=104, file naming is "forreri"

#	1.	Create a mamba env to run HHSD in
mamba create -c conda-forge -c bioconda -n hhsd3.9
mamba activate hhsd3.9
mamba install python=3.9.19

# needs python env >3.9 and <3.13
mamba create -c conda-forge -c bioconda -n hhsd3.12
mamba activate hhsd3.12
mamba install python=3.12.9

#	2.	Install HHSD v0.9.9
git clone https://github.com/abacus-gene/hhsd
cd hhsd
pip install .

#	3.	We need to reduce the sequence data for computational time's sake, and we need the
#		full locus (rather than just the SNPs) for each SNP that we want to include in 
#		this analysis. We also need to reduce the sequence data for computational time's 
#		sake -- but HHSD also performs better with more data so we don't want to go nuts
#		here -- so we'll remove any sites with missing data above 40%:
vcftools --vcf forreri_0.25miss_ldp.recode.vcf --max-missing 0.6 --recode --out forreri_0.6miss_ldp # kept 12453 out of a possible 65000 sites
vcftools --vcf ATL_MXPL_relaxed_0.25miss_ldp.recode.vcf --max-missing 0.6 --recode --out ATL_MXPL_0.6miss_ldp # kept 4684 out of a possible 50228 sites; 0.5 keeps 14332 / 50228 sites
vcftools --vcf new_PACMX_relaxed_0.25miss_ldp.recode.vcf --max-missing 0.6 --recode --out PACMX_0.6miss_ldp # kept 8226 out of a possible 45130 sites; 0.5 keeps 17616 / 45130 sites

#	4.	For foothills and mxpl, we need to generate new vcfs with subsets of the assemblies
#		containing relevant individuals. To do so:
vcftools --vcf ATL_MXPL_0.6miss_ldp.recode.vcf --keep mxpl_indv.txt --recode --out mxpl_06 # kept 116 inds
vcftools --vcf PACMX_0.6miss_ldp.recode.vcf --keep foothills_indv.txt --recode --out foothills_06 # kept 79 inds

#	5.	We need our phylip input file to have the full locus (rather than just the SNPs).
#		So, we'll get the RAD tag numbers and positions so we can retrieve the full loci
#		from the .loci file.
#			(a) Get a missing data report on the recoded vcf which we'll use to retrieve
#				locus and SNP names:
plink --vcf forreri_0.6miss_ldp.recode.vcf --missing --const-fid --allow-extra-chr --out forreri_0.6miss
plink --vcf foothills_06.recode.vcf --missing --const-fid --allow-extra-chr --out foothills_06
plink --vcf mxpl_06.recode.vcf --missing --const-fid --allow-extra-chr --out mxpl_06

#			(b) Using the lmiss file generated from Plink, run the `hhsd.R` script which 
#				will generate files named:
#					forreri_0.6miss_ldp_locusnames.txt
#					foothills_locusnames.txt & foothills_remove.txt & foothills_indv.txt
#					mxpl_locusnames.txt & mxpl_remove.txt & mxpl_indv.txt

#				It will also make an Imap file which specifies individual assignment into 
#				populations/lineages:
#			forreri-Imap.txt is for the ADMIXTURE results from K=5. Populations specified are:
#								forr=typotypic forreri, from N Mexico; flor=floresi+cora;
#								hilli=hillisi; miad=samples southward of hillisi, including LCI;
#								arce=Arcelia form

#		Pops from foothills: 	omig (omiltemana guerrero), omio (omiltemana oaxaca), 
#								magn (magnaocularis), yava (yavapaiensis), 
#								ates (Atenquique short form), atel (Atenquique long form)

#				Args for this script are `lmiss_prefix` and `dataset_name`
Rscript hhsd.R "foothills_06" "foothills"

#			(c) If you're using a subset of individuals (i.e., for the foothills or mxpl 
#				datasets), you need to select desired samples from the iPyrad *.loci file. 
#				To do so:
grep -vf foothills_remove.txt new_PACMX.loci > foothills.loci
grep -vf mxpl_remove.txt ATL_MXPL.loci > mxpl.loci

#			(c) Then, run the `process_loci_file.py` script which will generate separate txt 
#				files for each relevant locus from the iPyrad .loci file in the iPyrad 
#				output files directory:
mkdir loci_foothills
mkdir loci_mxpl
mkdir loci_forreri0.6
python process_loci_file.py loci_foothills/foothills_locus#.txt foothills_06_locusnames.txt foothills.loci
python process_loci_file.py loci_mxpl/mxpl_locus#.txt mxpl_06_locusnames.txt mxpl.loci
python process_loci_file.py loci_forreri0.6/forreri_0.6miss_ldp_locus#.txt forreri_0.6miss_ldp_locusnames.txt forreri.loci

#			(d)	Finally, run the following on all the files produced from `process_loci_file.py`
#				which will remove the unnecessary locus labeling lines in each of the txt 
#				files:
# 					If on linux, within each loci_* directory:
for file in *.txt; do sed -i '/|/d' $file; done

#			(e) Get the number of individuals (i.e., line numbers) and the number of sites
#				in each locus for the Phylip file.

#				To do so on a linux machine, run:
for file in loci_foothills/*.txt; do sed -i 's/^/^/' $file && chars=`head -n 1 $file | wc -m` && inds=`cat $file | wc -l` && tmp=`echo "$inds $(($chars-24))"` && sed -i "1s/^/$tmp\n/" $file && sed -i 's/^[ \t]*//' $file; done
for file in loci_mxpl/*.txt; do sed -i 's/^/^/' $file && chars=`head -n 1 $file | wc -m` && inds=`cat $file | wc -l` && tmp=`echo "$inds $(($chars-26))"` && sed -i "1s/^/$tmp\n/" $file && sed -i 's/^[ \t]*//' $file; done

#			(f) We're finally ready to make our Phylip file! Concatenate all the files 
#				together:
cat loci_forreri0.6/*.txt > forreri_0.6miss_ldp.phy
cat loci_foothills/foothills*.txt > foothills_06.phy
cat loci_mxpl/mxpl*.txt > mxpl_06.phy

#				At this stage, it's good to manually go in and verify that your phylip file
#				looks good.

#	6.	Make control files for all three HHSD runs. For each taxonomy (i.e., each of three),
#		we ran both split and merge algorithms and ran HHSD with three migration rate priors
#		to test the impact of priors on downstream results. Files are named according
#		to parameters run, e.g.:
#			"cf_mxpl_merge_0110.txt": control file for mxpl taxonomy, merge algorithm, mig rate prior (0.1, 10)
#			"cf_forreri_merge_220.txt": control file for mxpl taxonomy, merge algorithm, mig rate prior (2, 20)
#			"cf_foothills_split_0120.txt": control file for mxpl taxonomy, split algorithm, mig rate prior (0.1, 20)

#	7.	Run HHSD. Activate the mamba env and then run the control file, for example:
mamba activate hhsd3.9
hhsd --cfile cf_mxpl_merge_0110.txt

# ============= MIGRATION RATE PRIOR TESTING

# Even stricter missing data filters for migration rate prior testing:

vcftools --vcf new_PACMX_relaxed_0.25miss_ldp.recode.vcf --max-missing 0.7 --recode --out PACMX_0.7miss_ldp # kept 1835 out of a possible 45130 sites
vcftools --vcf PACMX_0.7miss_ldp.recode.vcf --keep foothills_indv.txt --recode --out foothills_07 # kept 79 inds
plink --vcf foothills_07.recode.vcf --missing --const-fid --allow-extra-chr --out foothills_07
Rscript hhsd.R "foothills_07" "foothills"
mkdir loci_foothills_07
python process_loci_file.py loci_foothills_07/foothills_locus#.txt foothills_07_locusnames.txt foothills.loci
for file in *.txt; do sed -i '/|/d' $file; done
for file in loci_foothills_07/*.txt; do sed -i 's/^/^/' $file && chars=`head -n 1 $file | wc -m` && inds=`cat $file | wc -l` && tmp=`echo "$inds $(($chars-24))"` && sed -i "1s/^/$tmp\n/" $file && sed -i 's/^[ \t]*//' $file; done
cat loci_foothills_07/foothills*.txt > foothills_07.phy

# Randomly sample same number of loci
vcftools --vcf new_PACMX_relaxed_0.25miss_ldp.recode.vcf --keep foothills_indv.txt --max-missing 0.25 --recode --out foothills_inds # kept 79 inds
plink --vcf foothills_inds.recode.vcf --const-fid --allow-extra-chr --make-bed --recode --out foothills_inds
sort -R foothills_inds.bim | head -1835 | awk '{print $1 "\t" $2}' > foot_rand.snps
plink --vcf foothills_inds.recode.vcf --extract foot_rand.snps --recode vcf --out foothills_rand --const-fid --allow-extra-chr
bcftools query -f '%POS\n' foothills_rand.vcf | wc -l # 1835
plink --vcf foothills_rand.vcf --missing --const-fid --allow-extra-chr --out foothills_rand
Rscript hhsd.R "foothills_rand" "foothills"
mkdir loci_foothills_rand
python process_loci_file.py loci_foothills_rand/foothills_locus#.txt foothills_rand_locusnames.txt foothills.loci
for file in *.txt; do sed -i '/|/d' $file; done
for file in loci_foothills_rand/*.txt; do sed -i 's/^/^/' $file && chars=`head -n 1 $file | wc -m` && inds=`cat $file | wc -l` && tmp=`echo "$inds $(($chars-24))"` && sed -i "1s/^/$tmp\n/" $file && sed -i 's/^[ \t]*//' $file; done
cat loci_foothills_rand/foothills*.txt > foothills_rand.phy


#vcftools --vcf forreri_0.25miss_ldp.recode.vcf --max-missing 0.7 --recode --out forreri_0.7miss_ldp # kept 2940 out of a possible 65000 sites
vcftools --vcf forreri_0.25miss_ldp.recode.vcf --max-missing 0.75 --recode --out forreri_0.75miss_ldp # kept 848 out of a possible 65000 sites
plink --vcf forreri_0.75miss_ldp.recode.vcf --missing --const-fid --allow-extra-chr --out forreri_075
Rscript hhsd.R "forreri_075" "forreri"
mkdir loci_forreri_075
python process_loci_file.py loci_forreri_075/forreri_locus#.txt forreri_075_locusnames.txt forreri.loci
for file in *.txt; do sed -i '/|/d' $file; done
for file in loci_forreri_075/*.txt; do sed -i 's/^/^/' $file && chars=`head -n 1 $file | wc -m` && inds=`cat $file | wc -l` && tmp=`echo "$inds $(($chars-24))"` && sed -i "1s/^/$tmp\n/" $file && sed -i 's/^[ \t]*//' $file; done
cat loci_forreri_075/forreri*.txt > forreri_075.phy


vcftools --vcf ATL_MXPL_relaxed_0.25miss_ldp.recode.vcf --max-missing 0.675 --recode --out ATL_MXPL_0.675miss_ldp # 1202 loci remain
vcftools --vcf ATL_MXPL_0.675miss_ldp.recode.vcf --keep mxpl_indv.txt --recode --out mxpl_0675 # kept 116 inds
plink --vcf mxpl_0675.recode.vcf --missing --const-fid --allow-extra-chr --out mxpl_0675
Rscript hhsd.R "mxpl_0675" "mxpl"
mkdir loci_mxpl_0675
python process_loci_file.py loci_mxpl_0675/mxpl_locus#.txt mxpl_0675_locusnames.txt mxpl.loci
# Move to loci_ dir
for file in *.txt; do sed -i '/|/d' $file; done
# Move back out of loci dir
for file in loci_mxpl_0675/*.txt; do sed -i 's/^/^/' $file && chars=`head -n 1 $file | wc -m` && inds=`cat $file | wc -l` && tmp=`echo "$inds $(($chars-26))"` && sed -i "1s/^/$tmp\n/" $file && sed -i 's/^[ \t]*//' $file; done
cat loci_mxpl_0675/mxpl*.txt > mxpl_0675.phy



# ============= TROUBLESHOOTING

#### Run within HHSD directory
hhsd/build/lib/hhsd/bpp/linux

bpp v4.7.0_linux_x86_64, 126GB RAM, 32 cores
https://github.com/bpp/bpp

Detected CPU features: mmx sse sse2 sse3 ssse3 sse4.1 sse4.2 popcnt avx avx2
Auto-selected SIMD ISA: AVX2


Starting timer..
Seed: 21676 (fixed by user)
Parsing species tree... Done
Unable to open file (/media/WangLab/WangLab/mxrana/hhsd/mxpl_06.phy)

../../hhsd/build/lib/hhsd/bpp/linux/bpp --cfile proposed_ctl.ctl

# ========================================================================================

wget https://github.com/bpp/bpp/releases/download/v4.8.4/bpp-4.8.4-linux-x86_64.tar.gz
tar zxvf bpp-4.8.4-linux-x86_64.tar.gz
rm bpp-4.8.4-linux-x86_64.tar.gz

# Check bpp control file generated by HHSD by running (from within bpp/bin dir):
./bpp --cfile ../../

bpp v4.8.4_linux_x86_64, 126GB RAM, 32 cores
https://github.com/bpp/bpp

*TLDR*

The 'mcmcfile' and 'outfile' options are now obsolete. Only the new 'jobname'
option needs to be specified, and all output files, including the MCMC sample
file, will use this prefix.

Details:
-------

Starting with BPP v4.8.0, the 'mcmcfile' and 'outfile' options have been removed.
Now, only the new 'jobname' option needs to be specified, serving as the prefix
for filenames, including the main output file and the MCMC sample file.

Example:
--------
  jobname = out

produces at minimum the following three files:

out.txt               : output file
out.mcmc.txt          : MCMC file
out.SeedUsed          : used seed

################

proposed_ctl.ctl:17: Error when processing option 'finetune'.

The syntax for the 'finetune' tag has changed since BPP v4.8.1, adopting a
dictionary-like key:val format where the key specifies the step length and
val the corresponding starting value.

Examples:

             +-> step lenghts are fixed to the specified values and are not auto adjusted during burnin
             |
  finetune = 0 Gage:5 Gspr:0.001 mix:0.3
  finetune = 1 Gage:5 Gspr:0.001 mix:0.3
             |
             +-> step lengths are set to the specified initial values and are adjusted during burnin

  finetune = 0    # (no fine tuning, use default values, *NOT RECOMMENDED*)
  finetune = 1    # (fine tuning, use default inital values, default setting)