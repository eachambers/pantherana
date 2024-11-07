# The following is a thorough walkthrough of all steps involved in the processing
# of Rana data

# ========================================================================================
# ================ running iPyrad to obtain pooled assembly (all samples) ================

# iPyrad v.0.9.85

# Ran each pool within each project separately for iPyrad steps 1-3:
#	JA19241 (n=130); pools 1-6
#	JA19242 (n=311); pools 1-13
#	JA20247 (n=20); pool 1
#	JA20248 (n=163); pools 1-7
#	Copied fastqs from Chambers et al. 2023 2bRAD project and re-ran steps 2-3 (n=9) ("epirana")

# Merged everything together with the following line (n=XX 633):
ipyrad -m all_rana params-JA19241_pool1_R1.txt params-JA19242_pool11_R1.txt params-JA19242_pool5_R1.txt  params-JA20248_pool2.txt params-JA19241_pool2_R1.txt params-JA19242_pool12_R1.txt params-JA19242_pool6_R1.txt params-JA20248_pool3.txt params-JA19241_pool3_R1.txt params-JA19242_pool13_R1.txt params-JA19242_pool7_R1.txt params-JA20248_pool4.txt params-JA19241_pool4_R1.txt params-JA19242_pool1_R1.txt params-JA19242_pool8_R1.txt params-JA20248_pool5.txt params-JA19241_pool5_R1.txt params-JA19242_pool2_R1.txt params-JA19242_pool9_R1.txt params-JA20248_pool6.txt params-JA19241_pool6_R1.txt params-JA19242_pool3_R1.txt params-JA20247_pool1_R1.txt params-JA20248_pool7.txt params-JA19242_pool10_R1.txt params-JA19242_pool4_R1.txt params-JA20248_pool1.txt params-epirana.txt

# Had to branch off and remove TF8608_Rber_CMX because there were no clusters formed; 
#	n now = 633 (called "rana_n-1")

# ========================================================================================
# =============== Getting basic data characteristics from pooled assembly ================

# 1.	Convert full phylip file (rana_n-1.phy) to Nexus for input into PAUP* (I did so manually).
#		Do the same for the .snps file for running pooled assembly in RAxML-ng (rana_n-1.snps.nexus).

# 2.	Import file into PAUP*: (this may take a while as the matrix size is large)
> execute ~/Desktop/rana_n-1.nexus

# Processing of file "~/Desktop/rana_n-1.nexus" begins...
# Data read in DNA format

# Data matrix has 632 taxa, 2262042 characters
# Valid character-state symbols: ACGT
# Missing data identified by 'N'
# Gaps identified by '-'
# "Equate" macros in effect:
#    R,r ==> {AG}
#    Y,y ==> {CT}
#    M,m ==> {AC}
#    K,k ==> {GT}
#    S,s ==> {CG}
#    W,w ==> {AT}
#    H,h ==> {ACT}
#    B,b ==> {CGT}
#    V,v ==> {ACG}
#    D,d ==> {AGT}

# Processing of input file "rana_n-1.nexus" completed.
	
# 3.	We can take a look at the number of invariant sites by doing:
> cstatus

# Character-status summary:
#   Current optimality criterion = parsimony
#   No characters are excluded
#   Of 2262042 total characters:
#     All characters are of type 'unord'
#     All characters have equal weight
#     2015255 characters are constant (proportion = 0.890901)
#     69627 variable characters are parsimony-uninformative
#     Number of parsimony-informative characters = 177160

# 4.	We also want to take a look at missing data proportions. To do so, we'll import the SNPs
# 		Nexus file into PAUP*, and do
> missdata
	
# Data matrix has 632 taxa, 350292 characters
	
# ========================================================================================
# ================== missing data calculations and removal of samples ====================

# **** BELOW ACTUALLY WASN'T DONE FOR THE POOLED ASSEMBLY BECAUSE RAXML WAS ONLY RUN ON THE SNP DATASET ****

# 1.	Remove individuals based on per-individual missing data proportions using
#		vcftools v0.1.13.

#	First remove all samples other than 100% missing data (n=624).

#	Generate a filtered dataset that retains individuals as long as they have a minimum of 
#	10K SNPs (n=555 / XXX)

#	Below is for 80p dataset (n=414); remove from already recoded min10K dataset ^
./vcftools --vcf /scratch/03123/eac3496/rana_n-1_min10K.recode.vcf --remove-indv TF7108_Rsie_RJO --remove-indv TF8616_Rjoh_RJO --remove-indv JHT4103_LCA --remove-indv V2752_PAC --remove-indv V2750_PAC --remove-indv JHT4106_LCA --remove-indv TF7112_Rzwe_RJO --remove-indv TF7113_Rzwe_RJO --remove-indv TF8605_Rjoh_RJO --remove-indv V2727_PAC --remove-indv TF7131_Rsie_RJO --remove-indv TF7130_Rsie_RJO --remove-indv JHT4110_LCA --remove-indv V2832_PAC --remove-indv TF8614_Rjoh_RJO --remove-indv JHT4095_LCA --remove-indv TF7109_Rsie_RJO --remove-indv V3076_PAC --remove-indv V2787_PAC --remove-indv TF8615_Rjoh_RJO --remove-indv TF7110_Rsie_RJO --remove-indv V2785_PAC --remove-indv JHT4096_LCA --remove-indv T3809_PAC --remove-indv T24985_Rmon_RJO --remove-indv V2784_PAC --remove-indv V2781_PAC --remove-indv N918_LCA --remove-indv T3591_LCA --remove-indv V2783_PAC --remove-indv T3042_LCA --remove-indv V2789_PAC --remove-indv V2788_PAC --remove-indv T3568_LCA --remove-indv MVZ264277_LCA --remove-indv T3698_PAC --remove-indv T1969_Rpus_RJO --remove-indv TF8603_Rjoh_RJO --remove-indv JHT4114_LCA --remove-indv T3555_LCA --remove-indv T3762_LCA --remove-indv T3782_PAC --remove-indv T14192_papa_PAC --remove-indv JHT2309_LCA --remove-indv T3688_PAC --remove-indv T3586_LCA --remove-indv T3678_CMX --remove-indv Rber_T1114_JA17105 --remove-indv V3082_PAC --remove-indv KEN177_LCA --remove-indv Rchi_T2034a_JA17105 --remove-indv CLP2570_LCA --remove-indv T2025_PAC --remove-indv TF8534_Rmon_RJO --remove-indv T3658_LCA --remove-indv T3600_LCA --remove-indv T11_CMX --remove-indv T2103_CMX --remove-indv JHT2312_LCA --remove-indv JHT4122_LCA --remove-indv T3616_PAC --remove-indv V1_PAC --remove-indv T134_Rspe_CMX --remove-indv T14316_Rpip_OUT --remove-indv JHT2760_LCA --remove-indv Rbla_D2865_JA17105 --remove-indv TF8540_Rber_CMX --remove-indv V2761_PAC --remove-indv TF8541_Rber_CMX --remove-indv TF8532_Rmon_RJO --remove-indv T14198_PAC --remove-indv TF7140_Rspe_CMX --remove-indv V5_PAC --remove-indv MVZ269990_LCA --remove-indv T2022_Aten_PAC --remove-indv JHT4091_LCA --remove-indv T3649_LCA --remove-indv T3651_LCA --remove-indv JA24801_PAC --remove-indv LK2_LCA --remove-indv T149_Rmon_RJO --remove-indv T141_Rspe_CMX --remove-indv T3554_LCA --remove-indv JHT2327_LCA --remove-indv TF8535_Rmon_RJO --remove-indv JHT4124_LCA --remove-indv V12_PAC --remove-indv Rbla_D2864_JA17105 --remove-indv T14197_Rspe_CMX --remove-indv T3571_LCA --remove-indv AH198_LCA --remove-indv T14236_Romi_CMX --remove-indv V2659_pap_PAC --remove-indv MVZ269989_LCA --remove-indv KEN592_LCA --remove-indv V2646_pap_PAC --remove-indv JHT2956_LCA --remove-indv JHT4100_LCA --remove-indv T3561_LCA --remove-indv T14318_Rare_OUT --remove-indv T210_Atoy_CMX --remove-indv MVZ264278_LCA --remove-indv V3077_PAC --remove-indv KEN590_LCA --remove-indv JHT3828_LCA --remove-indv V15_PAC --remove-indv V2695_Rfor_PAC --remove-indv UTA56347_Rspe_CMX --remove-indv MVZ263760_LCA --remove-indv MVZ145471_LCA --remove-indv JHT2314_LCA --remove-indv KEN591_LCA --remove-indv V2906_PAC --remove-indv UTA60803_Rfor_PAC --remove-indv T3653_LCA --remove-indv MVZ143381_LCA --remove-indv T14298_Ryav_OUT --remove-indv Rsph_T25870_JA17105 --remove-indv JHT2328_LCA --remove-indv JHT2308_LCA --remove-indv MVZ150254_LCA --remove-indv T208_CMX --remove-indv JHT2560_LCA --remove-indv V1311_PAC --remove-indv JHT4113_LCA --remove-indv V2764_PAC --remove-indv JHT2329_LCA --remove-indv Rsph_T26064_JA17105 --remove-indv JHT4107_LCA --remove-indv T3572_LCA --remove-indv JHT3762_LCA --remove-indv T3654_LCA --remove-indv TF8592_Rber_CMX --remove-indv JA22541_Rspe_CMX --remove-indv T147_Rmon_RJO --remove-indv T2057_Rchi_OUT --remove-indv JHT3830_LCA --remove-indv TF8533_Rmon_RJO --remove-indv JHT3829_LCA --remove-indv T136_Rspe_CMX --remove-indv JHT2141_LCA --recode --out /scratch/03123/eac3496/rana_n-1_max80p

# 3.	Convert from vcf to phylip (and nexus) to run RAxML using the vcf2phylip.py script 
#		(https://github.com/edgardomortiz/vcf2phylip)

#	Can also include the --nexus flag for it to output a nexus file too; can specify outgroup using -o flag so that it
#	places that taxon at the top of the file (RAxML typically roots using the first sequence)
#	Default is that min_inds_locus is 4 so just leave as is

python vcf2phylip.py -i input.vcf --nexus
	
# ========================================================================================
# ================== Pooled assembly: RAxML-ng with asc bias correction ===================

# 1.	We're going to run three separate trees using our three concatenated SNP matrices. There are
#		a few steps we need to take before doing so. Import the full (invariant+variant) file into PAUP*

# 2.	Now, ask PAUP* to calculate state frequencies (per sample and mean) using:
> stateFreq
	
#	It will output a table (which I've saved as pooled_statefreqs.txt).
#	Run State_freqs.R code to calculate state frequencies for each of the pooled assemblies.

# Total dataset size:		2,262,042 sites

### min_500
#Mean (A/C/G/T)			0.28  		0.22  		0.23  		0.27
#Number of A/C/G/Ts:		633371.8	497649.2	520269.7	610751.3
#Rounded					633372		497649		520270		610751

### min_10K
???

### max_80p
#Mean					0.280  		0.22  		0.22  		0.28
#Number of A/C/G/Ts:		633371.8	497649.2	497649.2	633371.8
#Rounded					633372		497649		497649		633372

# 3.	We now need to deal with any sites that may resolve to being invariant (due to ambiguous
#		bases). To do so, run Remove_invariant_sites.R on the three SNP datasets.
	
#	rana_n-1.snps_min10K.phy and rana_n-1.snps_max80p.phy both have 350,292 SNPs

#	Phylip files with potentially invariant sites removed saved as:
#	rana_n-1.snps_min10K_invarrem.phy # Removed XXX??? = 245,016 sites remain
#	rana_n-1.snps_max80p_invarrem.phy # Removed 139,073 of 350,267 sites = 211,194 sites remain <- didn't do this tree though
	
# 4.	Now, we're ready to run RAxML-ng. Let's use the STAM ascertainment bias correction.
#		The numbers specified in the ASC {}s are the number of *invariant* A/C/G/T sites.The 
#		line of code to run will be the following:

./raxml-ng --model GTR+G+ASC_STAM{631110/495387/508960/626585} --all --tree rand --bs-trees 200 --threads 24 --msa /scratch/03123/eac3496/rana_n-1.snps_min10K_invarrem.phy --prefix min10K

============================= Simple trees NJ trees in PAUP ============================

# Recall that you can always do >help commands in PAUP for it to spit out all available commands, and then
# you can do XXX to see further information/details on a specific command.

# 1.	Convert phy to nexus file (using Aliview, vcf2phylip above, or manually).

# 2.	Remove individuals from nexus file using Aliview (if using SNPs datasets; you'll have to
#		use vcf2phylip above if you're using allsites because the files are too large).

# 3.	Import desired SNPs nexus file into PAUP using execute command.

# 4.	Run >nj to calculate a NJ tree on your dataset.

# 5.	Once tree is done, can do the >saveTrees command and it'll save the tree in same wd as
#		the dataset (in my case, this is usually the desktop).

# ========================================================================================
# ============== Running iPyrad on separate assemblies (subsets of samples) ==============

# 1.	We now want to split up samples based on the phylogenetic tree, building assemblies with
#		subsets of closely related samples (to maximize phylo info obtained). To do so, we moved 
#		individual sorted/demultiplexed fastqs over to a separate directory and made a new iPyrad 
#		project from there, specifying the sorted fastq path (see params files in supp materials).
	
# example: ipyrad -p params-ATL_MXPL -s 123456

#	Please refer to the supplementary materials for which samples were assigned into each 
#	separate assembly. The separate assemblies are as follows:
	
#	PACMX		Pacific coast (lowlands and foothills) of Mexico. Includes described species:
#				R. forreri (MX), R. macroglossa (MX), R. spectabilis (Oaxaca only),
#				R. omiltemana, R. magnaocularis, R. yavapaiensis
#				Also includes unnamed species (Hillis & Wilcox 2005): sp. 7 (Jalisco) and 
#				sp. 8 (Puebla), as well as the "papagayo form" described by Hillis et al. 
#				(1983) and Arcelia and Colima forms.
				
#				Total of XXX samples
	
#	CENTAM		Samples occurring in Central America. Includes described species:
#				R. lenca, R. miadis, R. macroglossa (CENTAM), R. forreri (CENTAM), R. taylori
#				Also includes unnamed species (Hillis & Wilcox 2005): sp. 4 (PA), 
#				sp. 5 (CR), sp. 6 (CR)
				
#				Total of XXX samples

#	ATL_MXPL	Atlantic coast and the Mexican Plateau. Includes described species:
#				R. taylori, R. macroglossa, R. brownorum [not explicitly IDed], 
#				R. berlandieri, R. spectabilis, R. tlaloci, R. neovolcanica, 
#				R. chichicuahutla [not explicitly IDed].
#				Also includes unnamed species (Hillis & Wilcox 2005): sp. 3 (which is 
#				R. spectabilis).
	
#				Total of XXX samples

#	forreri		All members from the Pacific coastal lowlands. Includes described species:
#				R. miadis, R. forreri, R. cora, R. adleri, R. hillisi, R. floresi
				
#				Total of XXX samples

# 2.	For each of these assemblies, we generated the state frequency file, character status, and
#		missing data proportions in PAUP* using from a Nexus file for the SNP matrix. We also 
#		processed missing data and calculated average read depth using the Basic_stats.R script.

# ========================================================================================
# ======================================= ADMIXTURE ======================================

#	1.	For each of the separate assemblies, generate two subsets of data, one that removes 
#	inds with 
#		-	<=10K SNPs (relaxed)
#		-	<=100K SNPs (stringent)

# TODO can just use --keep or --remove with a list of inds; easier than doing below

# For ATL_MXPL, stringent (n=189)
vcftools --vcf ATL_MXPL.vcf --remove-indv T14179_Rber_CMX --remove-indv T3935_Rspe_CMX --remove-indv N919_LCA --remove-indv TF8679_Rspe_CMX --remove-indv T416_RTex_CMX --remove-indv TF8677_Rspe_CMX --remove-indv N997_LCA --remove-indv MVZ3720_LCA --remove-indv T3691_PAC --remove-indv N1020_LCA --remove-indv T130_Rspe_CMX --remove-indv MG66_LCA --remove-indv T2092_xochi_CMX --recode --out ATL_MXPL_relaxed

# For CENTAM, 18 inds removed because of high MD, leaving relaxed (n=140)
vcftools --vcf new_CENTAM.vcf --remove-indv T2757_LCA --remove-indv T3548_LCA --remove-indv T3594_LCA --remove-indv T3553_LCA --remove-indv N919_LCA --remove-indv N997_LCA --remove-indv MVZ207321_LCA --remove-indv MVZ3720_LCA --remove-indv T3592_LCA --remove-indv N1020_LCA --remove-indv T3595_LCA --remove-indv T3585_LCA --remove-indv MG66_LCA --remove-indv JHT4108_LCA --remove-indv JHT2139_LCA --remove-indv JHT4125_LCA --remove-indv KEN589_LCA --remove-indv JHT4104_LCA --recode --out new_CENTAM_relaxed

# For PACMX, 29 inds removed because of high MD, leaving relaxed (n=245)
vcftools --vcf new_PACMX.vcf --remove-indv T3662_PAC --remove-indv V2617_pap_PAC --remove-indv V2620_pap_PAC --remove-indv T3605_Llano_PAC --remove-indv V2621_pap_PAC --remove-indv V2656_pap_PAC --remove-indv T416_RTex_CMX --remove-indv V2618_pap_PAC --remove-indv TF8677_Rspe_CMX --remove-indv TF8679_Rspe_CMX --remove-indv V2613_pap_PAC --remove-indv V2611_pap_PAC --remove-indv T3691_PAC --remove-indv T14440_Rmag_CMX --remove-indv V2782_PAC --remove-indv V2612_pap_PAC --remove-indv T1962_PAC --remove-indv T3676_PAC --remove-indv V10_PAC --remove-indv V2709_Rfor_PAC --remove-indv V2830_PAC --remove-indv T14429_Aten_PAC --remove-indv V2644_pap_PAC --remove-indv T3624_PAC --remove-indv T24978_PAC --remove-indv V9_PAC --remove-indv V2751_PAC --remove-indv T15_PAC --remove-indv V2762_PAC  --recode --out new_PACMX_relaxed

# For forreri, XX
# Make min10K vcf:
vcftools --vcf rana_n-1.vcf --remove-indv T14179_Rber_CMX --remove-indv T2757_LCA --remove-indv T3548_LCA --remove-indv T3594_LCA --remove-indv T3662_PAC --remove-indv TF8598_Rber_CMX --remove-indv V2617_pap_PAC --remove-indv V2620_pap_PAC --remove-indv N919_LCA --remove-indv T3553_LCA --remove-indv T3605_Llano_PAC --remove-indv V2621_pap_PAC --remove-indv V2656_pap_PAC --remove-indv V2618_pap_PAC --remove-indv TF7106_Rsie_RJO --remove-indv T3935_Rspe_CMX --remove-indv N997_LCA --remove-indv TF8595_Rber_CMX --remove-indv MVZ3720_LCA --remove-indv T3592_LCA --remove-indv TF8563_Rber_CMX --remove-indv T416_RTex_CMX --remove-indv MVZ207321_LCA --remove-indv TF8677_Rspe_CMX --remove-indv T3595_LCA --remove-indv V2611_pap_PAC --remove-indv T130_Rspe_CMX --remove-indv T3585_LCA --remove-indv TF8679_Rspe_CMX --remove-indv N1020_LCA --remove-indv V2782_PAC --remove-indv T3691_PAC --remove-indv V2613_pap_PAC --remove-indv T14440_Rmag_CMX --remove-indv T1962_PAC --remove-indv MG66_LCA --remove-indv T3676_PAC --remove-indv V2830_PAC --remove-indv JHT2139_LCA --remove-indv V10_PAC --remove-indv JHT4108_LCA --remove-indv V2612_pap_PAC --remove-indv V2644_pap_PAC --remove-indv V2709_Rfor_PAC --remove-indv T14429_Aten_PAC --remove-indv T2092_xochi_CMX --remove-indv JHT4125_LCA --remove-indv T3624_PAC --remove-indv KEN589_LCA --remove-indv T24978_PAC --remove-indv V9_PAC --remove-indv V2751_PAC --remove-indv JHT4104_LCA --remove-indv V2762_PAC --remove-indv T15_PAC --remove-indv MVZ145484_LCA --remove-indv JHT3827_LCA --remove-indv JHT2761_LCA --remove-indv AH197_LCA --remove-indv MVZ149013_LCA --remove-indv JHT3719_LCA --remove-indv MVZ263872_LCA --remove-indv T3681_CMX --remove-indv T3807_CMX --remove-indv T3695_PAC --remove-indv V2786_PAC --remove-indv JHT4105_LCA --remove-indv T3593_LCA --remove-indv JHT4118_LCA --remove-indv TF8604_Rjoh_RJO --remove-indv TF7111_Rzwe_RJO --remove-indv TF8510_Rspe_CMX --remove-indv T14193_Rber_CMX --remove-indv IRL57_LCA --remove-indv T3791_CMX --remove-indv TF7114_Rzwe_RJO --remove-indv TF7107_Rsie_RJO --recode --out min10K

#	2.	Remove any sites that are more than 75% missing (for --max-missing, 0 is totally missing, 1 is not missing)
# kept 707436 out of a possible 718861 Sites
vcftools --vcf ATL_MXPL_relaxed.recode.vcf --max-missing 0.25 --recode --out ATL_MXPL_relaxed_0.25miss

# kept 434394 out of a possible 441907 Sites
vcftools --vcf new_CENTAM_relaxed.recode.vcf --max-missing 0.25 --recode --out new_CENTAM_relaxed_0.25miss

# kept 767517 out of a possible 775665 Sites
vcftools --vcf new_PACMX_relaxed.recode.vcf --max-missing 0.25 --recode --out new_PACMX_relaxed_0.25miss

# kept 546944 out of a possible 564916 Sites
TODO add forreri

#	3.	Use Plink 1.9 or vcftools to generate a variant-based missing data report (.lmiss):
plink --vcf ATL_MXPL_relaxed_0.25miss.recode.vcf --missing --const-fid --allow-extra-chr --out ATL_MXPL_relaxed_0.25miss
plink --vcf new_CENTAM_relaxed_0.25miss.recode.vcf --missing --const-fid --allow-extra-chr --out new_CENTAM_relaxed_0.25miss
plink --vcf new_PACMX_relaxed_0.25miss.recode.vcf --missing --const-fid --allow-extra-chr --out new_PACMX_relaxed_0.25miss

vcftools --vcf forreri_0.25miss.recode.vcf --missing-site --out forreri_0.25miss

#	4.	Identify which SNPs should be retained after LD-pruning using the LD-pruning.R script.

#	5.	Now, remove selected SNPs from the matrix using vcftools.
# For ATL_MXPL, LD-pruning kept 50228 out of a possible 707436 Sites
vcftools --vcf ATL_MXPL_relaxed_0.25miss.recode.vcf --positions ATL_MXPL_relaxed_0.25miss_ldpruned.txt --recode --out ATL_MXPL_relaxed_0.25miss_ldp

# For CENTAM, LD-pruning kept 36126 out of a possible 434394 Sites
vcftools --vcf new_CENTAM_relaxed_0.25miss.recode.vcf --positions new_CENTAM_relaxed_0.25miss_ldpruned.txt --recode --out new_CENTAM_relaxed_0.25miss_ldp

# For PACMX, LD-pruning kept 45130 out of a possible 767517 Sites
vcftools --vcf new_PACMX_relaxed_0.25miss.recode.vcf --positions new_PACMX_relaxed_0.25miss_ldpruned.txt --recode --out new_PACMX_relaxed_0.25miss_ldp

# For forreri, kept 65907 out of a possible 546944 Sites
vcftools --vcf forreri_0.25miss.recode.vcf --positions forreri_0.25miss_ldpruned.txt --recode --out forreri_0.25miss_ldp

#	6.	Make a .bed file (and associated bim and fam files) for ADMIXTURE using Plink
plink --vcf ATL_MXPL_relaxed_0.25miss_ldp.recode.vcf --out ATL_MXPL_relaxed_0.25miss_ldp --allow-extra-chr --make-bed --const-fid

#	7.	Change chromosome names (the first column of the bim file needs to be the chrom name,
#	which has a specified list of names of which none apply to our data). The following just
#	replaces the RAD tag numbers of the first column with 0s.

awk '{$1=0;print $0}' ATL_MXPL_relaxed_0.25miss_ldp.bim > ATL_MXPL_relaxed_0.25miss_ldp.bim.tmp
mv ATL_MXPL_relaxed_0.25miss_ldp.bim.tmp ATL_MXPL_relaxed_0.25miss_ldp.bim

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

### Rename filenames so they all have replicates appended

for f in *; do mv -v "$f" "${f%.*}.rep1.${f##*.}"; done

# ========================================================================================
# ================================== LANDGEN ON FORRERI ==================================

# ========= fEEMS

# The following code is largely adapted by code provided by Erik Enbody

# fEEMS requires four input files:
#			1. Plink files (bim, bed, fam) which we'll generate from a vcf file
#			2. Coordinates file that matches the vcf file (ENSURE ORDERING IS THE SAME!)
#			3. Custom grid file (discrete global grid) for your geographic region
#			4. Outer coordinates (a sequence of vertices that outline a closed polygon)
# Important notes: this will impute to the mean

# Remove one sample that has no coordinates
vcftools --vcf forreri_0.25miss_ldp.recode.vcf --remove-indv V2797_PAC --recode --out forreri_0.25miss_ldp_n103

# Prune more for missing data and ensure there aren't invariant sites
bgzip -c forreri_0.25miss_ldp_n103.recode.vcf > forreri_0.25miss_ldp_n103.recode.vcf.gz
tabix -p vcf forreri_0.25miss_ldp_n103.recode.vcf.gz
bcftools view -e 'AF==1 | AF==0 | AF<0.01 | ALT="*" | F_MISSING > 0.50' -O v -o forreri_FILT.vcf forreri_0.25miss_ldp_n103.recode.vcf.gz

# Create bim, bed, and fam files from your vcf
plink --vcf forreri_FILT.vcf --out forreri_FILT --allow-extra-chr --autosome-num 95 --const-fid --make-bed

# Copy feems.txt from here: https://github.com/NovembreLab/feems/issues/15
# Then do:
mamba create --name feems_new --file feems.txt python=3.8.3
mamba activate feems_new
git clone https://github.com/NovembreLab/feems
cd feems
pip install .

# 1.	Run first two steps of fEEMS.R script to generate input data files.
# 2.	Run run_feems.py which will generate the spatial graph and also run the cross-validation
#		analysis.
# 3.	Run rest of FEEMS.R script to generate figures.

# ========================================================================================
# ======================================== HHSD ==========================================

sudo snap install cmake --classic # doing sudo apt-get install cmake installs v3.16.3

sudo apt-get remove cmake
export PATH=/home/wanglab/.local/bin:$PATH
sudo apt-get install build-essential
wget https://github.com/Kitware/CMake/releases/download/v3.29.5/cmake-3.29.5-linux-x86_64.tar.gz
tar xf cmake-3.29.5-linux-x86_64.tar.gz
cd cmake-3.29.5-linux-x86_64
./configure
make

# Install HHSD v0.9.9
git clone https://github.com/abacus-gene/hhsd
cd hhsd
pip install .

# Having dependency issues now
# ERROR: Package 'hhsd' requires a different Python: 3.12.4 not in '<3.10,>=3.9'

mamba create -c conda-forge -c bioconda -n hhsd3.9
mamba activate hhsd3.9
mamba install python=3.9.19

# GdiParameterError: 'gdi_threshold' incorrectly specified as '<0.3'.

# We are going to run HHSD on four datasets that correspond to three groupings within our 
# phylogenetic tree that have unclear species boundaries. In all cases, individuals were
# assigned populations based on the majority cluster for their ancestry coefficient values.
# Guide trees were based on our phylogeny. The three groups are as follows:
#			(1) The foothills group (magnaocularis, yavapaiensis, omiltemana OA, omiltemana
#				GE, Atenquique short, and Atenquique long) n=79
#			(2) The northern Atlantic coast and Central Mexican Plateau (berlandieri,
#				neovolcanica, tlaloci, spectabilis, and chichicuahutla) n=116
#			(3) Pacific coastal lowlands, i.e., forreri complex (forreri+cora+adleri, floresi, 
#				hillisi, miadis, Arcelia) n=104

### GENERATE INPUT FILES

#		We need to vastly reduce the sequence data for computational time's sake, and we
#		need the full locus (rather than just the SNPs) for each SNP that we want to include
#		in this analysis. 

#	1.	First, reduce the sequence data for computational time's sake. Remove any sites with
#		missing data above XXXX%:

#		Let's remove any sites with missing data above 25%:
vcftools --vcf forreri_0.25miss_ldp.recode.vcf --max-missing 0.75 --recode --out forreri_0.75miss_ldp # kept 848 out of a possible 65000 sites
vcftools --vcf ATL_MXPL_relaxed_0.25miss_ldp.recode.vcf --max-missing 0.7 --recode --out ATL_MXPL_0.7miss_ldp # kept 633 out of a possible 50228 sites
vcftools --vcf new_PACMX_relaxed_0.25miss_ldp.recode.vcf --max-missing 0.75 --recode --out PACMX_0.75miss_ldp # kept 501 out of a possible 45130 sites

#		Let's make larger datasets by increasing allowable missing data, as HHSD performs
#		better with more data provided by removing sites with missing data above 40%:
vcftools --vcf forreri_0.25miss_ldp.recode.vcf --max-missing 0.6 --recode --out forreri_0.6miss_ldp # kept 12453 out of a possible 65000 sites
vcftools --vcf ATL_MXPL_relaxed_0.25miss_ldp.recode.vcf --max-missing 0.6 --recode --out ATL_MXPL_0.6miss_ldp # kept 4684 out of a possible 50228 sites; 0.5 keeps 14332 / 50228 sites
vcftools --vcf new_PACMX_relaxed_0.25miss_ldp.recode.vcf --max-missing 0.6 --recode --out PACMX_0.6miss_ldp # kept 8226 out of a possible 45130 sites; 0.5 keeps 17616 / 45130 sites

#		For foothills and mxpl, we need to generate new vcfs with subsets 
#		containing relevant individuals. To do so:
vcftools --vcf ATL_MXPL_0.6miss_ldp.recode.vcf --keep mxpl_indv.txt --recode --out mxpl_06 # kept 116 inds
vcftools --vcf PACMX_0.6miss_ldp.recode.vcf --keep foothills_indv.txt --recode --out foothills_06 # kept 79 inds

#	2.	We need our phylip input file to have the full locus (rather than just the SNPs).
#		So, we'll get the RAD tag numbers and positions so we can retrieve the full loci
#		from the .loci file.
#			(a) Get a missing data report on the recoded vcf which we'll use to retrieve
#				locus and SNP names:
plink --vcf forreri_0.6miss_ldp.recode.vcf --missing --const-fid --allow-extra-chr --out forreri_0.6miss
plink --vcf foothills_06.recode.vcf --missing --const-fid --allow-extra-chr --out foothills_06
plink --vcf mxpl_06.recode.vcf --missing --const-fid --allow-extra-chr --out mxpl_06

#			(b) Using the lmiss file generated from Plink, run the beginning portion of the 
#				BPP.R script which will generate files named:
#					forreri_0.6miss_ldp_locusnames.txt
#					foothills_locusnames.txt & foothills_remove.txt & foothills_indv.txt
#					mxpl_locusnames.txt & mxpl_remove.txt & mxpl_indv.txt
			
#			(c) If you're using a subset of individuals (i.e., for the foothills or mxpl 
#				datasets), you need to select desired samples from the iPyrad *.loci file. 
#				To do so:
grep -vf foothills_remove.txt new_PACMX.loci > foothills.loci
grep -vf mxpl_remove.txt ATL_MXPL.loci > mxpl.loci

#			(c) Then, run the process_loci_file.py script which will generate separate txt 
#				files for each relevant locus from the iPyrad .loci file in the iPyrad 
#				output files directory:
mkdir loci_foothills
mkdir loci_mxpl
mkdir loci_forreri0.6

python process_loci_file.py loci_foothills/foothills_locus#.txt foothills_06_locusnames.txt foothills.loci
python process_loci_file.py loci_mxpl/mxpl_locus#.txt mxpl_06_locusnames.txt mxpl.loci
python process_loci_file.py loci_forreri0.6/forreri_0.6miss_ldp_locus#.txt forreri_0.6miss_ldp_locusnames.txt forreri.loci

#			(d)	Finally, run the following on all the files produced from process_loci_file.py 
#				which will remove the unnecessary locus labeling lines in each of the txt 
#				files:
for file in loci_forreri0.6/*.txt; do sed -i '' '/|/d' $file; done # all 104 samples

# If on linux (within each loci_* directory):
for file in *.txt; do sed -i '/|/d' $file; done

#			(e) Get the number of individuals (i.e., line numbers) and the number of sites
#				in each locus for the Phylip file. To do so, run:
for file in loci_forreri0.6/*.txt; do sed -i '' 's/^/^/' $file && chars=`head -n 1 $file | wc -m` && inds=`cat $file | wc -l` && tmp=`echo "$inds $(($chars-24))"` && sed -i '' "1s/^/$tmp\n/" $file && sed -i '' 's/^[ \t]*//' $file; done

# If on linux (N.B.: COUNT DIFFERENTLY FOR EACH OF THE DATASETS TO ENSURE PHYLIPS ARE CORRECT):
for file in loci_foothills/*.txt; do sed -i 's/^/^/' $file && chars=`head -n 1 $file | wc -m` && inds=`cat $file | wc -l` && tmp=`echo "$inds $(($chars-24))"` && sed -i "1s/^/$tmp\n/" $file && sed -i 's/^[ \t]*//' $file; done
for file in loci_mxpl/*.txt; do sed -i 's/^/^/' $file && chars=`head -n 1 $file | wc -m` && inds=`cat $file | wc -l` && tmp=`echo "$inds $(($chars-26))"` && sed -i "1s/^/$tmp\n/" $file && sed -i 's/^[ \t]*//' $file; done

#			(f) We're finally ready to make our Phylip file! Concatenate all the files 
#				together:
cat loci_forreri0.6/*.txt > forreri_0.6miss_ldp.phy
cat loci_foothills/foothills*.txt > foothills_06.phy
cat loci_mxpl/mxpl*.txt > mxpl_06.phy

#	3.	Make an Imap file which specifies individual assignment into populations/lineages
#		using the BPP.R script.
#			forreri-Imap.txt is for the ADMIXTURE results from K=5. Populations specified are:
#								forr=typotypic forreri, from N Mexico; flor=floresi;
#								hilli=hillisi; miad=samples southward of hillisi, including LCI;
#								arce=Arcelia form

#		Pops from foothills: 	omig (omiltemana guerrero), omio (omiltemana oaxaca), 
#								magn (magnaocularis), yava (yavapaiensis), 
#								ates (Atenquique short form), atel (Atenquique long form)

#	3. Control file

### RUNNING HHSD
mamba activate hhsd3.9

hhsd --cfile cf_mxpl_merge_0110.txt

# WARNING WARNING INCORRECT STATEMENT USED IN LITERAL EVAL