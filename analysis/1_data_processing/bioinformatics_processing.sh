# ================ Running iPyrad to obtain pooled assembly (all samples) ================

# iPyrad v.0.9.85

# There were four sequencing projects (with varying numbers of pools within each) from
# our ddRAD runs. For iPyrad steps 1-3, we ran each pool within each project separately; 
# there are params files for each pool in our supplementary data files which are named
# according to the project and pool (e.g., "JA19241_pool2_R1.txt"):
#	JA19241 (n=130); pools 1-6
#	JA19242 (n=311); pools 1-13
#	JA20247 (n=20); pool 1
#	JA20248 (n=163); pools 1-7

# We also incorporated several samples from Chambers et al. (2023). We first copied fastqs
# from their SRA archive and re-ran iPyrad steps 2-3 (n=9). This assembly is named "epirana"
# in our supplementary files; the params file is called "params-epirana.txt".

# Once the above steps were completed, we merged everything together using the following:
ipyrad -m all_rana params-JA19241_pool1_R1.txt params-JA19242_pool11_R1.txt params-JA19242_pool5_R1.txt  params-JA20248_pool2.txt params-JA19241_pool2_R1.txt params-JA19242_pool12_R1.txt params-JA19242_pool6_R1.txt params-JA20248_pool3.txt params-JA19241_pool3_R1.txt params-JA19242_pool13_R1.txt params-JA19242_pool7_R1.txt params-JA20248_pool4.txt params-JA19241_pool4_R1.txt params-JA19242_pool1_R1.txt params-JA19242_pool8_R1.txt params-JA20248_pool5.txt params-JA19241_pool5_R1.txt params-JA19242_pool2_R1.txt params-JA19242_pool9_R1.txt params-JA20248_pool6.txt params-JA19241_pool6_R1.txt params-JA19242_pool3_R1.txt params-JA20247_pool1_R1.txt params-JA20248_pool7.txt params-JA19242_pool10_R1.txt params-JA19242_pool4_R1.txt params-JA20248_pool1.txt params-epirana.txt

# We also had to branch off and remove TF8608_Rber_CMX because there were no clusters formed; 
# to do so we ran:
XXX rana_n-1 -b <======================

# ============== Running iPyrad on separate assemblies (subsets of samples) ==============

# 1.	We now want to split up samples based on the phylogenetic tree, building assemblies with
#		subsets of closely related samples (to maximize phylo info obtained). To do so, we moved 
#		individual sorted/demultiplexed fastqs over to a separate directory and made a new iPyrad 
#		project from there, specifying the sorted fastq path (see params files in supp materials).
#		To do so, we ran:

ipyrad -p params-ATL_MXPL -s 123456
ipyrad -p params-forreri -s 123456
ipyrad -p params-PACMX -s 123456
ipyrad -p params-CENTAM -s 123456

#	Please refer to Dataset S1 for which samples were assigned into each separate assembly.
#	The separate assemblies are as follows:
	
#	PACMX		Pacific coast (lowlands and foothills) of Mexico. Includes described species:
#				R. forreri (MX), R. macroglossa (MX), R. spectabilis (Oaxaca only),
#				R. omiltemana, R. magnaocularis, R. yavapaiensis
#				Also includes unnamed species (Hillis & Wilcox 2005): sp. 7 (Jalisco) and 
#				sp. 8 (Puebla), as well as the "papagayo form" described by Hillis et al. 
#				(1983) and Arcelia and Colima forms.
				
#				Total of 274 samples
	
#	CENTAM		Samples occurring in Central America. Includes described species:
#				R. lenca, R. miadis, R. macroglossa (CENTAM), R. forreri (CENTAM), R. taylori
#				Also includes unnamed species (Hillis & Wilcox 2005): sp. 4 (PA), 
#				sp. 5 (CR), sp. 6 (CR)
				
#				Total of 158 samples

#	ATL_MXPL	Atlantic coast and the Mexican Plateau. Includes described species:
#				R. taylori, R. macroglossa, R. brownorum [not explicitly IDed], 
#				R. berlandieri, R. spectabilis, R. tlaloci, R. neovolcanica, 
#				R. chichicuahutla [not explicitly IDed].
#				Also includes unnamed species (Hillis & Wilcox 2005): sp. 3 (which is 
#				R. spectabilis).
	
#				Total of 202 samples

#	forreri		All members from the Pacific coastal lowlands. Includes described species:
#				R. miadis, R. forreri, R. cora, R. adleri, R. hillisi, R. floresi
				
#				Total of 104 samples

# 2.	For each of these assemblies, we generated missing data proportions in PAUP* using 
#		a Nexus file for the SNP matrix (see `basic_data_characteristics.sh`). We also 
#		processed missing data and calculated average read depth using the `basic_stats.R` 
#		script.