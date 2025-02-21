# ================================ FEEMS pre-processing ==================================

# The following code is largely adapted by code provided by Erik Enbody

# FEEMS requires four input files:
#			1. Plink files (bim, bed, fam) which we'll generate from a vcf file (`forreri_FILT.bed`)
#			2. Coordinates file that matches the vcf file (ENSURE ORDERING IS THE SAME! `forreri_coords.txt`)
#			3. Custom grid file (discrete global grid) for your geographic region (`forr_grid.shp`)
#			4. Outer coordinates (a sequence of vertices that outline a closed polygon; `forreri_outer.txt`)
# Important notes: this will impute to the mean

#	1.	Remove one sample that has no coordinates
vcftools --vcf forreri_0.25miss_ldp.recode.vcf --remove-indv V2797_PAC --recode --out forreri_0.25miss_ldp_n103

#	2.	Prune more for missing data and ensure there aren't invariant sites
bgzip -c forreri_0.25miss_ldp_n103.recode.vcf > forreri_0.25miss_ldp_n103.recode.vcf.gz
tabix -p vcf forreri_0.25miss_ldp_n103.recode.vcf.gz
bcftools view -e 'AF==1 | AF==0 | AF<0.01 | ALT="*" | F_MISSING > 0.50' -O v -o forreri_FILT.vcf forreri_0.25miss_ldp_n103.recode.vcf.gz

#	3.	Create bim, bed, and fam files from your vcf
#			17789 variants and 103 samples pass filters and QC
plink --vcf forreri_FILT.vcf --out forreri_FILT --allow-extra-chr --autosome-num 95 --const-fid --make-bed

#	4. Set up FEEMS. Copy feems.txt from here: https://github.com/NovembreLab/feems/issues/15
#		and then do:
mamba create --name feems_new --file feems.txt python=3.8.3
mamba activate feems_new
git clone https://github.com/NovembreLab/feems
cd feems
pip install .

#	5. Now, you're ready to do the following:
#		(a) Run first two steps of `fEEMS.R` script to generate input data files.
#		(b) Run `run_feems.py` which will generate the spatial graph and also run the 
#			cross-validation analysis.
#		(c) Run the rest of the `fEEMS.R` script to visualize results and generate figures.