# =============== Getting basic data characteristics from pooled assembly ================

# 1.	Convert full phylip file (rana_n-1.phy) to Nexus for input into PAUP* (I did so manually).
#		Do the same for the .snps file for running pooled assembly in RAxML-ng (rana_n-1.snps.nexus).

# 2.	Import file into PAUP*: (this may take a while as the matrix size is large)
> execute rana_n-1.nexus

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