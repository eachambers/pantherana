library(ape)
library(devtools)
library(here)

## This code does the following:
##     1. Removes potentially invariant sites from Phylip files for input into RAxML.

##    FILES REQUIRED:
##          rana_n-1.snps_min500.phy
##          rana_n-1.snps_min10K.phy
##          rana_n-1.snps_max80p.phy


# Load required functions -------------------------------------------------

#' Remove Invariant Sites
#'
#' https://github.com/bbanbury/phrynomics/blob/master/R/RemoveInvariantSites.R
#'
#' This function will determine if a site is variable (TRUE) or not (FALSE) and 
#' remove those sites that are not variable. Ambiguity codes are taken into account 
#' to mean either heterozygotes or uncertainty. For example, c("A", "A", "S") will 
#' return TRUE, because "S" reduces to "G" and "C" and so either of those is variable 
#' with "A"; c("A", "A", "M") will return FALSE, because M can be either A or C. If 
#' the "M" is uncertain and is an "A" then it is not variable.
#' 
#' @param SNPdataset SNP data in the class "matrix", "data.frame", or "snp"
#' @param chatty Optional print to screen messages
#' @export
#' @return Returns a subset dataset with only variable sites (i.e., SNPs).
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{IsVariable}
#' @examples
RemoveInvariantSites <- function(SNPdataset, chatty=FALSE){
  snpclass <- "table"
  if(class(SNPdataset) == "snp"){
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  snps <- sum(nchar(SNPdataset[1,])) # get no. of SNPs in original alignment
  initialLociLengths <- nchar(SNPdataset[1,])
  splitdata <- SplitSNP(SNPdataset) # split each SNP into a column; this takes a long time
  KeepVector <- apply(splitdata, 2, IsVariable) # determine which SNPs to keep depending on whether they're variable or not
  breaks <- which(splitdata[1,] == " ")
  newSNPdataset <- cSNP(splitdata, KeepVector=KeepVector, maintainLoci=TRUE) # make new SNP dataset with only variable sites
  newsnps <- sum(nchar(newSNPdataset[1,]))
  if(chatty)
    message(paste("removed", snps-newsnps, "of", snps, "sites"))
  if(snpclass == "snp")
    return(ReadSNP(newSNPdataset))
  else
    return(newSNPdataset)
}

#' Get lines to skip
#' For certain file formats, lines must be skipped (e.g., commented lines)
#'
#' @param file 
#'
#' @return
#' @export
#'
#' @examples
GetLinesToSkip <- function(file){
  a <- length(suppressWarnings(system(paste("grep 'ID:' ", file, sep=""), intern=TRUE)))
  b <- length(suppressWarnings(system(paste("grep ^[#+] ", file, sep=""), intern=TRUE)))
  c <- 0
  return(a+b+c)
}

#' Read in alignment file and convert to SNP file
#'
#' @param file input alignment file
#' @param row.names 
#' @param preprocess 
#' @param fileFormat specifies file format (options: "character", "phy", "nex", "data.frame")
#' @param extralinestoskip if the first line of file has something other than data, it must be skipped
#'
#' @return
#' @export
#'
#' @examples
ReadSNP <- function(file, row.names=1, preprocess=TRUE, fileFormat=NULL, extralinestoskip=0){
  if(class(file) == "character"){
    inputFileType <- fileFormat
    if(is.null(inputFileType))
      inputFileType <- FileFormat(file)
    if(inputFileType == "phy"){
      initializeTable <- read.table(file, row.names=row.names, skip=GetLinesToSkip(file)+extralinestoskip, stringsAsFactors=FALSE, colClasses=c("character"))
    }
    if(inputFileType == "nex")
      initializeTable <- ConvertNexDataToPhyData(read.nexus.data(file))
    if (is.na(inputFileType))
      stop("Having a hard time reading the file")
  }
  if(class(file) == "data.frame" || class(file) == "matrix"){
    initializeTable <- data.frame(lapply(file, as.character), stringsAsFactors=FALSE)
    rownames(initializeTable) <- rownames(file)
  }
  colnames(initializeTable) <- paste("locus", 1:ncol(initializeTable), sep="")
  if(preprocess){
    if(any(which(initializeTable == "_"))){
      underscores <- which(initializeTable == "_", arr.ind=T)
      initializeTable <- initializeTable[,-unique(underscores[,2])]
    }
  }
  ntax <- dim(initializeTable)[1]
  nloci <- dim(initializeTable)[2]
  nsites <- nchar(initializeTable[1,])
  snpdata <- list(data=as.data.frame(initializeTable), ntax=ntax, nloci=nloci, nsites=nsites)
  class(snpdata) <- "snp"
  return(snpdata)
}

#' Split each SNP into its own column
#'
#' @param SNPdataset 
#'
#' @return
#' @export
#'
#' @examples
SplitSNP <- function(SNPdataset){
  snpclass <- "table"
  if(class(SNPdataset) == "snp"){
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  loci <- dim(SNPdataset)[2]
  initialLociLengths <- nchar(SNPdataset[1,])
  splitSNP <- data.frame(matrix(nrow=dim(SNPdataset)[1], ncol=sum(initialLociLengths)+loci-1))
  for(j in sequence(dim(SNPdataset)[1])) {
    splitSNP[j,] <- strsplit(paste(SNPdataset[j,], collapse=" "), "")[[1]]
  }
  rownames(splitSNP) <- rownames(SNPdataset)
  if(snpclass == "snp")
    return(ReadSNP(splitSNP))
  else
    return(splitSNP)
}


cSNP <- function(splitSNP, KeepVector=NULL, maintainLoci=TRUE){
  snpclass <- "table"
  if(class(splitSNP) == "snp"){
    snpclass <- "snp"
    splitSNP <- splitSNP$data
  }
  if(is.null(KeepVector))
    KeepVector <- rep(TRUE, dim(splitSNP)[2])
  catSNP <- matrix(nrow=dim(splitSNP)[1], ncol=length(which(splitSNP[1,] == " "))+1)
  rownames(catSNP) <- rownames(splitSNP)
  for(j in sequence(dim(catSNP)[1])) {
    #totally annoying strsplit issue where it cleaves off one space at the end
    string <- strsplit(paste(splitSNP[j,][which(KeepVector=="TRUE")], collapse=""), " ")[[1]]
    if(length(string) +1 == dim(catSNP)[2]){
      string <- c(string, "")
    }
    catSNP[j,] <- string
  }
  if(any(which(catSNP[1,] == "")))
    catSNP <- catSNP[,-which(catSNP[1,] == "")]  #rewrite with lost loci
  if(!maintainLoci)
    catSNP <- apply(catSNP, 1, paste, collapse="")
  if(snpclass == "snp")
    return(ReadSNP(data.frame(catSNP, stringsAsFactors=FALSE)))
  else
    return(data.frame(catSNP, stringsAsFactors=FALSE))
}

#' Determine whether a site is variable or invariant
#'
#' @param SNP 
#'
#' @return
#' @export
#'
#' @examples
IsVariable <- function(SNP){
  var <- FALSE
  bases <- c("A", "C", "G", "T", "U")
  basesInSNP <- bases[which(bases %in% SNP)]
  if(all(SNP == " "))
    return(TRUE)
  if(length(basesInSNP) == 0)  #all "N"
    return(FALSE)
  if(length(basesInSNP) > 1)  #more than one base anyway, so skip checking ambigs
    return(TRUE)
  if(length(basesInSNP) == 1) { #if only one base plus ambigs
    for(i in 2:length(SNP)) {
      var <- c(var, !any(ReturnNucs(SNP[i]) %in% basesInSNP))
    }
  }
  return(any(var))
}

#' Returns nucleotide possibilities (resolves ambiguities)
#'
#' @param NucCode 
#' @param forSNAPP 
#'
#' @return
#' @export
#'
#' @examples
ReturnNucs <- function(NucCode, forSNAPP=FALSE) {
  possibilities <- NULL
  if(NucCode == "A" || NucCode == "G" || NucCode == "C" || NucCode == "T" || NucCode == "U")
    possibilities <- NucCode
  if(NucCode == "N" || NucCode == "-" || NucCode == "?") {
    if(forSNAPP)  possibilities <- "-"
    else  possibilities <- c("A", "G", "C", "T", "U")
  }
  if(NucCode == "R")  possibilities <- c("A", "G")
  if(NucCode == "Y")  possibilities <- c("C", "T")
  if(NucCode == "W")  possibilities <- c("A", "T")
  if(NucCode == "S")  possibilities <- c("G", "C")
  if(NucCode == "M")  possibilities <- c("A", "C")
  if(NucCode == "K")  possibilities <- c("G", "T")
  if(NucCode == "B")  possibilities <- c("G", "C", "T")
  if(NucCode == "H")  possibilities <- c("A", "C", "T")
  if(NucCode == "D")  possibilities <- c("A", "G", "T")
  if(NucCode == "V")  possibilities <- c("A", "G", "C")
  return(possibilities)
}

#' Write out SNP dataset with invariant sites removed
#'
#' @param SNPdataset output SNP dataset from RemoveInvariantSites
#' @param file file name to save
#' @param format file format (options: "phylip", "nexus")
#' @param missing how to code missing data (defaults to "N")
#'
#' @return
#' @export
#'
#' @examples
WriteSNP <- function(SNPdataset, file="", format="phylip", missing="N"){
  if(class(SNPdataset) == "snp")
    SNPdataset <- SNPdataset$data
  if(format == "phylip"){
    write(c(dim(SNPdataset)[1], sum(nchar(SNPdataset[1,]))), file=file)
    write.table(SNPdataset, file=file, quote=FALSE, append=TRUE, col.names=FALSE)
  }
  if(format == "nexus"){
    if(class(SNPdataset) == "data.frame" || class(SNPdataset) == "matrix")
      if(dim(SNPdataset)[2] > 1)
        SNPdataset <- data.frame(apply(SNPdataset, 1, paste, collapse=""))
    nchars <- min(apply(SNPdataset, 1, nchar))
    write(paste("#NEXUS"), file)
    write(paste("[Written ", Sys.Date(), " via phrynomics]", sep=""), file, append=TRUE)
    write(paste("BEGIN Data;"), file, append=TRUE)
    write(paste("   DIMENSIONS NTAX=", dim(SNPdataset)[1], " NCHAR=", nchars, ";", sep=""), file, append=TRUE)
    write(paste("   FORMAT DATATYPE=Standard INTERLEAVE=no missing=", missing, ";", sep=""), file, append=TRUE)
    write(paste("Matrix"), file, append=TRUE)
    write.table(SNPdataset, file, append=TRUE, quote=FALSE, col.names=FALSE)
    write(paste(""), file, append=TRUE)
    write(paste(";"), file, append=TRUE)
    write(paste("END;"), file, append=TRUE)
    write(paste(""), file, append=TRUE)
  }
}

# Run functions ------------------------------------------------------------

# ReadSNP will be run first, since we need to convert our nexus/phylip files into SNP ones
# Need to skip first line because it doesn't contain data
# This will take a little while
# MAKE SURE TO DO THE FOLLOWING ON ALL THREE FILES!
# dat <- ReadSNP(file = here("data", "rana_n-1.snps_min500.phy"), extralinestoskip = 1, fileFormat = "phy")
dat <- ReadSNP(file = here("data", "rana_n-1.snps_max80p.phy"), extralinestoskip = 1, fileFormat = "phy")

# Can check that ReadSNP worked:
dat$ntax # 595 / 414
dat$nsites # 350,292 / 350,267

# Remove invariant sites --------------------------------------------------

dat_invar <- RemoveInvariantSites(SNPdataset = dat, chatty = TRUE) # removed 103,616 / XXX of 350,292 sites

# Write out files ---------------------------------------------------------

# WriteSNP(dat_invar, file = here("data", "rana_n-1.snps_min500_invarrem.phy"), missing = "N")
# WriteSNP(dat_invar, file = here("data", "rana_n-1.snps_min10K_invarrem.phy"), missing = "N")
WriteSNP(dat_invar, file = here("data", "rana_n-1.snps_max80p_invarrem.phy"), missing = "N")
