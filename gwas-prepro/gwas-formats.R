#!/usr/bin/Rscript

options (width=300, stringsAsFactors=F)
args = commandArgs(trailingOnly=T)
library (parallel)

#------------------------------------------------------------------------------
## Format gwaspoly phenotype to plink format
#------------------------------------------------------------------------------
gwaspToPlinkPhenotype <- function (gwaspPhenoFile) 
{
	message ("Creating plink phenotype...")
	phenotype = read.csv (file=gwaspPhenoFile, header=T)
	idNames = as.character (phenotype [,1])

	samples = phenotype [,1]
	traits  = phenotype [,2]

	message (">>> Writing plink phenotype...")
	plinkPheno = cbind (FID=0,IID=samples, TRAIT= traits)
	outFile = paste0 (strsplit (gwaspPhenoFile, split="[.]")[[1]][1], "-plink.tbl")
	write.table (file=outFile, plinkPheno, col.names=T, row.names=F, quote=F, sep="\t")
}

#------------------------------------------------------------------------------
## Create plink MAP file from gwaspoly genotype 
#------------------------------------------------------------------------------
gwaspToPlinkMap <- function (gwaspGenoFile) 
{
	message (">>> Creating plink MAP file...")
	genotype    = read.csv (file=gwaspGenoFile, header=T,stringsAsFactors=F)
	markers     <- genotype [,1]
	chromosomes <- genotype [,2]
	positions   <- genotype [,3]

	plinkMap  = cbind (chr=chromosomes, iid=markers, dist=0, positions=positions)
	outFile   = paste0 (strsplit (gwaspGenoFile, split="[.]")[[1]][1], "-plink.map")
	write.table (file=outFile, plinkMap, col.names=F, row.names=F, quote=F, sep="\t")
}

#----------------------------------------------------------
## Create plink PED file from gwaspoly genotype 
#----------------------------------------------------------
gwaspToPlinkPed <- function (gwaspGenoFile) {
	message (">>> Creating plink PED file...")
	genotype    = read.csv (file=gwaspGenoFile, header=T)
	#namesGeno  = c("fid", "iid", "pid", "mid", "sex", "phe", rep ("X", ncol (talleles)))
	alleles    <- genotype [,-c(1,2,3)]

	message (">>> Creating transposed genotype...")
	markersIds        <- genotype [,1] 
	samplesIds        <- colnames (alleles)
	transposedAlleles <- t(alleles)
	rownames (transposedAlleles) = samplesIds
	colnames (transposedAlleles) = markersIds

	outFile   = paste0 (strsplit (gwaspGenoFile, split="[.]")[[1]][1], "-plink-transposed.tbl")
	write.table (file =outFile, transposedAlleles, col.names=T, row.names=T, quote=F, sep="\t")

	message (">>> Creating PED genotype...")
	plinkAlleles   = formatToPlinkAlleles (transposedAlleles)

	colnames (plinkAlleles) = markersIds
	genoPED    = cbind (0, samplesIds, 0,0,0,-9, plinkAlleles)

	outFile   = paste0 (strsplit (gwaspGenoFile, split="[.]")[[1]][1], "-plink.ped")
	write.table (file=outFile, genoPED, col.names=F, row.names=F, quote=F, sep="\t")

}

#----------------------------------------------------------
# Add tabs to alleels changign AG --> A	G
#----------------------------------------------------------
formatToPlinkAlleles <- function (transposedAlleles) {
	#alleles  [is.na (alleles)] = "00"
	message (">>> Formatting alleles to plink...")
	ncols = ncol (transposedAlleles)
	nrows = nrow (transposedAlleles)
	#alleles = t (alleles)
	tabs <- function (x) {return (sprintf("%s %s %s %s", substr(x,1,1),substr(x,2,2),substr(x,3,3),substr(x,4,4)))}
	tabbedAlleles = matrix (mclapply (transposedAlleles,tabs,mc.cores=4), ncol=ncols, nrow=nrows, byrow=F)
	return (tabbedAlleles)
}

#----------------------------------------------------------
# Main
#----------------------------------------------------------
args = c ("agrosavia-genotype-tetra-ACGT.tbl", "agrosavia-phenotype.tbl")

gwaspGenoFile  = args [1]
gwaspPhenoFile = args [2]

message (">>> Converting gwaspoly to plink formats...")
gwaspToPlinkPhenotype (gwaspPhenoFile)
gwaspToPlinkMap (gwaspGenoFile)
gwaspToPlinkPed (gwaspGenoFile)

	
