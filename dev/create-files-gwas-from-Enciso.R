#!/usr/bin/Rscript


args = commandArgs(trailingOnly = TRUE)
options (width=300)

args = c ("phenotype.tbl", "genotype.tbl", "genome.gff3")
phenotypeFile  = args [1]
genotypeFile   = args [2]
genomeFile     = args [3]


# Create genotypes
library(data.table)

genotypeAll    = read.table (file=genotypeFile, header=T)
genotypeBySNPs = transpose (genotypeAll)
colnames (genotypeBySNPs) = rownames (genotypeAll)
rownames (genotypeBySNPs) = colnames (genotypeAll)




break

# Create phenotypes for year 
phenotypeAll  = read.table (file=phenotypeFile, header=T)
for (year in c(2010:2015, 2017)) {
	print (year)
	phenoByYearAll = phenotypeAll [phenotypeAll$Year==year,]
	phenoByYear    = aggregate (phenoByYearAll$Score, list (phenoByYearAll$Genotype), mean)
	phenoByYear    = aggregate (phenoByYearAll$Score, list (phenoByYearAll$Genotype), mean)
	colnames (phenoByYear) = c ("Stub", "Score")

	outfile = sprintf ("phenotype-Enciso-LBlight-%s.tbl", year)
	write.table (phenoByYear, file=outfile, quote=F, row.names=F, sep=",")
	s = s + nrow (phenoByYearAll)
}

# Extract genome annotations 
genome    = as.matrix (read.table (file=genomeFile, header=F))
n         = nrow (genome)
genomeSnp = matrix (ncol=3, nrow=n) 
for (i in 1:nrow (genome)) {
	snp       = genome [i,]
	chr       = gsub ("chr(.+)", "\\1", snp[1])
	pos       = snp [4]
	name      = gsub (".*Name=(.+);.*", "\\1", snp [9])
	genomeSnp [i,] = c (name, chr, pos)
}
colnames (genomeSnp) = c("Marker", "Chrom", "Position")
write.table (file="genome-markers.tbl", genomeSnp, quote=F, sep="\t", row.names=F)
