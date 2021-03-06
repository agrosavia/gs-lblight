#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE)
#----------------------------------------------------------
# Functions for extracting samples and markers names
#----------------------------------------------------------
main <- function () {
	# Extract and sort the names with integer suffix from a table
	#args = "ClusterCall_prediction_CCC-fixed.csv"
	#args = "berdugo-genotype-s003-sample100.csv"
	print (args)
	inFile    = args [1]
	filename = strsplit (inFile, "[.]")[[1]][1] 
	print (sprintf ("\nExtracting names from %s to %s", inFile, filename))

	# Load and get names
	genoB.all = read.table (inFile, header=T, sep=",", row.names=1)
	dim (genoB.all)

	alleles = getAllelesTable (genoB.all, filename)
	getSampleNamesBerdugo (genoB.all, filename)
	getMarkerNames (genoB.all, filename)
}

#----------------------------------------------------------
# Filter extracting the Alleles (And_XXX)
#----------------------------------------------------------
getAllelesTable <- function (genoB.all, filename) {
	nms = names (genoB.all)
	genoB.alleles = genoB.all [,nms [grep ("^And[^.]*$", nms)]]
	outFile = sprintf ("%s-alleles.tbl", filename)
	write.table (file=outFile, genoB.alleles, quote=F)

	return (genoB.alleles)
}
#----------------------------------------------------------
# Get sample names from genotype
#----------------------------------------------------------
getSampleNamesBerdugo <- function (genoB.all, filename) {
	nms = names (genoB.all)

	# Filter to specific columns
	nms.andall = nms [grep ("And", nms)]
	nms.and = nms.andall [-grep ("[.]", nms.andall)]
	#nms.and = nms.andall

	# Create ids from suffixes
	la = unlist (strsplit(nms.and, "_"))
	nms.ids = sprintf ("%03d", strtoi (la [la != "And"]))

	# Create matrix ids + names, sort, and write
	m =  matrix(c(nms.ids,nms.and), ncol=2)
	markers.sorted = m[order (m[,1]),]
	write.table (file=sprintf ("%s-samples-names-ids.tbl", filename), markers.sorted,row.names=F, quote=F,col.names=F)
	write.table (file=sprintf ("%s-samples-names.tbl", filename), markers.sorted[,2],row.names=F, quote=F,col.names=F)
}

#----------------------------------------------------------
# Get Marker names from genotype
#----------------------------------------------------------
getMarkerNames <- function  (genoB.all, filename) {
	# Get, split, create ids
	markers = rownames (genoB.all)
	mc1     = markers [grep ("c1", markers)]
	mc2     = markers [grep ("c2", markers)]

	mls_c1  = unlist (strsplit (mc1, "c1_"))
	mls_c2  = unlist (strsplit (mc2, "c2_"))

	#mids = sprintf ("%05d", strtoi (mls [mls !="solcap"]))
	mids_c1 = sprintf ("c1%05d", strtoi (mls_c1[-grep ("solcap", mls_c1)]))
	mids_c2 = sprintf ("c2%05d", strtoi (mls_c2[-grep ("solcap", mls_c2)]))
	message (length(mids_c1))
	message (length(mids_c2))

	mcs  = c (mc1, mc2)
	mids = c (mids_c1, mids_c2)

	mmat  = matrix (c (mids, mcs), ncol=2)
	msort = mmat [order (mmat [,1]),]

	write.table (file=sprintf ("%s-markers-ids.tbl", filename), msort,row.names=F, quote=F,col.names=F)
	write.table (file=sprintf ("%s-markers.tbl", filename), msort[,2],row.names=F, quote=F,col.names=F)
}

#----------------------------------------------------------
# Call main
#----------------------------------------------------------
main ()

