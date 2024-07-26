#!/usr/bin/env Rscript

##
# Script to combine KO_map files
#
# This file requires a list 

library(openxlsx)

# This assumes that all the ko_map.txt files are in this directory
map_files <- list.files(pattern="\\.ko_map\\.txt$")

# Create a list to store the data from all the ko_map.txt files.
maps <- list()

# Load all the maps
for ( a_map in map_files ) { 
	sample <- sub("\\.ko_map\\.txt$", "", a_map)
	maps[[sample]] <- read.delim(a_map, row.names=1)
}

# Get a list of the kos from the first mapping file.  
# These should be the same for all the files.  
# In not, you would need to change this behavior
kos <- colnames(maps[[1]])

# Create an new excel workbook
wb <- createWorkbook()
k <- 0

# Loop over all the KOs
for ( ko in kos ) { 
	k <- k + 1
	# Merge the first two samples to create the initial table
	this_ko <- merge(maps[[1]][,ko, drop=F], maps[[2]][,ko, drop=F], by.x=0, by.y=0, all=T)
	colnames(this_ko) <- c("Taxon", names(maps)[1], names(maps)[2])
	# loop over the remaining samples to get the final table
	if ( length(maps) > 2 ) {
  	for ( m in 3:length(maps) ) { 
  		this_map <- maps[[m]][,ko,drop=F]
  		colnames(this_map) <- names(maps)[m]
  		this_ko <- merge(this_ko, this_map, by.x=1, by.y=0, all=T)	
  	}
	}
	# If there are any NA values set to zero (so that sum would work)
	this_ko[is.na(this_ko)] <- 0

	# Create a sum of the counts, except the first column as it is the taxonomic assignments
	this_ko$SUM <- rowSums(this_ko[,-1])

	# Remove any rows in which the SUM is zero
	this_ko <- subset(this_ko, SUM > 0)
	
	# Move the SUM column to after the taxon column
	this_ko <- this_ko[, c(1, ncol(this_ko), 2:(ncol(this_ko)-1))]

	# Create a worksheet tab for this KO and add the data
	addWorksheet(wb, ko)		
	writeDataTable(wb, k, this_ko, tableStyle="TableStyleLight1")
}

# Save the Excel worksheet, if you want to change the name of the file do so here.
saveWorkbook(wb, "ko_maps.xlsx", overwrite=T)
