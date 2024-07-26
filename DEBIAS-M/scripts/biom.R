library(dplyr)
library(biomformat)

#file = system.file("extdata", "./BIOM/131235/New folder/feature-table.biom", package = "biomformat")

#OTU_reads <- read_biom("./BIOM2/173205/feature-table.biom")
#OTU_reads <- read_hdf5_biom("./BIOM/131235/New folder/feature-table.biom")

OTU_reads <- read_hdf5_biom("./BIOM2/173205/feature-table.biom")
OTU_reads= biom(OTU_reads)


#outfile = tempfile()
#write_biom(OTU_reads, outfile)
#y = read_biom(outfile)



OTU_table <- as.data.frame(as.matrix(biom_data(OTU_reads)))


#### General ####
OTU_reads <- read_hdf5_biom("./BIOM/80310/otu_table.biom")
OTU_reads = biom(OTU_reads)

outfile = tempfile()
write_biom(OTU_reads, outfile)
y = read_biom(outfile)

OTU_table <- as.data.frame(as.matrix(biom_data(y)))



