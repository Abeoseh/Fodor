library(dplyr)
library(biomformat)
library(dada2)

#### General ####
OTU_reads <- read_hdf5_biom("173204_all.biom") #"C:/Users/brean/Downloads/masters/Fodor/DEBIAS-M/15016/173204_all.biom"


OTU_reads = biom(OTU_reads)

OTU_table <- as.data.frame(as.matrix(biom_data(OTU_reads)))

seqs <- row.names(OTU_table)


taxa <- assignTaxonomy(seqs, "C:/Users/brean/Downloads/masters/Fodor/Taxonomy_Dada2/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa_species <- addSpecies(taxa, "C:/Users/brean/Downloads/masters/Fodor/Taxonomy_Dada2/silva_species_assignment_v138.1.fa.gz")


df <- as.data.frame(taxa_species)

sum(is.na(df$Genus))

