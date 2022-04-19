# loading packages
PKGs <- c('ggplot2', 'ape', 'Biostrings', 'vegan', 
          'reshape', 'cowplot', 'microbiome', 'RColorBrewer',
          'ggpubr', 'dplyr')

lapply(PKGs, require, character.only = TRUE, warn.conflicts = FALSE)

# set work dorectory
setwd('E:/thermokast_lakes/water_microbes/')
#conduct a phyloseq project
#read in metadata
metadata <- read.csv("./meta_analysis/data/meta_data/sample_data.csv",
                     header = T, row.names = 1)

#read in otu table
meta.otu.table <- read.csv("./meta_analysis/data/meta_data/meta_otu_table.csv",
                           header = T, row.names = 1, stringsAsFactors = F)
meta.otu.table <- as.matrix(meta.otu.table)

#read in taxonomy
meta.taxonomy <- read.csv("./meta_analysis/data/meta_data/meta_taxonomy.csv", sep=",", row.names=1)
meta.taxonomy <- as.matrix(meta.taxonomy)

# read in tree
meta.phy.tree <- read_tree("./meta_analysis/data/meta_data/meta_tree.nwk")

#read in represent dna sequences
meta.ref.seqs <- readDNAStringSet(file = "./meta_analysis/data/meta_data/meta_ref_seqs.fasta",
                                  format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
meta.otu.table <- otu_table(meta.otu.table, taxa_are_rows = TRUE)
meta.tax.table <- tax_table(meta.taxonomy)
meta.table <- sample_data(metadata)
meta_physeq <- phyloseq(meta.tax.table, meta.otu.table, meta.table, meta.phy.tree, meta.ref.seqs)
meta_physeq
#Convert to relative abundance
#meta_physeq_rel = phyloseq::transform_sample_counts(meta_physeq, function(x){x / sum(x)})
#phyloseq::otu_table(meta_physeq)[1:5, 1:5]
meta_physeq_rel <- microbiome::transform(meta_physeq, "compositional")
meta_table <- as(sample_data(meta_physeq), "data.frame")

pa_phylo <- subset_samples(meta_physeq, Region == "Pan-Arctic")
pa_phylo <- prune_taxa(taxa_sums(pa_phylo) > 0, pa_phylo) 
tp_phylo <- subset_samples(meta_physeq, Region == "Tibetan Plateau")
tp_phylo <- prune_taxa(taxa_sums(tp_phylo) > 0, tp_phylo) 
