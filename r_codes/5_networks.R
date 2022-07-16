library(phyloseq)
#pruned the OTUs with the mean relative abundance less than 0.005%
physeqr = phyloseq::transform_sample_counts(meta_physeq, function(x) x / sum(x))
physeqrF = filter_taxa(physeqr, function(x) mean(x) < .005/100, TRUE)
rmtaxa = taxa_names(physeqrF)
alltaxa = taxa_names(meta_physeq)
myTaxa = alltaxa[!alltaxa %in% rmtaxa]
physeqaF <- prune_taxa(myTaxa,meta_physeq)
physeqaF

#divided the filtered OTU table into two groups (i.e., high-altitude vs. high-latitude)
pa_phylo_net <- subset_samples(physeqaF, Region == "Pan-Arctic")
pa_phylo_net <- prune_taxa(taxa_sums(pa_phylo_net) > 0, pa_phylo_net) 
tp_phylo_net <- subset_samples(physeqaF, Region == "Tibetan Plateau")
tp_phylo_net <- prune_taxa(taxa_sums(tp_phylo_net) > 0, tp_phylo_net) 

pa_phylo_net
tp_phylo_net

#randomly resample the samples from alpine thermokarst lakes 
set.seed(1234) #To ensure reproducibility, the seed was set as 1234
rand.sample <- sample(sample_data(tp_phylo_net)$Sample_Name, 118, replace = F)
tp_phylo_rand_net <- subset_samples(tp_phylo_net, Sample_Name %in% rand.sample)
tp_phylo_rand_net <- prune_taxa(taxa_sums(tp_phylo_rand_net) > 0, tp_phylo_rand_net)
tp_phylo_rand_net

comm.pa.net <- otu_table(pa_phylo_net)
comm.tp.net <- otu_table(tp_phylo_rand_net)

#prepare the format of input data for Molecular Ecological Network Analyses (MENA) pipeline (http://ieg4.rccc.ou.edu/mena)
comm.pa.net[comm.pa.net == 0] <- ""
comm.tp.net[comm.tp.net == 0] <- ""


write.table(comm.pa.net, file = 'E:/thermokast_lakes/water_microbes/meta_analysis/results/tables/networks_6.30/comm_pa_net.txt')
write.table(comm.tp.net, file = 'E:/thermokast_lakes/water_microbes/meta_analysis/results/tables/networks_6.30/comm_tp_net.txt')
#then import the otu table into MENA pipeline

#To compare the difference in empirical network 
#indices between high-altitude and high-latitude 
#thermokarst lakes, the student t-test was employed 
#using the standard deviations derived from corresponding random networks
library(BSDA)
#Average clustering coefficient
tsum.test(mean.x = 0.643, s.x = 0.006, n.x = 118,
          mean.y = 0.653, s.y = 0.003, n.y = 118)

#Average path distance
tsum.test(mean.x = 2.692, s.x = 0.018, n.x = 118,
          mean.y = 5.080, s.y = 0.015, n.y = 118)

tsum.test(mean.x = 0.544, s.x = 0.003, n.x = 118,
          mean.y = 0.868, s.y = 0.004, n.y = 118)
