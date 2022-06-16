#PERMANOVA analysis and NMDS plot
meta.ord <- ordinate(meta_physeq, "NMDS", "bray")
meta.ord
NMDS_plot <- plot_ordination(meta_physeq, meta.ord, type="samples", color="Region") +
  geom_point(size = 2.5) + 
  scale_color_manual(values = c('#d95f02', '#1b9e77')) +
  #stat_ellipse(type = "norm", linetype = 1) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill = Region)) +
  theme_bw() +
  theme(legend.position = c(0.85,0.88),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour='black',size=12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key = element_blank(),
        legend.background = element_rect(colour = "white"))
metadata <- as(sample_data(meta_physeq), "data.frame")
adonis2(phyloseq::distance(meta_physeq, method = "bray") ~ Region,
        data = metadata)

#diff class using microeco package
library(microeco)
meco_df <- phyloseq2meco(meta_physeq)
#calculate the abundance table
m1 <- meco_df$cal_abund()
#Lefse
ps_genus <- phyloseq::tax_glom(meta_physeq, taxrank = 'Genus')
meco_genus_df <- phyloseq2meco(ps_genus)
#calculate the abundance table
m1_genus <- meco_genus_df$cal_abund()

m1_genus <- trans_diff$new(dataset = meco_genus_df, method = "lefse", 
                           group = "Region", alpha = 0.01, 
                           lefse_subgroup = NULL)
# t1$res_lefse is the LEfSe result
# t1$res_abund is the abundance information
m1_genus$plot_lefse_bar(LDA_score = 4, color_values = c('#d95f02', '#1b9e77'), width = 0.5)
m1_genus$plot_diff_abund(use_number = 1:30, color_values = c('#d95f02', '#1b9e77'))

# we can format the data for Lefse analysis (http://huttenhower.sph.harvard.edu/galaxy)
ps_genus_sel <- subset_taxa(ps_genus, Family == "Sporichthyaceae" |
                              Genus == "Candidatus Aquiluna" |
                              Genus == "Candidatus Limnoluna" |
                              Genus == "Candidatus Planktoluna" |
                              Genus == "Sediminibacterium" |
                              Genus == "Algoriphagus" |
                              Genus == "Emticicia" |
                              Genus == "Flavobacterium" |
                              Genus == "Pedobacter" |
                              Family == "NS11_12marinegroup" |
                              Class == "Firmicutes" |
                              Family == "Caulobacteraceae" |
                              Genus == "Rhodoblastus" |
                              Genus == "Roseiarcus" |
                              Genus == "Sphingorhabdus" |
                              Genus == "Polynucleobacter" |
                              Genus == "Variovorax" |
                              Family == "Rubritaleaceae" |
                              Family == "Verrucomicrobiaceae" |
                              Genus == "uncultured Opitutaebacterium")

phyloseq2lefse <- function(
    ps,
    covars,
    file.name = "lefse_data.txt",
    taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
    transpose.otus = TRUE
) {
  
  # grab the taxa table from the phyloseq object and coerce into a matrix
  tax.tbl <- as(phyloseq::tax_table(ps), "matrix")
  tax.tbl <- tax.tbl[, taxa.levels]
  tax.tbl.join <- as.data.frame(do.call(paste, c(as.data.frame(tax.tbl), sep = "|")))
  
  row.names(tax.tbl.join) <- row.names(tax.tbl)
  colnames(tax.tbl.join) <- NULL
  
  # grab the otu table from the phyloseq object
  otu.tbl <- as(phyloseq::otu_table(ps), "matrix")
  colnames(otu.tbl) <- NULL
  tax.otu.tbl <- cbind(tax.tbl.join, otu.tbl)
  
  # grab the sample table from the phyloseq object
  smpl.data <- as(phyloseq::sample_data(ps), "data.frame")
  smpl.data$Sample <- row.names(smpl.data)
  t.smpl.data <- t(smpl.data)
  t.smpl.data <- as.matrix(t.smpl.data[c("Sample", covars), ])
  t.smpl.data <- cbind(rownames(t.smpl.data), t.smpl.data)
  colnames(t.smpl.data) <- colnames(tax.otu.tbl)
  
  final.data <- rbind(t.smpl.data, tax.otu.tbl)
  write.table(final.data, file = file.name, row.names = F, col.names = F, sep = "\t", quote = FALSE)
}

phyloseq2lefse(
  ps = ps_genus_sel,
  covars = 'Region',
  file.name = "./meta_analysis/results/tables/lefse_data_sel_genus.txt",
  taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  transpose.otus = TRUE
)


# alpha diversity
diversity <- estimate_richness(meta_physeq, measures = c("Chao1", 'Shannon', 'Simpson'))
meta_diversity <- cbind(diversity, Region = meta_table$Region, 
                        Sitegroup = meta_table$Sitegroup,
                        Site = meta_table$Site, MAT = metadata$MAT,
                        MAP = metadata$MAP, DOC = metadata$DOC,
                        SUVA254 = metadata$SUVA254, a320 = metadata$a320,
                        pH = metadata$pH)
meta_diversity %>% dplyr::select(c(1,3,4,5)) %>%
  group_by(Region) %>%
  dplyr::summarise_all(list(mean = mean, sd = sd, se = ~sd(./sqrt(.))))

# test the difference
par(mfrow = c(1, 3))
hist(log(meta_diversity$Chao1))
hist(meta_diversity$Shannon)
hist(log(meta_diversity$Simpson+1))
par(mfrow = c(1, 1))
library(lme4)
library(lmerTest)
library(multcomp)
mode1 <- lmer(log(Chao1) ~ Region + (1|Site), meta_diversity)
summary(mode1)
MuMIn::r.squaredGLMM(mode1)
mode2 <- lmer(Shannon ~ Region + (1|Site), meta_diversity)
summary(mode2)
MuMIn::r.squaredGLMM(mode2)
mode3 <- lmer(log(Simpson) ~ Region + (1|Site), meta_diversity)
summary(mode3)
MuMIn::r.squaredGLMM(mode3)
ggplot(meta_diversity, aes(x=MAT, y=Simpson)) + geom_point(alpha = 0.3)
mode4 <- lmer(log(Chao1) ~ MAT + (1|Site), meta_diversity)
summary(mode4)
MuMIn::r.squaredGLMM(mode4)
mode5 <- lmer(Shannon ~ MAT + (1|Site), meta_diversity)
summary(mode5)
MuMIn::r.squaredGLMM(mode5)
mode6 <- lmer(log(Simpson) ~ MAT + (1|Site), meta_diversity)
summary(mode6)
MuMIn::r.squaredGLMM(mode6)
ggplot(meta_diversity, aes(x=MAP, y=Chao1)) + geom_point(alpha = 0.3)
mode7 <- lmer(log(Chao1) ~ MAP + (1|Site), meta_diversity)
summary(mode7)
MuMIn::r.squaredGLMM(mode7)
mode8 <- lmer(Shannon ~ MAP + (1|Site), meta_diversity)
summary(mode8)
MuMIn::r.squaredGLMM(mode8)
mode9 <- lmer(log(Simpson) ~ MAP + (1|Site), meta_diversity)
summary(mode9)
MuMIn::r.squaredGLMM(mode9)
ggplot(meta_diversity, aes(x=DOC, y=Chao1)) + geom_point(alpha = 0.3)
mode10 <- lmer(log(Chao1) ~ DOC + (1|Site), meta_diversity)
summary(mode10)
MuMIn::r.squaredGLMM(mode10)
mode11 <- lmer(Shannon ~ DOC + (1|Site), meta_diversity)
summary(mode11)
MuMIn::r.squaredGLMM(mode11)
mode12 <- lmer(log(Simpson) ~ DOC + (1|Site), meta_diversity)
summary(mode12)
MuMIn::r.squaredGLMM(mode12)
mode13 <- lmer(log(Chao1) ~ SUVA254 + (1|Site), meta_diversity)
summary(mode13)
MuMIn::r.squaredGLMM(mode13)
mode14 <- lmer(Shannon ~ SUVA254 + (1|Site), meta_diversity)
summary(mode14)
MuMIn::r.squaredGLMM(mode14)
mode15 <- lmer(log(Simpson) ~ SUVA254 + (1|Site), meta_diversity)
summary(mode15)
MuMIn::r.squaredGLMM(mode15)
mode16 <- lmer(log(Chao1) ~ a320 + (1|Site), meta_diversity)
summary(mode16)
MuMIn::r.squaredGLMM(mode16)
mode17 <- lmer(Shannon ~ a320 + (1|Site), meta_diversity)
summary(mode17)
MuMIn::r.squaredGLMM(mode17)
mode18 <- lmer(log(Simpson) ~ a320 + (1|Site), meta_diversity)
summary(mode18)
MuMIn::r.squaredGLMM(mode18)
mode19 <- lmer(log(Chao1) ~ pH + (1|Site), meta_diversity)
summary(mode19)
MuMIn::r.squaredGLMM(mode19)
mode20 <- lmer(Shannon ~ pH + (1|Site), meta_diversity)
summary(mode20)
MuMIn::r.squaredGLMM(mode20)
mode21 <- lmer(log(Simpson) ~ pH + (1|Site), meta_diversity)
summary(mode21)
MuMIn::r.squaredGLMM(mode21)

# diversity plot
melted <- melt(meta_diversity[,c("Chao1", "Shannon", 'Simpson', "Region")], id.vars = c("Region"))
alpha_region_plot <- ggplot(melted, aes(x = Region, y = value, fill = Region)) +
  geom_violin(trim=T, width=0.5, aes(fill = Region), colour = "#000000") +
  scale_x_discrete(limits = c('Tibetan Plateau', 'Pan-Arctic')) +
  scale_fill_manual(values= c('#d95f02', '#1b9e77')) +
  geom_boxplot(width=0.1, fill="white", colour = "#000000") +
  labs(x = NULL, y = 'Alpha diversity') +
  facet_wrap(~variable, scales = 'free') +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12, 
                                   angle = 45, hjust = 1),
        panel.grid = element_blank(), legend.position = 'none')
alpha_region_plot

alpha_site_plot <- plot_richness(meta_physeq, x="Sitegroup", color = 'Sitegroup', 
                                 measures=c("Chao1", "Shannon", 'Simpson'))+
  geom_boxplot(aes(fill = Sitegroup), alpha=0.2)+
  theme(strip.text = element_text(size = 14),
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12,
                                   angle = 90, hjust = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
alpha_site_plot


