# melt to long format (for ggplot) 
otu_num_in_phylum <- data.frame(table(tax_table(meta_physeq)[,"Phylum"]))
otu_num_in_phylum <- otu_num_in_phylum %>% arrange(desc(Freq))

# pie plot for OTU number within each phylum
## data preparation
otu_num_in_phylum <- rbind(otu_num_in_phylum[1:11, ], data.frame(Var1 = c('Others'), Freq = sum(otu_num_in_phylum[-c(1:11), 2])))

otu_num_in_phylum <- data.frame(Phylum = otu_num_in_phylum$Var1, otu_num = otu_num_in_phylum$Freq ,
                                prop = otu_num_in_phylum$Freq/sum(otu_num_in_phylum$Freq)*100)
otu_count.data <- otu_num_in_phylum %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
otu_count.data$Phylum <- factor(otu_count.data$Phylum, ordered = T, levels = otu_num_in_phylum$Phylum)


library(RColorBrewer)
## Define the colors you want
mycols <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
            "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
            "#8A7C64", "#599861")
## Use scale_fill_manual
pie_for_otu_num_phylum <- ggplot(otu_count.data, aes(x="", y=prop, fill=reorder(Phylum, -lab.ypos))) +
  geom_bar(width = 1, stat = "identity", color = "white") + 
  coord_polar("y", start=0) +
  geom_text(aes(x = 1.35, y = lab.ypos, 
                label = paste0(otu_num, ' (', round(prop, 1), '%', ')', sep= '')),
                color = "white", size=3) +
  scale_fill_manual('Phylum', values = mycols) +
  guides(fill = guide_legend(reverse=T)) +
  theme_void()
pie_for_otu_num_phylum


##prune out Order below 1% in each sample and prevalence lower than 10/100 at class level
meta.com.ord <- microbiome::aggregate_rare(meta_physeq_rel, level = "Order", 
                                           detection = 1/100, prevalence = 10/100)
allmean <- rowMeans(otu_table(meta.com.ord))

regiontype <- as.factor(sample_data(meta.com.ord)$Region)
table(regiontype)
mean_region <- sapply(levels(regiontype),function(i){
  rowMeans(otu_table(meta.com.ord)[,region == i])
})
rel_abun_dat_ord <- data.frame(Phylum = rownames(mean_region), mean_region, All_mean = allmean)
rel_abun_dat_ord <- dplyr::arrange(rel_abun_dat_ord, desc(All_mean))
rel_abun_dat_ord

# pie plot for total community composition
## pie plot
rel_abun_dat_ord_plot <- data.frame(Order = rownames(rel_abun_dat_ord), prop = rel_abun_dat_ord$All_mean * 100)
rel_abun_dat_ord_plot$Order <- factor(rel_abun_dat_ord_plot$Order, ordered = T,
                                      levels = c(rel_abun_dat_ord_plot$Order[-2], 'Other'))
rel_abun_dat_ord_plot$ypos = cumsum(rel_abun_dat_ord_plot$prop)- 0.5 * rel_abun_dat_ord_plot$prop 

library(RColorBrewer)
## Define the number of colors you want
nb.cols <- 9
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
## Create a ggplot with 11 colors 
## Use scale_fill_manual
pie_for_taxa <- ggplot(rel_abun_dat_ord_plot, aes(x="", y=prop, fill=Order)) +
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) +
  scale_fill_manual(values = mycolors) +
  theme_void() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14)
  ) +
  geom_text(aes(y = ypos, label = scales::percent(prop/100),
                color = "white", size=6))
pie_for_taxa


# plot the community composition for TP and PA at order level
plot.composition.relAbun <- microbiome::plot_composition(meta.com.ord, 
                                                         average_by = "Region", 
                                                         otu.sort = "abundance") +
  scale_fill_manual(values = mycols) + 
  scale_x_discrete(limits = c('Tibetan Plateau', 'Pan-Arctic')) +
  scale_y_continuous(labels = function(x) paste0(x*100), expand = c(0, 0)) + 
  labs(x = NULL, y = 'Mean relative abundance (%)') +
  theme_bw() +
  guides(fill=guide_legend(ncol=1)) +
  theme(legend.position = 'right',
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12,
                                   angle = 45, hjust = 1),
        legend.title = element_text(size = 14, colour = 'black'),
        legend.text = element_text(size = 12, colour = 'black'))
plot.composition.relAbun 



## determine the composition within gammaproteobacteria
gammaproteobacteria_phylo <- subset_taxa(meta_physeq, Class == 'Gammaproteobacteria', level = "Order")
gammaproteobacteria_phylo <- tax_glom(gammaproteobacteria_phylo, taxrank="Order")
gammaproteobacteria_phylo_rel <- microbiome::transform(gammaproteobacteria_phylo, "compositional")
gammaproteobacteria_phylo_rel <- aggregate_rare(gammaproteobacteria_phylo_rel, level = "Order", detection = 1/100, prevalence = 10/100)

plot.gammaproteobacteria.relAbun <- microbiome::plot_composition(gammaproteobacteria_phylo_rel, 
                                                                 average_by = "Region", 
                                                                 otu.sort = "abundance") +
  scale_fill_brewer("Order", palette = "Paired") + 
  scale_x_discrete(limits = c('Tibetan Plateau', 'Pan-Arctic')) +
  scale_y_continuous(label = scales::percent, expand = c(0, 0)) + 
  labs(x = NULL, y = 'Relative abundance') +
  theme_bw() +
  theme(legend.position =  'bottom',
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12,
                                   angle = 45, hjust = 1),
        legend.title = element_text(size = 12, colour = 'black'),
        legend.text = element_text(size = 10, colour = 'black'))

## determine the composition within actinobacteria
actinobacteria_phylo <- subset_taxa(meta_physeq, Class == 'Actinobacteria', level = "Order")
actinobacteria_phylo <- tax_glom(actinobacteria_phylo, taxrank="Order")
actinobacteria_phylo_rel <- microbiome::transform(actinobacteria_phylo, "compositional")
actinobacteria_phylo_rel <- aggregate_rare(actinobacteria_phylo_rel, level = "Order", detection = 1/100, prevalence = 10/100)

plot.actinobacteria.relAbun <- microbiome::plot_composition(actinobacteria_phylo_rel, 
                                                            average_by = "Region", 
                                                            otu.sort = "abundance") +
  scale_fill_brewer("Order", palette = "Set1") + 
  scale_x_discrete(limits = c('Tibetan Plateau', 'Pan-Arctic')) +
  scale_y_continuous(label = scales::percent, expand = c(0, 0)) + 
  labs(x = NULL, y = 'Relative abundance') +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12,
                                   angle = 45, hjust = 1),
        legend.title = element_text(size = 12, colour = 'black'),
        legend.text = element_text(size = 10, colour = 'black'))

## determine the composition within bacteroidia
bacteroidia_phylo <- subset_taxa(meta_physeq, Class == 'Bacteroidia', level = "Order")
bacteroidia_phylo <- tax_glom(bacteroidia_phylo, taxrank="Order")
bacteroidia_phylo_rel <- microbiome::transform(bacteroidia_phylo, "compositional")
bacteroidia_phylo_rel <- aggregate_rare(bacteroidia_phylo_rel, level = "Order", detection = 1/100, prevalence = 10/100)

plot.bacteroidia.relAbun <- microbiome::plot_composition(bacteroidia_phylo_rel, 
                                                         average_by = "Region", 
                                                         otu.sort = "abundance") +
  scale_fill_brewer("Order", palette = "Accent") + 
  scale_x_discrete(limits = c('Tibetan Plateau', 'Pan-Arctic')) +
  scale_y_continuous(label = scales::percent, expand = c(0, 0)) + 
  labs(x = NULL, y = 'Relative abundance') +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12,
                                   angle = 45, hjust = 1),
        legend.title = element_text(size = 12, colour = 'black'),
        legend.text = element_text(size = 10, colour = 'black'))
## arrange the plots
domin_class_plot <- plot_grid(plot.gammaproteobacteria.relAbun, plot.actinobacteria.relAbun, plot.bacteroidia.relAbun,
                              labels = c('(a)', '(b)', '(c)'), 
                              ncol = 3, nrow = 1, 
                              label_x = .01, label_y = 1.01, hjust = 0, 
                              label_size=14, align = "v")
domin_class_plot
