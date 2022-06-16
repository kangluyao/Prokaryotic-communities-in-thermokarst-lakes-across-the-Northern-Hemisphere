##prune out Order below 1% in each sample and prevalence lower than 10/100 at class level
meta_physeq_rel <- microbiome::transform(meta_physeq, "compositional")
meta.com.ord <- microbiome::aggregate_rare(meta_physeq_rel, level = "Order", 
                                           detection = 1/100, prevalence = 10/100)

order_table <- otu_table(meta.com.ord)
order_table[1:5, 1:5]

allmean <- rowMeans(order_table)

regiontype <- as.factor(sample_data(meta.com.ord)$Region)
table(regiontype)
mean_region <- data.frame(Pan_Arctic = rowMeans(otu_table(subset_samples(meta.com.ord, Region == 'Pan-Arctic'))),Tibetan_Plateau = rowMeans(otu_table(subset_samples(meta.com.ord, Region == 'Tibetan Plateau'))))


# mean_region <- sapply(levels(regiontype), function(i) {
#   temp_phylo = subset_samples(meta.com.ord, Region == i)
#   rowMeans(otu_table(temp_phylo))
# })
rel_abun_dat_ord <- data.frame(Order = rownames(mean_region), mean_region, All_mean = allmean)
rel_abun_dat_ord <- dplyr::arrange(rel_abun_dat_ord, desc(All_mean))
Order_levels <- rel_abun_dat_ord$Order[!rel_abun_dat_ord$Order %in% c('Other', "uncultured")]

tax_table <- as.data.frame(t(data.frame(order_table))) %>% dplyr::select(Order_levels)
ncol(tax_table)
meta_env <-  metadata %>% dplyr::select(c('DOC', 'SUVA254', 'a320',
                                          'MAP', 'MAT', 'pH'))


meta_env[1:5, 1:6]

tax_env_table <-  data.frame(tax_table, meta_env)

write.csv(tax_env_table, './meta_analysis/results/tables/tax_env_table.csv')
outlierKD <- function(dt, var) {
  var_name <- eval(substitute(var),eval(dt))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  boxplot(var_name, main="With outliers")
  hist(var_name, main="With outliers", xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main="Without outliers")
  hist(var_name, main="Without outliers", xlab=NA, ylab=NA)
  title("Outlier Check", outer=TRUE)
  na2 <- sum(is.na(var_name))
  cat("Outliers identified:", na2 - na1, "n")
  cat("Propotion (%) of outliers:", round((na2 - na1) / sum(!is.na(var_name))*100, 1), "n")
  cat("Mean of the outliers:", round(mo, 2), "n")
  m2 <- mean(var_name, na.rm = T)
  cat("Mean without removing outliers:", round(m1, 2), "n")
  cat("Mean if we remove outliers:", round(m2, 2), "n")
  response <- readline(prompt="Do you want to remove outliers and to replace with NA? [yes/no]: ")
  if(response == "y" | response == "yes"){
    dt[as.character(substitute(var))] <- invisible(var_name)
    assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
    cat("Outliers successfully removed", "n")
    return(invisible(dt))
  } else{
    cat("Nothing changed", "n")
    return(invisible(var_name))
  }
}

tax_env_table <-  data.frame(tax_table, meta_env)
outlierKD(tax_env_table, Burkholderiales)
y
outlierKD(tax_env_table, SUVA254)
y
outlierKD(tax_env_table, DOC)
y
outlierKD(tax_env_table, a320)
y
outlierKD(tax_env_table, MAP)
y
outlierKD(tax_env_table, MAT)
y
# outlierKD(tax_env_table, SUVA254)
nrow(tax_env_table)
p_doc <- ggplot(tax_env_table, aes(x = DOC, y = Burkholderiales)) +
  geom_point(shape = 19, size = 2, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), size = 1.5, se = T, colour = 'black') +
  # stat_regline_equation(
  #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
  #   formula =y ~ poly(x, 2, raw = TRUE),
  #   label.x.npc = 0.45, label.y.npc = 0.98, size = 5) +
  ggpubr::stat_cor(aes(DOC, Burkholderiales, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                   method = "pearson", label.x.npc = 0.1, label.y.npc = 0.1, size = 5) +
  # xlab('SUVA254')+ylab('Burkholderiales') +
  # scale_y_continuous(expand = c(0,0), limits = c(0.2, 1)) +
  # scale_x_continuous(expand = c(0,0), limits = c(200, 550)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank())+ 
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_text(colour='black', size=14),
        axis.title.y = element_text(colour='black', size=14),
        axis.text.y = element_text(colour='black', size=12),
        axis.text.x = element_text(colour = "black", size = 12))

p_suv254 <- ggplot(tax_env_table, aes(x = SUVA254, y = Burkholderiales)) +
  geom_point(shape = 19, size = 2, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), size = 1.5, se = T, colour = 'black') +
  # stat_regline_equation(
  #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
  #   formula =y ~ poly(x, 2, raw = TRUE),
  #   label.x.npc = 0.45, label.y.npc = 0.98, size = 5) +
  ggpubr::stat_cor(aes(SUVA254, Burkholderiales, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                   method = "pearson", label.x.npc = 0.1, label.y.npc = 0.1, size = 5) +
  # xlab('SUVA254')+ylab('Burkholderiales') +
  # scale_y_continuous(expand = c(0,0), limits = c(0.2, 1)) +
  # scale_x_continuous(expand = c(0,0), limits = c(200, 550)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank())+ 
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_text(colour='black', size=14),
        axis.title.y = element_text(colour='black', size=14),
        axis.text.y = element_text(colour='black', size=12),
        axis.text.x = element_text(colour = "black", size = 12))

p_a320 <- ggplot(tax_env_table, aes(x = a320, y = Burkholderiales)) +
  geom_point(shape = 19, size = 2, colour ='tomato3', alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 1, raw = TRUE), size = 1.5, se = T, colour = 'black') +
  # stat_regline_equation(
  #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
  #   formula =y ~ poly(x, 2, raw = TRUE),
  #   label.x.npc = 0.45, label.y.npc = 0.98, size = 5) +
  ggpubr::stat_cor(aes(a320, Burkholderiales, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                   method = "pearson", label.x.npc = 0.1, label.y.npc = 0.1, size = 5) +
  # xlab('SUVA254')+ylab('Burkholderiales') +
  # scale_y_continuous(expand = c(0,0), limits = c(0.2, 1)) +
  # scale_x_continuous(expand = c(0,0), limits = c(200, 550)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank())+ 
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_text(colour='black', size=14),
        axis.title.y = element_text(colour='black', size=14),
        axis.text.y = element_text(colour='black', size=12),
        axis.text.x = element_text(colour = "black", size = 12))


## arrange the plots
plot_DOM_Burkho <- plot_grid(p_doc, p_suv254, p_a320,
                              labels = c('A', 'B', 'C'), 
                              ncol = 3, nrow = 1, 
                              label_x = .01, label_y = 1.01, hjust = 0, 
                              label_size=14, align = "v")
plot_DOM_Burkho


#correlation between diversity and environmental variables
x<-log(tax_table)
y<-scale(meta_env)
groups<-c(rep('All',306))

#You can use kendall, spearman, or pearson below:
method<-"spearman"


#Now calculate the correlation between individual Taxa and the environmental data
df<-NULL
for(i in colnames(x)){
  for(j in colnames(y)){
    for(k in unique(groups)){
      a<-x[groups==k,i,drop=F]
      b<-y[groups==k,j,drop=F]
      tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],
                     use="everything",method=method),
             cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value,k)
      if(is.null(df)){
        df<-tmp  
      }
      else{
        df<-rbind(df,tmp)
      }    
    }
  }
}

df<-data.frame(row.names=NULL,df)
colnames(df)<-c("Index","Env","Correlation","Pvalue","Type")
df$Pvalue<-as.numeric(as.character(df$Pvalue))
df$AdjPvalue<-rep(0,dim(df)[1])
df$Correlation<-as.numeric(as.character(df$Correlation))
#You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
# 1 -> donot adjust
# 2 -> adjust Env + Type (column on the correlation plot)
# 3 -> adjust Taxa + Type (row on the correlation plot for each type)
# 4 -> adjust Taxa (row on the correlation plot)
# 5 -> adjust Env (panel on the correlation plot)
adjustment_label<-c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")
adjustment<-5

if(adjustment==1){
  df$AdjPvalue<-df$Pvalue
} else if (adjustment==2){
  for(i in unique(df$Env)){
    for(j in unique(df$Type)){
      sel<-df$Env==i & df$Type==j
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==3){
  for(i in unique(df$Taxa)){
    for(j in unique(df$Type)){
      sel<-df$Taxa==i & df$Type==j
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==4){
  for(i in unique(df$Taxa)){
    sel<-df$Taxa==i
    df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
  }
} else if (adjustment==5){
  for(i in unique(df$Env)){
    sel<-df$Env==i
    df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
  }
}

#Now we generate the labels for signifant values
df$Significance<-cut(df$Pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

#We ignore NAs
df<-df[complete.cases(df),]

#We want to reorganize the Env data based on they appear
df$Env<-factor(df$Env,levels=c("MAT","MAP",'DOC', 'SUVA254', 'a320', 'pH'))
df$Index<-factor(df$Index,levels = rev(c(Order_levels)))
#df$Taxa<-factor(df$Taxa,levels=rev(taxa_list))
#df$Correlation[which(abs(df$AdjPvalue) >= 0.05)]<-0
#We use the function to change the labels for facet_grid in ggplot2
#


#plot
p_env_taxa_fun <- ggplot(aes(x=Env, y=Index, fill=Correlation), data=df)+
  geom_tile() +
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C")+
  geom_text(aes(label=Significance), color="black", size=4)+labs(y=NULL, x=NULL, fill=method)+
  theme(axis.title =element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))


p_env_taxa_fun <-ggplot(aes(x=Env, y=Index, col=Correlation), data=df) +
  geom_point(aes(size = abs(Correlation)))+ #,shape=1
  scale_color_gradient2(low='#1b9e77', mid='white', high='#d95f02') +
  theme_bw(base_size = 12)+
  geom_text(aes(label=Significance), color="black", size = 5)+
  scale_size(range = c(1, 8),breaks=c(0.1,0.2,0.4,0.6,0.8), limits=c(0,0.7))+ 
  #theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = 'Environmental factors', y = '', col = "Spearman's r")+ #legend title
  theme(panel.background = element_rect(fill='white', colour='black'),
        #panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.25,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12, angle = 45, hjust = 1),
        legend.position = 'right',
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))

print(p_env_taxa_fun)


#
env_tax_lcbd_plot <- ggdraw() +
  draw_plot(p_env_taxa_fun, x = 0, y = 0, width = 0.55, height = 1) +
  draw_plot(meta_lcbd_map, x = 0.55, y = 0.5, width = 0.4, height = 0.5) +
  draw_plot(meta_c_lcbd_map, x = 0.55, y = 0, width = 0.4, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), size = 14,
                  x = c(0, 0.55, 0.55), y = c(1, 1, 0.5))

env_tax_lcbd_plot









#edgeR analysis
# Glom OTUs to genus level for further statistical analysis & reasonable power
ps1 <- tax_glom(meta.com.ord, "Order", NArm = TRUE)
dat = data.frame(sample_data(ps1))
i <- rownames(dat[complete.cases(dat), ]) # remove the samples with NA
ps2 <- prune_samples(i, ps1)
#Use code snippet straight provided by authors of PhyloSeq to export phyloseq object to an EdgeR object
phyloseq_to_edgeR = function(physeq, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}

# Make normalized phyloseq object (ps6) into an edgeR object. It needs a grouping factor. We use location.
dge = phyloseq_to_edgeR(ps2)
#The crunching to follow is much easier if metadata is pulled out this way into object "a"
dat = data.frame(sample_data(ps2))

#So many things don't understand what kind of variables they need to be. Make sure they understand
chitinolysis <- as.numeric(a$chitinolysis)
cellulolysis <- as.numeric(a$cellulolysis)
fermentation <- as.numeric(a$fermentation)
methanogenesis <- as.numeric(a$methanogenesis)
methanotrophy <- as.numeric(a$methanotrophy)
methylotrophy <- as.numeric(a$methylotrophy)
DOC <- as.numeric(a$DOC)
SUVA254 <- as.numeric(a$SUVA254) 
a320 <- as.numeric(a$a320)
pH <- as.numeric(a$pH)
MAT <- as.numeric(a$MAT)
MAP <- as.numeric(a$MAP)

# ##Test and forward selection of environmental variables
# Env_select <- function (Dist_Matrix,Env,Number_Permutations=999) {
#   mod1<-capscale(Dist_Matrix~.,Env,add=T,na.action = na.omit)
#   mod0<-capscale(Dist_Matrix~1,Env,add=T,na.action = na.omit)
#   mod<-ordiR2step(mod0,scope=formula(mod1),perm.max=999,trace = F)
#   Env.P<-anova.cca(mod,permutations=Number_Permutations)[[4]][1]
#   Env_se<-mod$CCA$biplot
#   Env_se<-Env[,rownames(Env_se)]
#   return(Env_se)
# }
# carbon_asvs <- otu_table(carbon_phy_rel)
# carbon.dist <- vegdist(sqrt(t(carbon_asvs)))
# env.table <- sample_data(carbon_phy_rel)
# env_vars <- data.frame(scale(env.table[ ,-c(1:6)]))
# env_se_carbon <- Env_select (carbon.dist, env_vars)
# colnames(env_se_carbon)
# Design for my linear model
design <- model.matrix(~ MAP + MAT + pH + SUVA254 + DOC + a320, model.frame(~ ., dat, na.action=na.pass))

# EdgeR needs to calculate dispersion again after you've fed it the design.
#This step projectile vomits if chloride or conductivity aren't sqrt transformed.
x = calcNormFactors(dge, method="RLE")
x = estimateGLMCommonDisp(dge, design)
x = estimateGLMTrendedDisp(dge, design)
x = estimateGLMTagwiseDisp(dge, design)

fit <-glmFit(x, design)

# grab the coefficients I care about
lrt <- glmLRT(fit, coef=2:7)

# lrt to z-scores
table<-lrt$table
table<-apply(table, 2, function(x) scale(x, center = TRUE, scale = TRUE))
rownames(table) = rownames(lrt$table)
table<-table[,1:6]
table<-as.data.frame(table)

q<-lrt$genes
table1<-cbind(table, Order = q$Order)

# Sometimes want to do this to check what has at least one z score beyond threshold
table2 <- table1[apply(table1[,1:6], 1, function(x) any(abs(x)>1.96)), ]

#write.csv(table2, "./tables/EdgeR-Zscores-Genera.csv")
melted<-melt(table2)
melt2<-subset(melted, abs(melted$value) > 1.96 )

filter<-unique(melt2$OTU)
mini<-prune_taxa((rownames(otu_table(water_physeq)) %in% filter), water_physeq)
f<-phy_tree(mini)
tax = as.data.frame(tax_table(mini))

tax$Genus = gsub("D_5__", "", tax$Genus)
tax$Phylum = gsub("D_1__", "", tax$Phylum)

# change tip labels to genera / need to gsub
library(ggtree)
mytree = phy_tree(f)$tip.label
foo = data.frame(label=mytree, label2=paste(tax$Genus, tax$ASV, sep='_'), label3=tax$Phylum, stringsAsFactors = F)
#foo$label2[foo$label2 == NA] <- 'uncultured'
# Make a less putrid color palette

library(RColorBrewer)

colourCount <- length(unique(tax$Phylum))
Mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(colourCount)
names(Mycolors) <- levels(as.factor(tax$Phylum))
colScale <- scale_colour_manual(name = "grp",values = Mycolors)

p<-ggtree(f) %<+% foo + 
  geom_tiplab(size=2, align=TRUE, linesize=0.5, aes(label=label2), hjust =0) +
  geom_tippoint(size=2, aes(label=label2, group=label3, color=label3) ) + 
  theme_tree2() + 
  scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(colourCount))
p + xlim (0, 3) 

test<-as.matrix(table2[1:20])
tested<-apply(test, 2, function(x){cut(x, br=c(-14, -4, -2, 2, 4, 14))})
rownames(tested) = rownames(table2)
colnames(tested) = gsub("logFC.", "", colnames(tested))
colnames(tested) = gsub("logCPM", "Rotorua", colnames(tested))

heatmapcols <-colorRampPalette(brewer.pal(6, "Spectral"))(6)
names(heatmapcols) <- levels(as.factor(tested[1:10000]))
heatmapcols <- c(heatmapcols[5],heatmapcols[4],heatmapcols[2],heatmapcols[3],heatmapcols[1],heatmapcols[6])
# Make sure to check that the order of these is correct if this is rerunhe
heatmapcols[1] = "#881111"
heatmapcols[2] = "#FFCCCF"
heatmapcols[3] = "#F2F2F2"
heatmapcols[4] = "#AABBDD"
heatmapcols[5] = "#112288"
heatmapcols[6] = "#112288"

mycols <- c("#881111","#D9444D","#FFCCCF","#F2F2F2","#AABBDD","#112288")
tested <- tested[,c('chitinolysis', 'cellulolysis', 'fermentation', 'methanogenesis', 
                    'methanotrophy', 'methylotrophy',  'pH', 'MAP', 'MAT',
                    'Comp4', 'Comp3', 'S275_295', 'SUVA254', 'DOC', 'NH4_N',
                    'Na', 'K', 'Ca', 'Conductivity', 'Depth')] 
p1<-gheatmap(p, tested, offset = 0.25, width = 3, font.size=2.5, 
             colnames_angle=-90, colnames_position = "top", hjust=1) + 
  scale_fill_manual(values=heatmapcols)
p1 <-p1 + theme(legend.position="right")
p1











































































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
## arrange data
rel_abun_dat_ord_plot <- data.frame(Order = rownames(rel_abun_dat_ord), prop = rel_abun_dat_ord$All_mean * 100)

rel_abun_dat_ord_plot <- rel_abun_dat_ord_plot[match(c(rel_abun_dat_ord_plot$Order[-2], 'Other'),
                                                     rel_abun_dat_ord_plot$Order),]
rel_abun_dat_ord_plot <- rel_abun_dat_ord_plot %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
rel_abun_dat_ord_plot$Order <- factor(rel_abun_dat_ord_plot$Order, ordered = T,
                                      levels = rel_abun_dat_ord_plot$Order)

## pie plot
pie_for_taxa <- ggplot(rel_abun_dat_ord_plot, aes(x="", y=prop, fill=reorder(Order, -lab.ypos))) +
  geom_bar(width = 1, stat = "identity", color = "white") + 
  coord_polar("y", start=0) +
  geom_text(aes(x = 1.35, y = lab.ypos, 
                label = paste0(round(prop, 1), '%',  sep= '')),
            color = "white", size=3) +
  scale_fill_manual('Order', values = mycols) +
  guides(fill = guide_legend(reverse=T)) +
  theme_void()
pie_for_taxa
