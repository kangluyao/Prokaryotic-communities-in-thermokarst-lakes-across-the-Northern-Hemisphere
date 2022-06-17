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