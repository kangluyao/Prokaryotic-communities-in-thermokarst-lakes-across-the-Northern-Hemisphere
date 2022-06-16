## test the relationship between LCBD and MAP across TP and PA
library(adespatial)
###total community
beta_meta_div <- beta.div(t(as.matrix(otu_table(meta_physeq))), 
                          method = "hellinger", sqrt.D = FALSE, samp = TRUE, 
                          nperm = 999, adj = TRUE, save.D = FALSE, clock = FALSE)

env_div <- data.frame(LCBD = beta_meta_div$LCBD, sample_data(meta_physeq))

### calculate the average LCBD for each site
library(dplyr)
env_div_agg_meta <-  env_div %>% 
  dplyr::select(c(1, 5, 6, 7, 8, 19)) %>%
  group_by(Sitegroup, Site, Sitegroup1) %>%
  dplyr::summarise(across(, mean, na.rm = TRUE))

ggplot(env_div, aes(x=pH, y=LCBD)) + geom_point(alpha = 0.3)
lm1 <- lm(LCBD ~  pH + I(pH^2), data = env_div)
summary(lm1)
ggplot(env_div, aes(x=DOC, y=LCBD)) + geom_point(alpha = 0.3)
lm2 <- lm(LCBD ~  DOC, data = env_div)
summary(lm2)
ggplot(env_div, aes(x=a320, y=LCBD)) + geom_point(alpha = 0.3)
lm3 <- lm(LCBD ~  a320 + I(a320^2), data = env_div)
summary(lm3)
ggplot(env_div, aes(x=SUVA254, y=LCBD)) + geom_point(alpha = 0.3)
lm4 <- lm(LCBD ~  SUVA254, data = env_div)
summary(lm4)
ggplot(env_div, aes(x=log(MAT+15), y=LCBD)) + geom_point(alpha = 0.3)
lm5 <- lm(LCBD ~ log(MAT+15), data = env_div)
summary(lm5)
ggplot(env_div, aes(x=MAP, y=LCBD)) + geom_point(alpha = 0.3)
lm6 <- lm(LCBD ~  MAP + I(MAP^2), data = env_div)
summary(lm6)

## random forest analysis
### Regression:
env_div_rf <- env_div %>%
  select(c('LCBD',  'MAT', 'MAP', 'DOC', 'SUVA254', 'a320', 'pH'))
library(randomForest)
library(rfPermute)
set.seed(131)
lcbd.rf <- rfPermute(LCBD ~ ., data=env_div_rf, ntree=1000, num.rep = 100,
                     importance=TRUE, na.action=na.omit)

impor.dat <-data.frame(variables = rownames(importance(lcbd.rf)), IncMSE=importance(lcbd.rf)[,1])

library(dplyr)
impor_meta_plot <- impor.dat %>% 
  arrange(IncMSE)  %>%
  mutate(variables=factor(variables, levels=variables)) %>%
  ggplot(aes(x = variables, y = IncMSE)) +
  geom_bar(stat="identity", fill="#d95f02", colour=NA, width = 0.65) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 45)) +
  ylab('Increases in MSE %') + xlab(NULL) +
  theme_classic() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size =12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.6, 0.75)) +
  coord_flip()

## Show "importance" of variables: higher value mean more important:
mean(lcbd.rf$rf$rsq)
sd(lcbd.rf$rf$rsq)

meta_lcbd_map <- ggplot(env_div, aes(MAP, LCBD, col = Sitegroup)) +
  geom_smooth(method="lm", formula = y ~ poly(x, 2),
              size=1, se=T, colour='black') +
  geom_point(size=1, alpha=0.8) +
  xlab('MAP (mm)') + ylab('LCBD') +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.6, 0.75))

rf_plot <- cowplot::plot_grid(impor_meta_plot, meta_lcbd_map,
                              labels = c('(a)', '(b)'), ncol = 2, 
                              label_x = .01, label_y = 1, 
                              hjust = 0, label_size = 14, align = "v")

#ggsave(rf_plot, file = "./meta_analysis/results/figs/rf_plot.pdf",
#       width = 13, height = 6, units = 'in', device='pdf', dpi=300)

#alpine dataset analysis
# beta diversity analysis
## local contribution to beta diversity (LCBD) analysis for total community
library(adespatial)
beta_tax_div <- beta.div(t(as.matrix(otu_table(tp_physeq))), 
                         method = "hellinger", sqrt.D = FALSE, samp = TRUE, 
                         nperm = 999, adj = TRUE, save.D = FALSE, clock = FALSE)
env_div <- data.frame(LCBD = beta_tax_div$LCBD, sample_data(tp_physeq))
### calculate the average LCBD for each site
library(dplyr)
env_div_agg <-  env_div %>% 
  dplyr::select(-c(2:6, 9, 12)) %>%
  group_by(Site, Sitegroup1) %>%
  dplyr::summarise(across(, mean, na.rm = TRUE))
#write.csv(env_div_agg, file = './tibet_dada2_asv/results/tables/env_div_agg.csv')

### determine the relationships between LCBD and envs using linear regression modes
vars <- c("MAT", "MAP", "DOC", "S275_295", "SUVA254", "a300", "BIX", "HIX",
          "TN", "NH4_N", "DO", "pH", "Conductivity", "Salinity",
          "K", "Ca", "Na", "Mg")
mode <- lapply(vars, function(x) {
  lm(substitute(LCBD ~ i, list(i = as.name(x))), data = env_div_agg)})
sum.mode <- lapply(mode, broom::glance)
### normality test using Shapiro-Wilk test 
res <- lapply(mode, residuals)
norm_test <- lapply(res, shapiro.test)
norm_results <- data.frame(
  variables = vars, 
  w = sapply(norm_test, "[[", "statistic"), 
  pvalue = sapply(norm_test, "[[", "p.value")
)
norm_results

### extract the standardized regression coefficients
sd.coeff <- lapply(mode, QuantPsyc::lm.beta)
### arrange the table for plot
LCBD <- c(rep('LCBD', length(vars)))
sd.coeff <- sapply(sd.coeff, function(x){as.numeric(x[1])})
r.squared <- sapply(sum.mode, "[[", "r.squared")
adj.r.squared <- sapply(sum.mode, "[[", "adj.r.squared")
pvalue <- sapply(sum.mode, "[[", "p.value")
sig <- cut(pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
results <- data.frame(vars, LCBD, sd.coeff, r.squared, adj.r.squared, pvalue, sig)
results

#model selection
library(MASS)
library(glmulti)
A1 <- glmulti(LCBD ~ MAP + MAT + S275_295 + SUVA254 + a300 + BIX + HIX +
                TN + Conductivity + Salinity + Mg + K + Na, data=env_div_agg,
              level=1, fitfunction=lm, crit="aicc", confsetsize= 2^13, plotty = F, trace = 0)
top <- weightable(A1)
###  models with values more than 2 units away are considered substantially 
### less plausible than those with AICc values closer to that of the best model. 
### refrence:Anderson, D. R. (2007). Model based inference in the life sciences: A primer on evidence. New York: Springer. 
top_1 <- top[top$aicc <= min(top$aicc) + 2,] # 
top_1

modes_inf <- NULL
for(i in 1:nrow(top_1)){
  rse_sum <- summary(A1@objects[[i]])
  adj.r.squared <- rse_sum$adj.r.squared # obtain the adjust r squared
  multicollinearity <- any(car::vif(A1@objects[[i]]) > 2) # check the multicollinearity
  tmp <- data.frame(adj.r.squared, multicollinearity)
  if(is.null(modes_inf)){
    modes_inf<-tmp
  } else {
    modes_inf <- rbind(modes_inf,tmp)
  } 
}
modes_inf <- cbind(top_1, modes_inf)
modes_inf

vpa.mod <- varpart(env_div_agg$LCBD, ~ env_div_agg$HIX,
                   ~ env_div_agg$MAP)
plot(vpa.mod)

### heatmap using standardized regression coefficients to explore the relationship between the LCBD and environment factors
results_plot_data <- data.frame(group = c(rep('Total community', nrow(results)), rep('Carbon cycling', nrow(results_carbon))),
                                rbind(results, results_carbon))

results_plot_data$vars <- factor(results_plot_data$vars,levels = rev(vars))
results_plot_data$group <- factor(results_plot_data$group,
                                  levels=c('Total community', 'Carbon cycling'))
p_env_div <- ggplot(aes(x=LCBD, y=vars, fill=sd.coeff), data=results_plot_data) +
  geom_tile() +
  scale_fill_gradient2(low='#1b9e77', mid='white', high='#d95f02') +
  geom_text(aes(label=sig), color="black", size=6) +
  labs(y=NULL, x=NULL, fill='Standardized regression coefficients') +
  facet_wrap( .~ group, ncol = 2) +
  theme_bw()+
  theme(legend.position="bottom", 
        panel.border = element_blank(),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        legend.title=element_text(size = 12),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "white"))
###  plot the linear regression relationships between the LCBD and best explained variables
#### total community
p_linear <- env_div_agg %>%
  dplyr::select(LCBD, MAP, HIX) %>%
  tidyr::gather(varibales, value, MAP:HIX, factor_key=TRUE) %>%
  ggplot(aes(value, LCBD)) +
  geom_point(size=3.5, alpha=0.8, aes(colour = as.factor(varibales))) +
  geom_smooth(method="lm", size=1, se=T, colour='black') +
  scale_color_manual(values = c('#1b9e77', '#d95f02')) +
  scale_y_continuous(limits = c(0, 0.01)) +
  facet_wrap( .~ varibales, scales="free_x", ncol = 2) +
  ylab('LCBD')+xlab('Values') +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        strip.text = element_text(size = 14),
        legend.position='none')










## partial mantel test function for Tibatan Plateau and Pan Arctic
partial.mantel.fun <- function(phylo) {
  env.table <- data.frame(sample_data(phylo))
  #env.table <- env.table[complete.cases(env.table), ]
  otu_table <- as.matrix(t(otu_table(phylo)))
  otu_table_hel <- decostand(otu_table, 'hellinger')
  otu_table_hel_dist <- vegdist(otu_table_hel, 'bray',upper=F)
  df <- NULL
  vars <- c("DOC", "SUVA254", "a320", "MAP", "MAT", "pH")
  for (x in vars) {
    x.dist <- vegdist(scale(env.table[,x]), 'euclidean', na.rm = T)
    z.dist <- vegdist(scale(env.table[ , setdiff(vars, x)]), 'euclidean', na.rm = T)
    mode <- mantel.partial(x.dist, otu_table_hel_dist, z.dist, 
                           method = "pearson", permutations = 999, na.rm = T)
    r <- mode$statistic
    p <- mode$signif
    tmp <- data.frame(env = x, r = r, p.value = p)
    if(is.null(df))
      df <- tmp
    else
      df <- rbind(df ,tmp)
  }
  return(df)
}

all.total.commun.par.mant <- partial.mantel.fun(meta_physeq)
all.carbon.commun.par.mant <- partial.mantel.fun(meta_carbon_phy)
## PLOT
## devtools::install_github('hannet91/ggcor')
library(ggcor)
set.seed(123456)

all.par.man.tibble <- tibble(spec = c(rep('all.total.commun.par.mant', nrow(all.total.commun.par.mant)), 
                                      rep('all.carbon.commun.par.mant', nrow(all.carbon.commun.par.mant))), 
                             rbind(all.total.commun.par.mant, all.carbon.commun.par.mant))
#par.man.tibble <- tibble(spec = c(rep('Total community composition', nrow(total.commun.par.mant))), 
#                         total.commun.par.mant)
vars <- c("MAT", "MAP", "DOC", "SUVA254", "a320", "pH")
env.table <- sample_data(meta_physeq)[ , vars]

mantel02 <- all.par.man.tibble %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.3, 0.4, Inf), 
                 labels = c("<0.30", "0.30-0.40", ">0.40"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">0.05"),
                       right = T))
quickcor(env.table[complete.cases(env.table),], type = "upper") + geom_square() + 
  add_link(mantel02, mapping = aes(colour = p.value, size = r),
           diag.label = TRUE) +
  scale_color_manual(values = c('#d95f02', 'grey', '#1b9e77', '#3C5488FF')) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  geom_diag_label() + remove_axis("x")

