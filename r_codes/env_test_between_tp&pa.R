# test the difference of the environmental factors between TP and PA
##standard error function
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
meta_env <- metadata %>% dplyr::select(c('Region', 'DOC', 'SUVA254', 'a320',
                                         'MAP', 'MAT', 'pH')) %>% group_by(Region) %>%
  dplyr::summarise_all(list(mean = mean, sd = sd, se = stderr), na.rm = TRUE)
write.csv(meta_env, file = 'E:/thermokast_lakes/water_microbes/meta_analysis/results/tables/meta_env.csv')
# test the difference
library(lme4)
library(lmerTest)
library(multcomp)
par(mfrow = c(2, 3))
hist(metadata$DOC)
hist(metadata$SUVA254)
hist(metadata$a320)
hist(metadata$MAP)
hist(metadata$MAT)
hist(metadata$pH)
par(mfrow = c(1, 1))

par(mfrow = c(2, 3))
hist(log(metadata$DOC))
hist(log(metadata$SUVA254))
hist(log(metadata$a320))
hist(log(metadata$MAP))
hist(log(metadata$MAT+15))
hist(log(metadata$pH))
par(mfrow = c(1, 1))

env_mode1 <- lmer(log(DOC) ~ Region + (1|Site), metadata)
summary(env_mode1)
env_mode2 <- lmer(log(SUVA254) ~ Region + (1|Site), metadata)
summary(env_mode2)
env_mode3 <- lmer(log(a320) ~ Region + (1|Site), metadata)
summary(env_mode3)
env_mode4 <- lmer(log(MAP) ~ Region + (1|Site), metadata)
summary(env_mode4)
env_mode5 <- lmer(MAT ~ Region + (1|Site), metadata)
summary(env_mode5)
env_mode6 <- lmer(log(pH) ~ Region + (1|Site), metadata)
summary(env_mode6)

par(mfrow = c(2, 3))
qqnorm(resid(env_mode1))
qqnorm(resid(env_mode2))
qqnorm(resid(env_mode3))
qqnorm(resid(env_mode4))
qqnorm(resid(env_mode5))
qqnorm(resid(env_mode6))
par(mfrow = c(1, 1))