library(dplyr)
library(ggplot2)

results=read.csv('./prelim-results.csv') %>%
  filter(!is.na(po_auc))


ggplot(results, aes(po_auc, pa_auc)) +
  geom_point(aes(size=sp_area)) +
  scale_size(range=c(0,20)) +
  geom_abline(intercept=0, slope=1) +
  theme_set(theme_grey(base_size = 28)) 
