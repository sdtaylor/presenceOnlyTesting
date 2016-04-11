library(dplyr)
library(ggplot2)
library(RColorBrewer)

results=read.csv('./results/results.csv') %>%
  filter(!is.na(po_auc), po_kappa>0) %>%
  mutate(point_density=present_sites/sp_area) %>%
  mutate(auc_diff = pa_auc - po_auc, sens_diff = pa_sensitivity - po_sensitivity) %>%
  mutate(sp_area=sp_area/1000)


basegraph=ggplot(results) +
          theme_bw(base_size=28)+
          ylim(0,1.0) + xlim(0,1.0) +
          geom_abline(intercept=0, slope=1) +
          ylab('Presence/Absence Data') +
          xlab('Presence Only Data') +
          scale_size_continuous(name='Species Range\n (1000 sq. km)') 

#Kappa
basegraph + 
  geom_point(aes(x=po_kappa, y=pa_kappa, size=sp_area)) +
  ggtitle('Kappa Statistic')

#Sensitivity
basegraph + 
  geom_point(aes(x=po_sensitivity, y=pa_sensitivity, size=sp_area)) +
  ggtitle('Sensitivity')

#Specificity
basegraph + 
  geom_point(aes(x=po_specificity, y=pa_specificity, size=sp_area)) +
  ggtitle('Specificity')

#Specificity
basegraph + 
  geom_point(aes(x=po_auc, y=pa_auc, size=sp_area)) +
  ggtitle('AUC')

logit=function(vec, reverse=FALSE){
  if(!reverse){return(log1p(vec/(1-vec)))}
     else {return( exp(vec)/(1+exp(vec)) )}
  }


auc_lm=lm(logit(pa_auc - po_auc) ~ sp_area*present_sites, data=results)
summary(auc_lm)

kappa_lm=lm(logit(pa_kappa - po_kappa) ~ sp_area*present_sites, data=results)
summary(auc_lm)

sens_lm=lm(logit(pa_sensitivity - po_sensitivity) ~ sp_area*present_sites, data=results)
summary(sens_lm)

spec_lm=lm(logit(pa_specificity - po_specificity) ~ sp_area*present_sites, data=results)
summary(spec_lm)












