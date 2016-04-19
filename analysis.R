library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(lme4)

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

#auc
basegraph + 
  geom_point(aes(x=po_auc, y=pa_auc, size=sp_area)) +
  ggtitle('AUC')

###########################################################
#Meta analysis

logit=function(vec, reverse=FALSE){
  if(!reverse){return(log1p(vec/(1-vec)))}
     else {return( exp(vec)/(1+exp(vec)) )}
  }



metaResults = results %>%
  select(pa=pa_auc, po=po_auc, Aou, sp_area, present_sites, point_density) %>%
  gather(model_type,auc, -Aou, -sp_area, -present_sites, -point_density) %>%
  mutate(point_density=scale(point_density), sp_area=scale(sp_area))
auc_lm=lm(auc~ model_type +sp_area + point_density + as.factor(Aou) , data=metaResults)
summary(auc_lm)


metaResults = results %>%
  select(pa=pa_kappa, po=po_kappa, Aou, sp_area, present_sites, point_density) %>%
  gather(model_type, kappa , -Aou, -sp_area, -present_sites, -point_density) %>%
  mutate(point_density=scale(point_density), sp_area=scale(sp_area))
kappa_lm=lm(kappa~ model_type +sp_area + point_density + as.factor(Aou), data=metaResults)
summary(kappa_lm)

metaResults = results %>%
  select(pa=pa_sensitivity, po=po_sensitivity, Aou, sp_area, present_sites, point_density) %>%
  gather(model_type, sensitivity, -Aou, -sp_area, -present_sites, -point_density) %>%
  mutate(point_density=scale(point_density), sp_area=scale(sp_area))
sense_lm=lm(sensitivity~ model_type +sp_area + point_density + as.factor(Aou), data=metaResults)
summary(sense_lm)

metaResults = results %>%
  select(pa=pa_specificity, po=po_specificity, Aou, sp_area, present_sites, point_density) %>%
  gather(model_type, specificity, -Aou, -sp_area, -present_sites, -point_density) %>%
  mutate(point_density=scale(point_density), sp_area=scale(sp_area))
spec_lm=lm(specificity~ model_type +sp_area + point_density + as.factor(Aou), data=metaResults)
summary(spec_lm)
