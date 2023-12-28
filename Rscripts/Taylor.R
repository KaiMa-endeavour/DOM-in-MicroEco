# @Project: DOM-in-MicroEco
# @File: Taylor.R
# @Author: Kai Ma
# Date: 13/08/2022


library(tidyfst)
library(dplyr)

data <- fread('./DOM_scaled_rarefied_int_tab.csv')
proj <- reshape2::melt(data, value.name = 'count', variable.name = 'ASVs') %>% mutate(nreads = rowSums(data[, -1])[1]) %>% filter(count >0) %>% group_by( Samples, ASVs ) %>% mutate( tf = mean(count/nreads), o = n(), tvpf = mean( (count^2 - count)/nreads^2 ) ) %>% ungroup() %>% group_by( Samples ) %>% mutate( f = o*tf, vf = o*tvpf ) #%>% mutate(vf = vf - f^2 )
gamma_pars_dom <- proj %>% select( Samples, ASVs, o, f, vf ) %>% 
  mutate( cv = sqrt(vf/f^2) ) %>% distinct() %>% 
  mutate( beta = 1./cv^2, theta = f/beta ) %>% mutate(group = 'DOM') 

data <- fread('./ASVstable_rarefied_Bacteria.csv')
proj <- reshape2::melt(data, value.name = 'count', variable.name = 'ASVs') %>% mutate(nreads = rowSums(data[, -1])[1]) %>% filter(count >0) %>% group_by( Samples, ASVs ) %>% mutate( tf = mean(count/nreads), o = n(), tvpf = mean( (count^2 - count)/nreads^2 ) ) %>% ungroup() %>% group_by( Samples ) %>% mutate( f = o*tf, vf = o*tvpf ) #%>% mutate(vf = vf - f^2 )
gamma_pars_bac <- proj %>% select( Samples, ASVs, o, f, vf ) %>% 
  mutate( cv = sqrt(vf/f^2) ) %>% distinct() %>% 
  mutate( beta = 1./cv^2, theta = f/beta ) %>% mutate(group = 'Bacteria') 

data <- fread('./ASVstable_rarefied_Fungi.csv')
proj <- reshape2::melt(data, value.name = 'count', variable.name = 'ASVs') %>% mutate(nreads = rowSums(data[, -1])[1]) %>% filter(count >0) %>% group_by( Samples, ASVs ) %>% mutate( tf = mean(count/nreads), o = n(), tvpf = mean( (count^2 - count)/nreads^2 ) ) %>% ungroup() %>% group_by( Samples ) %>% mutate( f = o*tf, vf = o*tvpf ) #%>% mutate(vf = vf - f^2 )
gamma_pars_fun <- proj %>% select( Samples, ASVs, o, f, vf ) %>% 
  mutate( cv = sqrt(vf/f^2) ) %>% distinct() %>% 
  mutate( beta = 1./cv^2, theta = f/beta ) %>% mutate(group = 'Fungi') 


gamma_pars <- rbind(gamma_pars_dom, gamma_pars_bac, gamma_pars_fun)
