# @Project: DOM-in-MicroEco
# @File: LDG.R
# @Author: Kai Ma
# Date: 24/04/2022


library(tidyfst)

dom <- fread('./DOM_scaled_rarefied_int_tab.csv') %>% arrange_dt(Samples)
bac <- fread('./ASVstable_rarefied_Bacteria.csv') %>% arrange_dt(Samples)
fun <- fread('./ASVstable_rarefied_Fungi.csv') %>% arrange_dt(Samples)

geo <- fread('./geo_site.csv')

lat_div <- function(data, sample_names) {
  AOR <- taxa_partition(data, 'frequency')
  AOR$whole_community <- data
  
  res <- lapply(1:length(AOR), function(i) {
    sample_id <- which(rowSums(select_dt(data, cols = colnames(AOR[[i]]))) > 0)
    geo_exp <- lapply(1:nrow(geo), function(r) {
      data.table(station = rep(geo$Station[r], sum(grepl(geo$Station[r], sample_names[sample_id]))),
                 Longitude = rep(geo$Longitude[r], sum(grepl(geo$Station[r], sample_names[sample_id]))),
                 Latitude = rep(geo$Latitude[r], sum(grepl(geo$Station[r], sample_names[sample_id])))      )
    }) %>% rbindlist() %>% arrange_dt(station)
    
    data.table(alphaDiversity(comm = AOR[[i]], sample_names = sample_names[sample_id], method = c('Richness', 'Chao1')), Taxa = c('Abundant', 'Occasional', 'Rare', 'Whole')[i], lat = geo_exp$Latitude)
  }) %>% rbindlist()
  res
}

latdiv_dom <- lat_div(dom[, -1], dom$Samples) %>% mutate_dt(group = 'DOM')
latdiv_bac <- lat_div(bac[, -1], bac$Samples) %>% mutate_dt(group = 'Bacteria')
latdiv_fun <- lat_div(fun[, -1], fun$Samples) %>% mutate_dt(group = 'Fungi')
