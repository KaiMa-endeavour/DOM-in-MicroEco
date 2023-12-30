# @Project: DOM-in-MicroEco
# @File: Trans_turnover.R
# @Author: Kai Ma
# Date: 14/05/2023


library(tidyfst)

trans <- fread('./molecular_transformations_profile.csv') %>% arrange_dt(Samples)

geo <- fread('./geo_site.csv')
geo_exp <- lapply(1:nrow(geo), function(r) {
  data.table(Station = rep(geo$Station[r], sum(grepl(geo$Station[r], unlist(trans[, 1])))), 
             Longitude = rep(geo$Longitude[r], sum(grepl(geo$Station[r], unlist(trans[, 1])))), 
             Latitude = rep(geo$Latitude[r], sum(grepl(geo$Station[r], unlist(trans[, 1]))))      )
}) %>% rbindlist() %>% arrange_dt(Station)

lat_abun <- cbind(geo_exp[, c(1, 3)], trans[, -1]) %>% reshape2::melt(id.vars = c('Station', 'Latitude'))

captrue_trans <- names(table(lat_abun$variable))
slope_t <- matrix(NA, length(captrue_trans), 6)

for (i in 1:length(captrue_trans)) {
  model <- lm(value ~ Latitude , data = filter_dt(lat_abun, variable == captrue_trans[i]))
  summ <- summary(model)
  res <- summ$coefficients
  if (res[2, 4] <= 0.05) {
    message(captrue_trans[i])
    mse <- sum(summ$residuals**2)/length(summ$residuals)
    if (res[2, 4] <= 0.001) {
      star <- '***'
    } else if (res[2, 4] <= 0.01) {
      star <- '**'
    } else if (res[2, 4] <= 0.05) {
      star <- '*'
    } else if (res[2, 4] <= 0.1) {
      star <- '.'
    } else {
      star <- ''
    }
    slope_t[i, ] <- c(captrue_trans[i], res[2, 1], res[2, 3], res[2, 4], star, mse)
  }
}
colnames(slope_t) <- c('Name', 'slope', 't_value', 'P_value', 'star', 'MSE')
slope_t <- as.data.frame(slope_t[!is.na(slope_t[, 1]), ])
bio <- fread('./Biotic-abiotic-transfromation-classification.csv')
slope_t_bio <- merge(slope_t, bio)
slope_t_bio$slope <- as.numeric(slope_t_bio$slope)
slope_t_bio$MSE <- as.numeric(slope_t_bio$MSE)
slope_t_bio <- arrange_dt(slope_t_bio, slope)
slope_t_bio
nrow(slope_t_bio)
fwrite(slope_t_bio, 'trans-slope-P.csv')
