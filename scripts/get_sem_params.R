# FUNCTION to process model fits
get_sem_params <- function(fit_model){
  ITS.samples <- data.frame(fixef(fit_model, summary=FALSE)) %>%
   dplyr::select(starts_with("ITS")) # get fixed effects
  colnames(ITS.samples) <- sub("ITSotus_","",colnames(ITS.samples))
  samples <- ITS.samples %>% 
    mutate(veg.AMF = vegS  
           ,veg.ECM.ERM = vegS + symbiontECM.ERM.vegS 
           ,veg.Nfix = vegS + symbiontNfix.vegS 
# calculate intercepts
           ,AMF.shrub = Intercept + vegS
           ,ECMERM.shrub = Intercept + symbiontECM.ERM + vegS + symbiontECM.ERM.vegS
           ,Nfix.shrub = Intercept+symbiontNfix + vegS + symbiontNfix.vegS
           
           ,AMF.herb = Intercept
           ,ECMERM.herb = Intercept + symbiontECM.ERM 
           ,Nfix.herb = Intercept+symbiontNfix 

# Do the shrub soil communities differ among symbiont types? (regardless of their 'effects' relative to herb communities)
           ,dif.AMF.ECMERM = AMF.shrub - ECMERM.shrub
           ,dif.AMF.Nfix = AMF.shrub - Nfix.shrub
           ,dif.ECMERM.Nfix =  ECMERM.shrub - Nfix.shrub) 
    
    
  temp <- summary(as.mcmc(samples), quantiles = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))
  fit.stats <- as.data.frame(signif(cbind(temp$statistics[,1:2], temp$quantiles),3))
  colnames(fit.stats) <- c("mean","sd","x2.5","x5","x10","x50","x90","x95","x97.5")
  p.temp <- pnorm(0,fit.stats$mean, fit.stats$sd)
  fit.stats$p <- round(sapply(p.temp, function(x) min(x, 1-x)), 5)
  
  # Get sem results
  soil.samples <- data.frame(fixef(fit_model, summary=FALSE)) %>%
    dplyr::select(starts_with("ph"), starts_with("soilCN"))

  temp <- summary(as.mcmc(soil.samples), quantiles = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))
  soil.stats <- as.data.frame(signif(cbind(temp$statistics[,1:2], temp$quantiles),3))
  colnames(soil.stats) <- c("mean","sd","x2.5","x5","x10","x50","x90","x95","x97.5")
  p.temp <- pnorm(0,soil.stats$mean, soil.stats$sd)
  soil.stats$p <- round(sapply(p.temp, function(x) min(x, 1-x)), 5)
  
  return(model.stats <- list(fit.stats, soil.stats))
}

get_sem_params_rich <- function(fit_model){
  ITS.samples <- data.frame(fixef(fit_model, summary=FALSE)) %>%
    dplyr::select(starts_with("ITS")) # get fixed effects
  colnames(ITS.samples) <- sub("ITSotus_","",colnames(ITS.samples))
  samples <- ITS.samples %>% 
    mutate(veg.AMF = vegS  
           ,veg.ECM.ERM = vegS + symbiontECM.ERM.vegS 
           ,veg.Nfix = vegS + symbiontNfix.vegS 
      # calculate intercepts
           ,AMF.shrub = Intercept + vegS
           ,ECMERM.shrub = Intercept + symbiontECM.ERM + vegS + symbiontECM.ERM.vegS
           ,Nfix.shrub = Intercept+symbiontNfix + vegS + symbiontNfix.vegS
           
           ,AMF.herb = Intercept
           ,ECMERM.herb = Intercept + symbiontECM.ERM 
           ,Nfix.herb = Intercept+symbiontNfix 
           
      # Do the shrub soil communities differ among symbiont types? (regardless of their 'effects' relative to herb communities)
           ,dif.AMF.ECMERM = AMF.shrub - ECMERM.shrub
           ,dif.AMF.Nfix = AMF.shrub - Nfix.shrub
           ,dif.ECMERM.Nfix =  ECMERM.shrub - Nfix.shrub) 
  
  
  temp <- summary(as.mcmc(samples), quantiles = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))
  fit.stats <- as.data.frame(signif(cbind(temp$statistics[,1:2], temp$quantiles),3))
  colnames(fit.stats) <- c("mean","sd","x2.5","x5","x10","x50","x90","x95","x97.5")
  p.temp <- pnorm(0,fit.stats$mean, fit.stats$sd)
  fit.stats$p <- round(sapply(p.temp, function(x) min(x, 1-x)), 5)
  
  # Get sem results
  soil.samples <- data.frame(fixef(fit_model, summary=FALSE)) %>%
    dplyr::select(starts_with("ph"), starts_with("soilCN"))
  
  temp <- summary(as.mcmc(soil.samples), quantiles = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))
  soil.stats <- as.data.frame(signif(cbind(temp$statistics[,1:2], temp$quantiles),3))
  colnames(soil.stats) <- c("mean","sd","x2.5","x5","x10","x50","x90","x95","x97.5")
  p.temp <- pnorm(0,soil.stats$mean, soil.stats$sd)
  soil.stats$p <- round(sapply(p.temp, function(x) min(x, 1-x)), 5)
  
  return(model.stats <- list(fit.stats, soil.stats))
}
