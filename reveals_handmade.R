library(zipfR)
library(rioja)
setwd('~/Reveals_NEUS/')
#--------------------------------------------------------------------------------------------------------------]
#function to estimate species dependent depositional coefficient
Ki <- function(b, R, zmax=400000){
  ul1 <- b*(zmax-R)^(1/8)
  ul2 <- b*(zmax+R)^(1/8)
  ul3 <- b*(2*R)^(1/8)
  gamma_Ki <- Igamma(8, ul1, lower=TRUE) -  
    Igamma(8, ul2, lower=TRUE) +  
    Igamma(8, ul3, lower=TRUE)
  
  return(4*pi*R/b^8*gamma_Ki)
}

#---------------------------------------------------------------------------------
#function to estimate vegetation proportion
REVEALS_gamma <- function(pollen,fall_speed,ppes,R,zmax,n,u,c){
  
  #calculate parameter b later on used to estimate deposition coeffcient 
  b <- 4/sqrt(pi) * fall.speed/(n*u*c)  
  
  #species specific depositional coefficient
  K_species <- 
    sapply(b,function(x){
      Ki(x[[1]],R=R,zmax = zmax)
    })
  
  #eq (5) in Sugita (2007a)
  weighted_pollen <- pollen/(ppes*K_species)
  veg_proportion <- weighted_pollen/sum(weighted_pollen)  
  return(as.numeric(veg_proportion))
}


#------------------------------------------------------------------------------------------------------------
#load data
pollen <- read.csv("data/reveals_input.csv")
pollen <- pollen[-1]
names(pollen) <- unlist(strsplit(names(pollen),'[.]'))[seq(2,(2*ncol(pollen)),2)] 


reveals.params <- read.csv('data/reveals_input_params_variable.csv')
taxa <- reveals.params$species

#fall speed
fall.speed <- reveals.params$fallspeed
names(fall.speed) <- taxa
fall.speed <- as.data.frame(t(fall.speed))


#ppe
ppes <- reveals.params$PPEs
names(ppes) <- taxa

fall.speed <- fall.speed[names(fall.speed)%in%names(pollen)]
ppes <- ppes[names(ppes)%in%names(pollen)] 

#massive difference depending on use of meters or km
REVEALS_gamma(pollen,fall_speed,ppes,R=1,zmax=100,n=0.25,u=3,c=0.12)
REVEALS_gamma(pollen,fall_speed,ppes,R=1000,zmax=100000,n=0.25,u=3,c=0.12)


#----------------------------------------------------------------------------------------------------
#test sensitivity to maximum dispersal
distances <- c(seq(10,100,10),seq(200,1000,100))

sensitivity_reveals_gamma <- 
  sapply(distances,function(x){
    REVEALS_gamma(pollen,fall_speed,ppes,R=1,zmax=x,n=0.25,u=3,c=0.12)
  })    

colnames(sensitivity_reveals_gamma) <- distances


#distances
dist(t(sqrt(sensitivity_reveals_gamma[,c("50","100","400")])))^2
paldist(t(sensitivity_reveals_gamma[,c("50","100","400")]))

#-----------------------------------------------------------------------------------------------------
# load new pollen data
#-----------------------------------------------------------------------------------------------------
load('~/workflow_stepps_calibration/vegetation/data_nb/prediction_13_taxa_6796_cells_120_knots_cal_pl_Ka_Kgamma_EPs_79_sites_final.rdata')
pollen <- y  
colnames(pollen)[grep("Other",colnames(pollen))] <- c('Other_conifer','Other_hardwood') 
pollen <- pollen[,colnames(pollen)%in%names(ppes)]

sensitivity_all <-   
lapply(distances,function(distan){
  pred_reveals <- 
    apply(pollen,1,function(x){
      REVEALS_gamma(x,fall_speed,ppes,R=1,zmax=distan,n=0.25,u=3,c=0.12)
    })
  pred_reveals <- t(pred_reveals)
  colnames(pred_reveals) <- colnames(pollen)
  round(pred_reveals,3)
})

names(sensitivity_all) <- distances
