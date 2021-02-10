# Gurue Paper - Agricultural expansion in African savannas – dissimilar effects on diversity and composition of trees and mammals
# hgtripathi05@gmail, Script - 2021
# Hemant G. Tripathi com 
# January 2021

# 1. Prepare workspace ----------
# install jags software from http://www.sourceforge.net/projects/mcmc-jags/files
# rm(list = ls(all=TRUE)) # clean up space
.libPaths("E:/rlib") # Set library
#setwd(choose.dir(getwd())) # set working directory to where .RData file is
setwd("E:/res_projs/gurue_frag/gurue_analysis/biodiversity_conservation/revision_20sept/analysis/mldi3_trldi2")

# 2. Load Packages and functions ----------
packages<-(c("plyr","dplyr", "vegan", "adespatial", "ape", "CommEcol",
             "FactoMineR","picante","FD","extrafont","extrafontdb" , "psych",  "jagsUI", "coda",
             "ggplot2", "gridExtra", "ade4", "fpc", "party")) # List of packages required for analysis
lapply(packages, library, character.only = TRUE) # Load packages

loadfonts(quiet = T)
antilogit <- function(x) { exp(x) / (1 + exp(x) ) }

# 3. Analysis -  Mammal Communuities ----------
# 3.1 Prepare data -------
m<-read.csv("mam.csv", header=TRUE); head(m) # mammal dataset 
m$engname<-as.factor(as.character(m$engname)); m$sciname<-as.factor(as.character(m$sciname))

(spec.name.list <- tapply(m$spid, m$sciname, mean)) # species ID
(ordered.spec.name.list <- spec.name.list[order(spec.name.list)]) # ID-order list
DET<- as.matrix(m[,17:86]); DET[DET > 1] <- 1; sum(DET) # should be 986  #detection
nsite<-length(unique(m$grid_code)) # number of camera trap grids
nspec<-length(unique(m$sciname)) # number of species 
nrep<-dim(DET)[2] # number of day replicates 

Y <- array(NA, dim = c(nsite, nrep, nspec)) # Put detection data into 3D array: site x rep x species
for(i in 1:nspec){
  Y[,,i] <- DET[((i-1)*nsite+1):(i*nsite),]
}

# get mammal grid covariates 
orig.ldi <- m$ldi[1:nsite] # land division index
orig.hl<-1-m$loss_intensity[1:nsite] # intensity of change between 07 and 14
orig.me <- m$miombo_extent[1:nsite]# total woodland area
orig.map <- m$MAP[1:nsite] # mean annual precipitation

# standardise
ldi <- as.numeric(scale(orig.ldi, center = TRUE, scale = TRUE)) # centered ldi
hl<-as.numeric(scale(orig.hl, center=TRUE, scale=TRUE)) # cenetered int 
me <- as.numeric(scale(orig.me, center = TRUE, scale = TRUE)) # centered wc
map<-as.numeric(scale(orig.map, center=TRUE, scale=TRUE)) # cenetered map 

# extract residuals of predictors - importance sequence ldi  + wc + int 
me.res<-summary(lm(me~ldi))$residuals 
hl.res<-summary(lm(hl~ldi+me.res))$residuals 

# 3.2 Mammal Model 1 - Mammal occupancy model -------
# DR community model with data augmentation

nz <- 30 - nspec                  # Number of potential species in superpopulation
M <- nspec + nz          # Size of augmented data set ('superpopulation')
Ysum <- apply(Y, c(1,3), sum, na.rm = T) # Collapse to detection frequency
Yaug <- cbind(Ysum, array(0, dim=c(nsite, nz))) # Add all zero histories
dim(Yaug)

# Bundle and summarize data set
str( win.data <- list(Yaug = Yaug, nsite = nrow(Ysum), nrep = m$camdays[1:nsite]+1, M = M, nspec = nspec, nz = nz, 
                      var1=ldi, var2=me.res, var3=hl.res, var4=map)) 
# Initial values
wst <- rep(1, nspec+nz)                   # Simply set everybody at 'occurring'
zst <- array(1, dim = c(nsite, nspec+nz)) # ditto for z
inits <- function() list(z = zst, w = wst, lpsi = rnorm(n = nspec+nz), lp = rnorm(n = nspec+nz))

# Parameters monitored
params <- c("mu.lpsi", "mu.betalpsi", "sd.betalpsi", "sd.lpsi", "betalpsi", 
            "lpsi", "mu.lp", "sd.lp", "p", "Nsite", "Ntotal", "omega", "n0")
ni <- 75000   ;   nt <- 50  ;   nb <- 25000   ;   nc <- 3 # MCMC settings

# part 1 of model 1 
#mmod11 <- jags(win.data, inits, params, "mmod1.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb) ; save.image()
# part 2 of model 1 
#params2 <- c("z") # detection corrected presence 
#mmod12<- jags.basic(win.data, inits, params2,"mmod1.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
#mmod12.1 <- as.matrix(mmod12) ; save.image()# Put output from 3 chains into a matrix

# 3.3 Mammal Model 2 - Mammal Diversity Cubic Regression models -------
# Extract diversity variables
pm<-mmod11$mean$Nsite # Species Richness
psd<-mmod11$sd$Nsite # posterior sd's of species richness estimates 
cri<-cbind(mmod11$q2.5$Nsite, mmod11$q97.5$Nsite) # CRL's of species richness
SR.pm<-as.data.frame(cbind('SR.pm' = pm, 'SR.sd' = psd, 'SR-2.5%' = cri[,1], 'SR-97.5%' = cri[,2])) 

# Species Beta Diversity 
nsamp<-dim(mmod12.1)[1] # Sps matrix 
Z<-array(NA, dim = c(nsite, nspec+nz, nsamp))
for (j in 1: nsamp) { 
  cat (paste ("nMCMC sample", j, "\n"))
  Z[,,j]<- mmod12.1[j, 2:1111]}
N <- Z[,1:nspec,]; dim(N) # Restrict computations to observed species

# Turnover component of Baselga's Multipart Beta Diversity
MCD <- array(NA, dim = c(nsite, nsamp)) # array for mean distance between communities in composition of species 
for(k in 1:nsamp){
  print(k)
  MCD[,k]<-apply((as.matrix(beta.div.comp(N[,,k], coef = "BS", quant = F)$repl)), MARGIN = 1, mean)
} 

pm <- apply(MCD, 1, mean, na.rm = TRUE) # Post. mean of Jsite wrt. site 1
psd <- apply(MCD, 1, sd, na.rm = TRUE) # Post. sd of Jsite wrt. site 1
cri <- apply(MCD, 1, function(x) quantile(x, prob = c(0.025, 0.975), na.rm =TRUE)) # CRI
MCD.pm<-as.data.frame(cbind('MCD.pm' = pm, 'MCD.sd' = psd, 'MCD-2.5%' = cri[1,], 'MCD-97.5%' = cri[2,])) 

# Nestedness component of Baselga's Multipart Beta Diversity
MND <- array(NA, dim = c(nsite, nsamp)) # array for mean distance between communities in composition of species 
for(k in 1:nsamp){
  print(k)
  MND[,k]<-apply((as.matrix(beta.div.comp(N[,,k], coef = "BS", quant = F)$rich)), MARGIN = 1, mean)
} 

pm <- apply(MND, 1, mean, na.rm = TRUE) # Post. mean of Jsite wrt. site 1
psd <- apply(MND, 1, sd, na.rm = TRUE) # Post. sd of Jsite wrt. site 1
cri <- apply(MND, 1, function(x) quantile(x, prob = c(0.025, 0.975), na.rm =TRUE)) # CRI
MND.pm<-as.data.frame(cbind('MND.pm' = pm, 'MND.sd' = psd, 'MND-2.5%' = cri[1,], 'MND-97.5%' = cri[2,])) 

com.pm<-as.data.frame(cbind(SR.pm["SR.pm"], MCD.pm["MCD.pm"], MND.pm["MND.pm"]))
com.sd<-as.data.frame(cbind(SR.pm["SR.sd"], MCD.pm["MCD.sd"], MND.pm["MND.sd"]))

# Bundle and summarize data set
pred.ldi<- (seq(0.001, 0.99,,500)-mean(orig.ldi))/sd(orig.ldi)
str(win.data <- list(pm=as.matrix(com.pm), psd=as.matrix(com.sd) ,n = dim(com.pm)[1], 
                     var1=ldi, var2=me.res, var3=hl.res, var4=map, nv=dim(com.pm)[2], 
                     pred.ldi = pred.ldi,npred = length(pred.ldi)))
inits <- function() list(beta0=runif(dim(com.pm)[2], 0,1)) # initial values for chain start 
params <- c("beta0", "beta", "Npred")
#mmod2com<- jags(win.data, inits, params, "mmod2com.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb) 

# 3.4 Mammal Results -------

# Detection corrected estimates of observed diversity 
# 3.4.1 Mammal Species richness  -----
obsldi<-as.data.frame(cbind(x=orig.ldi, y=SR.pm$SR.pm, ymin = SR.pm$`SR-2.5%`, ymax =  SR.pm$`SR-97.5%`))
ldisr<-as.data.frame(cbind(x=seq(0, 0.99,,500), y= mmod2com$mean$Npred[1,], ymin = mmod2com$q2.5$Npred[1,], ymax = mmod2com$q97.5$Npred[1,]))


# Thresholds Hooper et al 
base<-max(mmod2com$mean$Npred[1,])
band1<-as.data.frame(cbind(x=seq(0, 0.99,,500),y=base,  ymin =0.60*base, ymax =0.79*base))
band1$thtype<-"Intermediate"
band2<-as.data.frame(cbind(x=seq(0, 0.99,,500),y=base,  ymin =0.40*base, ymax =0.59*base))
band2$thtype<-"High"
th<-rbind(band2, band1)
th$thtype<-factor(th$thtype, levels=c('Intermediate', "High"))

pred.sr<-ldisr
pred.sr$y1 <- ((pred.sr$y/base)-1)*100
pred.sr$ymin1 <- ((pred.sr$ymin/base)-1)*100
pred.sr$ymax1 <- ((pred.sr$ymax/base)-1)*100
pred.sr$restype<-"Richness"

fig5b_mamsr_ldi<-ggplot(data = ldisr, aes(y = y, x = x, ymin = ymin, ymax = ymax)) +
  geom_line(colour="darkblue", size=0.75, alpha=0.4) +
  geom_ribbon(alpha = 0.175, fill = "blue") + 
  geom_point(data = obsldi, aes(x = x, y = y), size=4, alpha=0.4, shape=19, colour="grey70")+ 
  geom_linerange(data=obsldi, aes(ymin=ymin, ymax=ymax), size=.5,  colour="grey70", alpha=0.4) + 
  scale_x_continuous(limits = c(0, 0.955), breaks=seq(0,1,0.25), "Land Division Index") + 
  theme(legend.position="bottom",legend.title=element_text(colour="gray25", size = 11, family = "Calibri"), 
        legend.text=element_text(colour="gray30", size = 11, family = "Calibri"),
        axis.text=element_text(colour="grey50", size = rel(0.9), family = "Calibri"), 
        axis.title.x = element_text(colour = "grey40", family = "Calibri", size=12), 
        axis.title.y = element_text(colour = "grey40", family = "Calibri", size=12),panel.background = element_blank()) + 
  theme(axis.line = element_line(color = 'grey75', size=0.2)) + 
  scale_y_continuous(limits = c(-3.5, 16), "Estimated species richness of mammals", 
                     sec.axis = sec_axis(~ ((./base)-1)*100 , name = "% reduction in species richness", breaks = c(-75, -50, -25, 0)))+ 
  geom_ribbon(data=th, aes(colour = factor(thtype), fill = factor(thtype)),                                                                                                        linetype = 3, size=0.1, alpha= 0.125) + 
  scale_fill_manual(values = c('orange','red')) +
  guides(colour=FALSE, fill = guide_legend(override.aes = list(size=4))) +  
  labs(fill = "Levels of species loss (Hooper et al. 2012)") + 
  geom_smooth(data = obsldi, method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE, colour="grey50", linetype=2, alpha=0.75, size=0.35)
fig5b_mamsr_ldi
ggsave(plot=fig5b_mamsr_ldi, "fig5b_mamsr_ldi.png",width = 5, height = 6) 

# 3.4.2 Mammal Occupancy and Detection  -----
par(mfrow = c(1,1))
psi.sample <- plogis(rnorm(10^6, mean = mmod11$mean$mu.lpsi, sd = mmod11$mean$sd.lpsi))
p.sample <- plogis(rnorm(10^6, mean = mmod11$mean$mu.lp, sd = mmod11$mean$sd.lp))
summary(psi.sample)   ;   summary(p.sample) # check psi and p
samp<-as.data.frame(cbind(psi=psi.sample, p=p.sample))
psi.hist<-ggplot(samp, aes(x=psi)) + geom_histogram(aes(y=..density..), col="white", bins = 50, fill="grey30", alpha=0.75) +
  ylab("Density") + xlab( expression(Occupancy~(psi))) + 
  theme(axis.text=element_text(colour="gray40", size = 11, family = "Calibri"), 
        axis.title.x = element_text(colour = "grey30", family = "Calibri", size=11), 
        axis.title.y = element_text(colour = "grey30", family = "Calibri", size=11),
        axis.line.x = element_line(color="grey30", size = 0.2),
        axis.line.y = element_line(color="grey30", size = 0.2),
        strip.text.x = element_text(colour = "gray20", family = "Calibri", size=11),  panel.background = element_blank())

p.hist<-ggplot(samp, aes(x=p)) + geom_histogram(aes(y=..density..), col="white", bins = 50, fill="grey30", alpha=0.75) +
  ylab("") + xlab( expression(Detection~(p))) + 
  theme(axis.text=element_text(colour="gray40", size = 11, family = "Calibri"), 
        axis.title.x = element_text(colour = "grey30", family = "Calibri", size=11), 
        axis.title.y = element_text(colour = "grey30", family = "Calibri", size=11),
        axis.line.x = element_line(color="grey30", size = 0.2),
        axis.line.y = element_line(color="grey30", size = 0.2),
        strip.text.x = element_text(colour = "gray20", family = "Calibri", size=11),  panel.background = element_blank())
occuhist<-grid.arrange(psi.hist, p.hist, ncol=2)
ggsave(plot=occuhist, "figure3_community_occu_det_hist.png", width = 4.9, height = 2.75) 

# 3.4.3 Mammal Species Level Occupancy and Detection  -----
# Occupancy
mam.occu<-as.data.frame(cbind(mmod11$mean$lpsi))
colnames(mam.occu)<-"occu"
mam.occu$c1<-mmod11$q2.5$lpsi
mam.occu$c2<-mmod11$q97.5$lpsi
mam.occu<-mam.occu[1:nspec,] # select only observed species 
mam.occu$resname<-rownames(ordered.spec.name.list)
mam.occu$occu1<-antilogit(mam.occu$occu)
mam.occu$c11<-antilogit(mam.occu$c1)
mam.occu$c21<-antilogit(mam.occu$c2)

# Plot of mammal occupancy 
mamsps.occu<-ggplot(mam.occu, aes(x = reorder(resname, occu1) , y = occu1, ymin=c11, ymax = c21)) + 
  geom_point(aes(size =occu), alpha = 0.5, shape = 19, colour = "black") + scale_size_continuous(range = c(2, 5)) + 
  geom_linerange(size=0.25, alpha = 0.75, color="grey25")+ coord_flip()+ 
  ylab("Occupancy probability") + xlab("Mammal species") + 
  theme(legend.position="none", legend.title=element_blank(), 
        legend.text=element_text(colour="gray30", size = 10, family = "Calibri"), 
        axis.text=element_text(colour="grey50", size = rel(0.85), family = "Calibri", face="italic"), 
        axis.title.x = element_text(colour = "grey30", family = "Calibri", size=13), 
        axis.ticks.y=element_blank(),
        axis.title.y = element_text(colour = "grey30", family = "Calibri", size=13),
        strip.text.x = element_text(colour = "grey40", family = "Calibri"),  panel.background = element_blank())
mamsps.occu
ggsave(plot=mamsps.occu, "mam_sps_ocu.png", width = 4.2, height = 4.75)  

# Plot of mammal Detection probablity
mam.det<-as.data.frame(cbind(mmod11$mean$p)); colnames(mam.det)<-"det"
mam.det$c1<-mmod11$q2.5$p
mam.det$c2<-mmod11$q97.5$p
mam.det<-mam.det[1:nspec,] # select only observed species 
mam.det$resname<-rownames(ordered.spec.name.list)
mamsps.det<-ggplot(mam.det, aes(x = reorder(resname, det) , y = det, ymin=c1, ymax = c2)) + 
  geom_point(aes(size =det), alpha = 0.5, shape = 19, colour = "black") + scale_size_continuous(range = c(2, 5)) + 
  geom_linerange(size=0.25, alpha = 0.75, color="grey25")+ coord_flip()+ 
  ylab("Detection probability") + xlab("Mammal species") + 
  theme(legend.position="none", legend.title=element_blank(), 
        legend.text=element_text(colour="gray30", size = 10, family = "Calibri"), 
        axis.text=element_text(colour="grey50", size = rel(0.85), family = "Calibri", face="italic"), 
        axis.title.x = element_text(colour = "grey30", family = "Calibri", size=13), 
        axis.ticks.y=element_blank(),
        axis.title.y = element_text(colour = "grey30", family = "Calibri", size=13),
        strip.text.x = element_text(colour = "grey40", family = "Calibri"),  panel.background = element_blank()); mamsps.det
ggsave(plot=mamsps.det, "mam_sps_det.png", width = 4.2, height = 4.75) 

# 3.4.4 Extract predicted species richness, turnover and nestedness of Mammals -----
# Richness
pred.sr<-as.data.frame(cbind(x=seq(0, 0.99,,500), y= mmod2com$mean$Npred[1,], ymin = mmod2com$q2.5$Npred[1,], ymax = mmod2com$q97.5$Npred[1,]))
base<-max(mmod2com$mean$Npred[1,])
pred.sr$y1 <- ((pred.sr$y/base)-1)*100
pred.sr$ymin1 <- ((pred.sr$ymin/base)-1)*100
pred.sr$ymax1 <- ((pred.sr$ymax/base)-1)*100
pred.sr$restype<-"Richness"
head(pred.sr)

# Turnover
pred.turn<-as.data.frame(cbind(x=seq(0, 0.99,,500), y= mmod2com$mean$Npred[2,], ymin = mmod2com$q2.5$Npred[2,], ymax = mmod2com$q97.5$Npred[2,]))
base<-max(mmod2com$mean$Npred[2,])
pred.turn$y1 <- ((pred.turn$y/base)-1)*100
pred.turn$ymin1 <- ((pred.turn$ymin/base)-1)*100
pred.turn$ymax1 <- ((pred.turn$ymax/base)-1)*100
pred.turn$restype<-"Turnover"
head(pred.turn)

# Nestedness
pred.nest<-as.data.frame(cbind(x=seq(0, 0.99,,500), y= mmod2com$mean$Npred[3,], 
                               ymin = mmod2com$q2.5$Npred[3,], ymax = mmod2com$q97.5$Npred[3,]))
base<-min(mmod2com$mean$Npred[3,])
pred.nest$y1 <- ((pred.nest$y/base)-1)*100
pred.nest$ymin1 <- ((pred.nest$ymin/base)-1)*100
pred.nest$ymax1 <- ((pred.nest$ymax/base)-1)*100
pred.nest$restype<-"Nestedness"
head(pred.nest)

pred.mam<-as.data.frame(rbind(pred.sr, pred.turn, pred.nest))
pred.mam$taxa<-"Mammals"
head(pred.mam)


# Observed raw values 

obs.sr<-as.data.frame(cbind(x=orig.ldi, y=SR.pm$SR.pm, ymin = SR.pm$`SR-2.5%`, ymax =  SR.pm$`SR-97.5%`)); head(obssr)
obs.sr$restype<-"Richness"

obs.turn<-as.data.frame(cbind(x=orig.ldi, y=MCD.pm$MCD.pm, ymin = MCD.pm$`MCD-2.5%`, ymax =  MCD.pm$`MCD-97.5%`)); head(obs.turn) 
obs.turn$restype<-"Turnover"

obs.nest<-as.data.frame(cbind(x=orig.ldi, y=MND.pm$MND.pm, ymin = MND.pm$`MND-2.5%`, ymax =  MND.pm$`MND-97.5%`)); head(obs.nest) 
obs.nest$restype<-"Nestedness"

obs.mam<-as.data.frame(rbind(obs.sr, obs.turn, obs.nest))
obs.mam$taxa<-"Mammals"
head(obs.mam)


# 3.4.4 Extract Coefficent table for mammals -----
# Model 1 - Occupancy
slope.w<-as.data.frame(cbind(mmod11$mean$betalpsi))
slope.w$intercept<-mmod11$mean$lpsi # Intercept 
slope.w$intercept.sd<-mmod11$sd$lpsi # sd of intercept 

pred.var<-colnames(slope.w)[1:6] # number predictor veriables
pred.names<-c("LDI", "LDI2","LDI3", "ME.resid", "HL.resid", "MAP")

slope<-reshape(slope.w, varying = pred.var, v.names = "slope",timevar = "predictor",times = pred.names, direction = "long")
prob.pos<-apply(ifelse(mmod11$sims.list$betalpsi>0, 1, 0),c(2,3), function (x) sum(x)/3000) # probablity of +ive
slope$prob.pos<-reshape(as.data.frame(prob.pos),varying = pred.var, v.names = "prob.pos", direction = "long")[,2]
prob.neg<-apply(ifelse(mmod11$sims.list$betalpsi<0, 1, 0),c(2,3), function (x) sum(x)/3000) # probablity of -ive
slope$prob.neg<-reshape(as.data.frame(prob.neg), varying = pred.var,v.names = "prob.neg", direction = "long")[,2]
slope$c1<-reshape(as.data.frame(cbind(mmod11$q2.5$betalpsi)),varying = pred.var, v.names = "c1",  direction = "long")[,2]
slope$c2<-reshape(as.data.frame(cbind(mmod11$q97.5$betalpsi)), varying = pred.var, v.names = "c2", direction = "long")[,2]
slope$sd<-reshape(as.data.frame(cbind(mmod11$sd$betalpsi)), varying = pred.var, v.names = "sd", direction = "long")[,2]
mod1.res<-slope

mod1.res$model<-rep("mod1", length(mod1.res$slope))
mod1.res<-mod1.res[mod1.res$id<=length(unique(m$sciname)),] # select only observed species 
mod1.res$res.var<-rep("pa", length(mod1.res$slope)) 
mod1.res$res.sd<-1  
mod1.res$resname<-rownames(ordered.spec.name.list)
head(mod1.res)
dim(mod1.res) # should be 21* 6 *14

# Community mean coefficients from model 1
slope.w<-as.data.frame(cbind(t(mmod11$mean$mu.betalpsi)))
slope.w$intercept<-mmod11$mean$mu.lpsi # 1st level  of Alpha1  ie. boom/ intensity class preboom is the intercept 
slope.w$intercept.sd<-mmod11$sd$mu.lpsi 
slope<-reshape(slope.w, varying = pred.var, v.names = "slope",timevar = "predictor", times = pred.names, direction = "long")
slope$prob.pos<-apply(ifelse(mmod11$sims.list$mu.betalpsi>0, 1, 0),2, function (x) sum(x)/3000) # probablity of +ive
slope$prob.neg<-apply(ifelse(mmod11$sims.list$mu.betalpsi<0, 1, 0),2, function (x) sum(x)/3000) # probablity of +ive
slope$c1<-t(mmod11$q2.5$mu.betalpsi)[1,]
slope$c2<-t(mmod11$q97.5$mu.betalpsi)[1,]
slope$sd<-t(mmod11$sd$mu.betalpsi)[1,]

com.mean.res<-slope
com.mean.res$model<-rep("com.mean", length(com.mean.res$slope))
com.mean.res$res.var<-rep("pa", length(com.mean.res$slope)) 
com.mean.res$res.sd<-1; com.mean.res$resname<-NA; com.mean.res[is.na(com.mean.res)]<-0
dim(com.mean.res); head(com.mean.res)

# model 2 - diversity model coefficients
slope.w<-as.data.frame(cbind(mmod2com$mean$beta))
slope.w$intercept<-mmod2com$mean$beta0 # 1st level  of Alpha1  ie. boom/ intensity class preboom is the intercept 
slope.w$intercept.sd<-mmod2com$sd$beta0
slope<-reshape(slope.w, varying = pred.var, v.names = "slope",timevar = "predictor",times = pred.names, direction = "long")

prob.pos<-apply(ifelse(mmod2com$sims.list$beta>0, 1, 0),c(2,3), function (x) sum(x)/3000) # probablity of +ive
slope$prob.pos<-reshape(as.data.frame(prob.pos),varying = pred.var, v.names = "prob.pos", direction = "long")[,2]

# probbality of negative effect 
prob.neg<-apply(ifelse(mmod2com$sims.list$beta<0, 1, 0),c(2,3), function (x) sum(x)/3000) # probablity of -ive

slope$prob.neg<-reshape(as.data.frame(prob.neg),varying = pred.var,v.names = "prob.neg", direction = "long")[,2]
slope$c1<-reshape(as.data.frame(cbind(mmod2com$q2.5$beta)),varying = pred.var, v.names = "c1",  direction = "long")[,2]
slope$c2<-reshape(as.data.frame(cbind(mmod2com$q97.5$beta)),varying = pred.var, v.names = "c2", direction = "long")[,2]
slope$sd<-reshape(as.data.frame(cbind(mmod2com$sd$beta)),varying = pred.var, v.names = "sd", direction = "long")[,2]
mod2.res<-slope
mod2.res$model<-rep("mod2", length(mod2.res$slope))
mod2.res$res.var<-rep("Com", length(mod2.res$slope)) 
mod2.res$res.sd<-apply(com.pm, 2,  sd)
mod2.res$resname<-c("SR", "MCD", "MND")
dim(mod2.res)
mmod.res<-rbind(mod1.res, mod2.res, com.mean.res) # all mammal model results 
mmod.res$species_group<-"Mammals"
dim(mmod.res); save.image()

# 4. Analysis -  Tree Communuities ----------
# 4.1 Prepare data -------
tr<-read.csv("tr.csv", header=TRUE); head(tr) # tree abundance dataset
clus<-c("G497", "G543", "G585", "G627", "G628", "G673", "G611", "G654", "G653", "G609","G652", "G651") 
tr1<-tr[!tr$grid_code %in% clus,]
tr1$grid_code<-factor(tr1$grid_code)
tr<-tr1

# remove cluster of trees
counts <- as.matrix(tr[,13:17]); head(counts) ; sum(counts, na.rm = T) # should be 1885 # abundance
nsite<-length(unique(tr$grid_code)); nsite # number of grids
nspec<-length(unique(tr$sciname)); nspec # number of species
nrep<-dim(counts)[2]; nrep # number of plot replicates

Yc <- array(NA, dim = c(nsite, nrep, nspec)) # Put detection data into 3D array: site x rep x species
for(i in 1:nspec){
  Yc[,,i] <- counts[((i-1)*nsite+1):(i*nsite),]  # 'c' for counts
}
sum(Yc, na.rm = T) 
# get tree grid covariates 
orig.ldi <- tr$ldi[1:nsite] # land division index
orig.hl<-1-tr$loss_intensity[1:nsite] # intensity of change between 07 and 14
orig.me <- tr$miombo_extent[1:nsite]# total woodland area
orig.map <- tr$MAP[1:nsite] # mean annual precipitation

# standardise
ldi <- as.numeric(scale(orig.ldi, center = TRUE, scale = TRUE)) # centered ldi
hl<-as.numeric(scale(orig.hl, center=TRUE, scale=TRUE)) # cenetered int 
me <- as.numeric(scale(orig.me, center = TRUE, scale = TRUE)) # centered wc
map<-as.numeric(scale(orig.map, center=TRUE, scale=TRUE)) # cenetered map 

# extract residuals of predictors - importance sequence ldi  + wc + int 
me.res<-summary(lm(me~ldi))$residuals 
hl.res<-summary(lm(hl~ldi+me.res))$residuals 

# 4.2 Tree Model 1 - Abundance -------
nz <- 150 - nspec         # Use for informative prior, we expected about 150 species in our study area
Yaug <- array(0, dim=c(nsite, nrep, nspec+nz)) # array with only zeroes
Yaug[,,1:nspec] <- Yc      # copy into it the observed data

# Bundle and summarize data set
str(win.data <- list(Yc = Yaug, nsite = dim(Yc)[1], nrep = dim(Yc)[2],nspec = dim(Yc)[3],
                     M=nspec+nz, var1=ldi, var2=me.res, var3=hl.res, var4=map))

wst <- rep(1, nspec+nz) # Initial values
some.more <- 1       # May have to play with this until JAGS is happy
Nst <- apply(Yaug, c(1,3), max, na.rm = T) + some.more
Nst[Nst == '-Inf'] <- 10          # May have to play with this, too
Nst <- Nst
inits <- function()list(N = Nst, w=wst)

# Parameters monitored
params <- c("omega", "mlambda",  
            "det0", "mu.det0", "sd.det0",  
            "beta0", "mu.beta0", "sd.beta0",
            "beta",  "mu.beta", "sd.beta",
            "Nsite", "Ntotal", "n0")
# part 1 of model 1 
#trmod1C<- jags(win.data, inits, params, "trmod1C.txt", n.chains = nc,
               n.thin = nt, n.iter = ni, n.burnin = nb)

# part 2 of model 1 
#params2 <- c( "N")
#trmodC12 <- jags.basic(win.data, inits, params2, "trmod1C.txt", n.chains = nc, 
                       n.thin = nt, n.iter = ni, n.burnin = nb)
#trmodC12.1<-as.matrix(trmodC12); save.image()

# 4.3 Tree Model 2 - Tree Diversity Cubic Regression models -------
pm<-trmod1C$mean$Nsite # Species richness estimates 
psd<-trmod1C$sd$Nsite # posterior sd's of species richness estimates 
cri<-cbind(trmod1C$q2.5$Nsite, trmod1C$q97.5$Nsite) # CRL's of species richness
SR.pm<-as.data.frame(cbind('SR.pm' = pm, 'SR.sd' = psd, 'SR-2.5%' = cri[,1], 'SR-97.5%' = cri[,2])) 

# Species Beta Diversity 
nsamp<-dim(trmodC12.1)[1]
N<-array(NA, dim = c(nsite, nspec+nz, nsamp))
for (j in 1: nsamp) { 
  cat (paste ("nMCMC sample", j, "\n"))
  N[,,j]<- trmodC12.1[j, 1:(4050)]}

#5850 if you have used all species 

# Restrict computations to observed species
N <- N[,1:nspec,] # Species 1 to 98, keep only observed species 
dim(N) # check, should be 135 X 98 X 3000 

# Species Beta Diversity 
# Turnover component of Baselga's Multipart Beta Diversity
MCD <- array(NA, dim = c(nsite, nsamp)) # array for mean distance between communities in composition of species 
for(k in 1:nsamp){
  print(k)
  MCD[,k]<-apply((as.matrix(beta.div.comp(N[,,k], coef = "BS", quant = T)$repl)), MARGIN = 1, mean)
} 

pm <- apply(MCD, 1, mean, na.rm = TRUE) # Post. mean of Jsite wrt. site 1
psd <- apply(MCD, 1, sd, na.rm = TRUE) # Post. sd of Jsite wrt. site 1
cri <- apply(MCD, 1, function(x) quantile(x, prob = c(0.025, 0.975), na.rm =TRUE)) # CRI
MCD.pm<-as.data.frame(cbind('MCD.pm' = pm, 'MCD.sd' = psd, 'MCD-2.5%' = cri[1,], 'MCD-97.5%' = cri[2,])) 

# Nestedness component of Baselga's Multipart Beta Diversity
MND <- array(NA, dim = c(nsite, nsamp)) # array for mean distance between communities in composition of species 
for(k in 1:nsamp){
  print(k)
  MND[,k]<-apply((as.matrix(beta.div.comp(N[,,k], coef = "BS", quant = T)$rich)), MARGIN = 1, mean)
} 

pm <- apply(MND, 1, mean, na.rm = TRUE) # Post. mean of Jsite wrt. site 1
psd <- apply(MND, 1, sd, na.rm = TRUE) # Post. sd of Jsite wrt. site 1
cri <- apply(MND, 1, function(x) quantile(x, prob = c(0.025, 0.975), na.rm =TRUE)) # CRI
MND.pm<-as.data.frame(cbind('MND.pm' = pm, 'MND.sd' = psd, 'MND-2.5%' = cri[1,], 'MND-97.5%' = cri[2,])) 
com.pm<-as.data.frame(cbind(SR.pm["SR.pm"], MCD.pm["MCD.pm"], MND.pm["MND.pm"]))
com.sd<-as.data.frame(cbind(SR.pm["SR.sd"], MCD.pm["MCD.sd"], MND.pm["MND.sd"]))

# Community Model 
# Bundle and summarize data set
pred.ldi<-(seq(0.01, 0.99,,500)-mean(orig.ldi))/sd(orig.ldi) # predicted ldi
str(win.data <- list(pm=as.matrix(com.pm), psd=as.matrix(com.sd) ,n = dim(com.pm)[1], 
                     var1=ldi, var2=me.res, var3=hl.res,  var4=map, 
                    nv=dim(com.pm)[2], pred.ldi=pred.ldi, npred=length(pred.ldi)))
# initial values for chain start 
inits <- function() list(beta0=runif(dim(com.pm)[2], 0,0.01)) 
params <- c("beta0",  "beta", "Npred")
#trmod2com<- jags(win.data, inits, params, "trmod2com.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb) 

# 4.4 Tree Results -------
# 4.4.1 Tree Species Richness  -----
obsldi<-as.data.frame(cbind(x=orig.ldi, y=SR.pm$SR.pm, ymin = SR.pm$`SR-2.5%`, ymax =  SR.pm$`SR-97.5%`))
ldisr<-as.data.frame(cbind(x=seq(0, 0.99,,500), y= trmod2com$mean$Npred[1,], ymin = trmod2com$q2.5$Npred[1,], ymax = trmod2com$q97.5$Npred[1,]))
base<-max(trmod2com$mean$Npred[1,])# Thresholds 
band1<-as.data.frame(cbind(x=seq(0, 0.99,,500),y=base,  ymin =0.60*base, ymax =0.79*base))
band1$thtype<-"Intermediate"
band2<-as.data.frame(cbind(x=seq(0, 0.99,,500),y=base,  ymin =0.40*base, ymax =0.59*base))
band2$thtype<-"High"
th<-rbind(band2, band1); th$thtype<-factor(th$thtype, levels=c('Intermediate', "High"))

pred.sr<-ldisr
pred.sr$y1 <- ((pred.sr$y/base)-1)*100 ; pred.sr$ymin1 <- ((pred.sr$ymin/base)-1)*100 ; pred.sr$ymax1 <- ((pred.sr$ymax/base)-1)*100
pred.sr$restype<-"Richness"

fig5a_trsr_ldi<-ggplot(data = ldisr, aes(y = y, x = x, ymin = ymin, ymax = ymax)) +
  geom_line(colour="darkblue", size=0.75, alpha=0.4) + geom_ribbon(alpha = 0.175, fill = "blue") + 
  geom_point(data = obsldi, aes(x = x, y = y), size=4, alpha=0.4, shape=19, colour="grey70")+ 
  geom_linerange(data=obsldi, aes(ymin=ymin, ymax=ymax), size=.2,  colour="grey70", alpha=0.2) + 
  scale_x_continuous(limits = c(0, 1), breaks=seq(0,1,0.25), "Land Division Index") + 
theme(legend.position="bottom",legend.title=element_text(colour="gray25", size = 11, family = "Calibri"), 
        legend.text=element_text(colour="gray30", size = 11, family = "Calibri"),
        axis.text=element_text(colour="grey50", size = rel(0.99), family = "Calibri"), 
        axis.title.x = element_text(colour = "grey40", family = "Calibri", size=12), 
        axis.title.y = element_text(colour = "grey40", family = "Calibri", size=12),panel.background = element_blank()) + 
  theme(axis.line = element_line(color = 'grey75', size=0.2)) + 
  scale_y_continuous(limits = c(10, 75),"Estimated species richness of trees", 
              sec.axis = sec_axis(~ ((./base)-1)*100 , name = "% reduction in species richness", breaks = c(-75, -50, -25, 0))) +   
  geom_ribbon(data=th, aes(colour = factor(thtype), fill = factor(thtype)), linetype = 3, size=0.1, alpha= 0.125) + 
  scale_fill_manual(values = c('orange','red')) + guides(colour=FALSE, fill = guide_legend(override.aes = list(size=4))) +  
  labs(fill = "Levels of species loss (Hooper et al. 2012)") + geom_smooth(data = obsldi, method = "lm", formula = y ~ splines::bs(x, 3), 
                                                                           se = FALSE, colour="grey50", linetype=2, alpha=0.75, size=0.35)
fig5a_trsr_ldi
ggsave(plot=fig5a_trsr_ldi, "fig5a_trsr_ldi.png",width = 5, height = 6) 

# 4.4.2 Tree Species Abundance  -----
tree.abund<-as.data.frame(cbind(trmod1C$mean$mlambda))
colnames(tree.abund)<-"abund"
tree.abund$c1<-trmod1C$q2.5$mlambda
tree.abund$c2<-trmod1C$q97.5$mlambda
tree.abund<-tree.abund[1:nspec,] # select only observed species 
tree.abund$resname<-unique(tr$sciname)
tree.abund$res.sd<-apply(Yc[,1,], 2, sd)

# Plot of Tree Abundance
tr.abund<-ggplot(tree.abund, aes(x = reorder(resname, abund) , y = abund, ymin=c1, ymax = c2)) + 
  geom_point(aes(size =abund), alpha = 0.5, shape = 19, colour = "black") + scale_size_continuous(range = c(2, 5)) + 
  geom_linerange(size=0.25, alpha = 0.75, color="grey25")+ coord_flip()+ 
  ylab("Abundance") + xlab("Tree species") + 
  theme(legend.position="none", legend.title=element_blank(), 
        legend.text=element_text(colour="gray30", size = 10, family = "Calibri"), 
        axis.text=element_text(colour="grey50", size = rel(0.65), family = "Calibri", face="italic"), 
        axis.title.x = element_text(colour = "grey30", family = "Calibri", size=13), 
        axis.ticks.y=element_blank(),
        axis.title.y = element_text(colour = "grey30", family = "Calibri", size=13),
        strip.text.x = element_text(colour = "grey40", family = "Calibri"),  panel.background = element_blank())
tr.abund
ggsave(plot=tr.abund, "tr_sps_abund.png", width = 6.1, height = 9)  

# 4.4.3 Extract predicted species richness, turnover and nestedness of tree species -----
# Richness
pred.sr<-as.data.frame(cbind(x=seq(0, 0.99,,500), y= trmod2com$mean$Npred[1,], ymin = trmod2com$q2.5$Npred[1,], ymax = trmod2com$q97.5$Npred[1,]))
base<-max(trmod2com$mean$Npred[1,])
pred.sr$y1 <- ((pred.sr$y/base)-1)*100
pred.sr$ymin1 <- ((pred.sr$ymin/base)-1)*100
pred.sr$ymax1 <- ((pred.sr$ymax/base)-1)*100
pred.sr$restype<-"Richness"
head(pred.sr)


# Turnover
pred.turn<-as.data.frame(cbind(x=seq(0, 0.99,,500), y= trmod2com$mean$Npred[2,], ymin = trmod2com$q2.5$Npred[2,], ymax = trmod2com$q97.5$Npred[2,]))
base<-min(trmod2com$mean$Npred[2,])
pred.turn$y1 <- ((pred.turn$y/base)-1)*100
pred.turn$ymin1 <- ((pred.turn$ymin/base)-1)*100
pred.turn$ymax1 <- ((pred.turn$ymax/base)-1)*100
pred.turn$restype<-"Turnover"
head(pred.turn)

# Nestedness
pred.nest<-as.data.frame(cbind(x=seq(0, 0.99,,500), y= trmod2com$mean$Npred[3,], 
                               ymin = trmod2com$q2.5$Npred[3,], ymax = trmod2com$q97.5$Npred[3,]))
base<-max(trmod2com$mean$Npred[3,])
pred.nest$y1 <- ((pred.nest$y/base)-1)*100
pred.nest$ymin1 <- ((pred.nest$ymin/base)-1)*100
pred.nest$ymax1 <- ((pred.nest$ymax/base)-1)*100
pred.nest$restype<-"Nestedness"
head(pred.nest)

pred.tree<-as.data.frame(rbind(pred.sr, pred.turn, pred.nest))
pred.tree$taxa<-"Trees"
head(pred.tree)

# combine all predicted diversity 
pred.div<-as.data.frame(rbind(pred.mam, pred.tree))
mean(pred.mam$y1)


# Observed raw values 

obs.sr<-as.data.frame(cbind(x=orig.ldi, y=SR.pm$SR.pm, ymin = SR.pm$`SR-2.5%`, ymax =  SR.pm$`SR-97.5%`)); head(obssr)
obs.sr$restype<-"Richness"

obs.turn<-as.data.frame(cbind(x=orig.ldi, y=MCD.pm$MCD.pm, ymin = MCD.pm$`MCD-2.5%`, ymax =  MCD.pm$`MCD-97.5%`)); head(obs.turn) 
obs.turn$restype<-"Turnover"

obs.nest<-as.data.frame(cbind(x=orig.ldi, y=MND.pm$MND.pm, ymin = MND.pm$`MND-2.5%`, ymax =  MND.pm$`MND-97.5%`)); head(obs.nest) 
obs.nest$restype<-"Nestedness"

obs.tree<-as.data.frame(rbind(obs.sr, obs.turn, obs.nest))
obs.tree$taxa<-"Trees"
head(obs.tree)

# combine all obsicted diversity 
obs.div<-as.data.frame(rbind(obs.mam, obs.tree)); head(obs.div)

# make a plot 
sr<-pred.div[pred.div$restype=="Richness",]
turn<-pred.div[pred.div$restype=="Turnover",]
nest<-pred.div[pred.div$restype=="Nestedness",]




ggplot(data = sr, aes(y = y1, x = x, ymin = ymin, ymax = ymax)) +
  geom_line(colour="darkblue", size=0.75, alpha=0.4) + geom_ribbon(alpha = 0.175, fill = "blue") + 
  geom_point(data = obsldi, aes(x = x, y = y), size=4, alpha=0.4, shape=19, colour="grey70")+ 
  geom_linerange(data=obsldi, aes(ymin=ymin, ymax=ymax), size=.2,  colour="grey70", alpha=0.2) + 
  scale_x_continuous(limits = c(0, 1), breaks=seq(0,1,0.25), "Land Division Index") + 
  theme(legend.position="bottom",legend.title=element_text(colour="gray25", size = 11, family = "Calibri"), 
        legend.text=element_text(colour="gray30", size = 11, family = "Calibri"),
        axis.text=element_text(colour="grey50", size = rel(0.99), family = "Calibri"), 
        axis.title.x = element_text(colour = "grey40", family = "Calibri", size=12), 
        axis.title.y = element_text(colour = "grey40", family = "Calibri", size=12),panel.background = element_blank()) + 
  theme(axis.line = element_line(color = 'grey75', size=0.2)) + 
  scale_y_continuous(limits = c(10, 75),"Estimated species richness of trees", 
                     sec.axis = sec_axis(~ ((./base)-1)*100 , name = "% reduction in species richness", breaks = c(-75, -50, -25, 0))) +   
  geom_ribbon(data=th, aes(colour = factor(thtype), fill = factor(thtype)), linetype = 3, size=0.1, alpha= 0.125) + 
  scale_fill_manual(values = c('orange','red')) + guides(colour=FALSE, fill = guide_legend(override.aes = list(size=4))) +  
  labs(fill = "Levels of species loss (Hooper et al. 2012)") + geom_smooth(data = obsldi, method = "lm", formula = y ~ splines::bs(x, 3), 
                                                                           se = FALSE, colour="grey50", linetype=2, alpha=0.75, size=0.35)






p1<-ggplot() +geom_line(data = sr, aes(x = x, y = y1, color =  taxa), size = 0.75)+
  geom_point(data = obs.div, aes(x = x, y = y), size=4, alpha=0.4, shape=19, colour="grey70") + 
  theme(legend.position="bottom",legend.title=element_blank(), 
        legend.text=element_text(colour="gray50", size = 11, family = "Calibri"),
        axis.text=element_text(colour="grey55", size = rel(0.8), family = "Calibri"), 
        axis.title.x = element_text(colour = "grey40", family = "Calibri", size=12), 
        axis.title.y = element_text(colour = "grey40", family = "Calibri", size=12), 
        strip.text.x = element_text(colour = "gray35", family = "Calibri", size=12),panel.background = element_blank())+
  theme(axis.line = element_line(color = 'grey75', size=0.2)) + 
  scale_y_continuous(limits = c(-50, 10), breaks = c(-40, -20,  0, 20, 40),  "Percentage change") + 
  scale_x_continuous(limits = c(0.1, 0.99), breaks=seq(0,1,0.2), "Land Division Index") + 
  scale_colour_manual(breaks = c('Trees', 'Mammals'), values=c( "coral3", "springgreen3")) + 
  guides(fill=FALSE, colour = guide_legend(override.aes = list(size=2))) ; p1
ggsave(plot=p1, "sr_mam_tree.png", width = 4, height = 3.5) 


p2<-ggplot() +geom_line(data = turn, aes(x = x, y = y1, color =  taxa), size = 0.75)+
  theme(legend.position="bottom",legend.title=element_blank(), 
        legend.text=element_text(colour="gray50", size = 11, family = "Calibri"),
        axis.text=element_text(colour="grey55", size = rel(0.8), family = "Calibri"), 
        axis.title.x = element_text(colour = "grey40", family = "Calibri", size=12), 
        axis.title.y = element_text(colour = "grey40", family = "Calibri", size=12), 
        strip.text.x = element_text(colour = "gray35", family = "Calibri", size=12),panel.background = element_blank())+
  theme(axis.line = element_line(color = 'grey75', size=0.2)) + 
  scale_y_continuous(limits = c(-30, 20), breaks = c(-40, -20,  0, 20, 40),  "Percentage change") + 
  scale_x_continuous(limits = c(0.1, 0.99), breaks=seq(0,0.8,0.2), "Land Division Index") + 
  scale_colour_manual(breaks = c('Trees', 'Mammals'), values=c( "coral3", "springgreen3")) + 
  guides(fill=FALSE, colour = guide_legend(override.aes = list(size=2))) ; p2
ggsave(plot=p2, "turn_mam_tree.png", width = 4, height = 3.5) 


p3<-ggplot() +geom_line(data = nest, aes(x = x, y = y1, color =  taxa), size = 0.75)+
  theme(legend.position="bottom",legend.title=element_blank(), 
        legend.text=element_text(colour="gray50", size = 11, family = "Calibri"),
        axis.text=element_text(colour="grey55", size = rel(0.8), family = "Calibri"), 
        axis.title.x = element_text(colour = "grey40", family = "Calibri", size=12), 
        axis.title.y = element_text(colour = "grey40", family = "Calibri", size=12), 
        strip.text.x = element_text(colour = "gray35", family = "Calibri", size=12),panel.background = element_blank())+
  theme(axis.line = element_line(color = 'grey75', size=0.2)) + 
  scale_y_continuous(limits = c(-60, 40), breaks = c(-40, -20,  0, 20, 40),  "Percentage change") + 
  scale_x_continuous(limits = c(0.1, 0.99), breaks=seq(0,0.8,0.2), "Land Division Index") + 
  scale_colour_manual(breaks = c('Trees', 'Mammals'), values=c( "coral3", "springgreen3")) + 
  guides(fill=FALSE, colour = guide_legend(override.aes = list(size=2))) ; p3
ggsave(plot=p3, "nest_mam_tree.png", width = 4, height = 3.5) 


pred.div$restype<-factor(pred.div$restype, levels=c("Richness", "Turnover", "Nestedness"))

p4<-ggplot() +
  geom_line(data = pred.div, aes(x = x, y = y1, color = taxa), size = 1)+
  geom_ribbon(data =pred.div, aes(x = x, ymax = ymax1, ymin=ymin1, fill = taxa), alpha=0.12)+ 
  theme(legend.position="bottom",legend.title=element_blank(), 
        legend.text=element_text(colour="gray50", size = 11, family = "Calibri"),
        axis.text=element_text(colour="grey55", size = rel(0.8), family = "Calibri"), 
        axis.title.x = element_text(colour = "grey40", family = "Calibri", size=12), 
        axis.title.y = element_text(colour = "grey40", family = "Calibri", size=12), 
        strip.text.x = element_text(colour = "gray35", family = "Calibri", size=12),panel.background = element_blank())+
  theme(axis.line = element_line(color = 'grey75', size=0.2)) + 
  scale_y_continuous( limits = c(-100, 160), breaks = c(-40, -20, 0, 20, 40),  "Percentage change") + 
  scale_x_continuous(limits = c(0.15, 0.95),breaks=seq(0,0.8,0.2), "Land Division Index") + 
  scale_colour_manual(breaks = c('Trees', 'Mammals'), 
                      values=c("palegreen3", "dodgerblue3")) + 
  scale_fill_manual(breaks = c('Trees', 'Mammals'), 
                    values=c("palegreen3", "dodgerblue3")) + 
  facet_grid( ~ restype, scales = "free", shrink=T, drop=TRUE) + 
  guides(fill=FALSE, colour = guide_legend(override.aes = list(size=2)))  ; p4

ggsave(plot=p4, "tree_mammal_diversity_panel.png", width = 7.5, height = 4.25) 


# 4.4.4 Results - Extract Coefficent table for Tree models  -----
# Model 1 - Abundance
slope.w<-as.data.frame(cbind(trmod1C$mean$beta))
slope.w$intercept<-trmod1C$mean$beta0 # Intercept 
slope.w$intercept.sd<-trmod1C$sd$beta0 # sd of intercept 

pred.var<-colnames(slope.w)[1:6] # number predictor veriables
pred.names<-c("LDI", "LDI2","LDI3", "ME.resid", "HL.resid", "MAP")

slope<-reshape(slope.w, varying = pred.var, v.names = "slope", timevar = "predictor", times = pred.names, direction = "long")

prob.pos<-apply(ifelse(trmod1C$sims.list$beta>0, 1, 0), c(2,3), function (x) sum(x)/3000) # probablity of +ive
slope$prob.pos<-reshape(as.data.frame(prob.pos), varying = pred.var, v.names = "prob.pos", direction = "long")[,2]
prob.neg<-apply(ifelse(trmod1C$sims.list$beta<0, 1, 0),c(2,3), function (x) sum(x)/3000) # probablity of -ive
slope$prob.neg<-reshape(as.data.frame(prob.neg),varying = pred.var,v.names = "prob.neg", direction = "long")[,2]

slope$c1<-reshape(as.data.frame(cbind(trmod1C$q2.5$beta)),varying = pred.var, v.names = "c1",  direction = "long")[,2]
slope$c2<-reshape(as.data.frame(cbind(trmod1C$q97.5$beta)),varying = pred.var, v.names = "c2", direction = "long")[,2]
slope$sd<-reshape(as.data.frame(cbind(trmod1C$sd$beta)),varying = pred.var, v.names = "sd", direction = "long")[,2]
mod1.res<-slope

mod1.res$model<-rep("mod1", length(mod1.res$slope))
mod1.res<-mod1.res[mod1.res$id<=nspec,] # select only observed species 
mod1.res$res.var<-rep("C", length(mod1.res$slope)) 
mod1.res$res.sd<-1  
mod1.res$resname<-unique(tr$sciname)
head(mod1.res); dim(mod1.res) # should be nspec* 6 *14

# Community mean coefficients from model 1
slope.w<-as.data.frame(cbind(t(trmod1C$mean$mu.beta)))
slope.w$intercept<-trmod1C$mean$mu.beta0 
slope.w$intercept.sd<-trmod1C$sd$mu.beta0
slope<-reshape(slope.w, varying = pred.var, v.names = "slope", timevar = "predictor",times = pred.names, direction = "long")
slope$prob.pos<-apply(ifelse(trmod1C$sims.list$mu.beta>0, 1, 0),2, function (x) sum(x)/3000) # probablity of +ive
slope$prob.neg<-apply(ifelse(trmod1C$sims.list$mu.beta<0, 1, 0),2, function (x) sum(x)/3000) # probablity of -ive

slope$c1<-t(trmod1C$q2.5$mu.beta)[1,]
slope$c2<-t(trmod1C$q97.5$mu.beta)[1,]
slope$sd<-t(trmod1C$sd$mu.beta)[1,]

com.mean.res<-slope
com.mean.res$model<-rep("com.mean", length(com.mean.res$slope))
com.mean.res$res.var<-rep("C", length(com.mean.res$slope)) 
com.mean.res$res.sd<-1  
com.mean.res$resname<-NA; com.mean.res[is.na(com.mean.res)]<-0
head(com.mean.res); dim(com.mean.res) # 

# model 2 - diversity model coefficients
slope.w<-as.data.frame(cbind(trmod2com$mean$beta))
slope.w$intercept<-trmod2com$mean$beta0 # 1st level  of Alpha1  ie. boom/ intensity class preboom is the intercept 
slope.w$intercept.sd<-trmod2com$sd$beta0
slope<-reshape(slope.w, varying = pred.var, v.names = "slope",timevar = "predictor", times = pred.names, direction = "long")

prob.pos<-apply(ifelse(trmod2com$sims.list$beta>0, 1, 0),c(2,3), function (x) sum(x)/3000) # probablity of +ive
slope$prob.pos<-reshape(as.data.frame(prob.pos), varying = pred.var, v.names = "prob.pos", direction = "long")[,2]

# probbality of negative effect 
prob.neg<-apply(ifelse(trmod2com$sims.list$beta<0, 1, 0), c(2,3), function (x) sum(x)/3000) # probablity of -ive
slope$prob.neg<-reshape(as.data.frame(prob.neg), varying = pred.var,v.names = "prob.neg", direction = "long")[,2]
slope$c1<-reshape(as.data.frame(cbind(trmod2com$q2.5$beta)),varying = pred.var, v.names = "c1",  direction = "long")[,2]
slope$c2<-reshape(as.data.frame(cbind(trmod2com$q97.5$beta)),varying = pred.var, v.names = "c2", direction = "long")[,2]
slope$sd<-reshape(as.data.frame(cbind(trmod2com$sd$beta)), varying = pred.var, v.names = "sd", direction = "long")[,2]
mod2.res<-slope
mod2.res$model<-rep("mod2", length(mod2.res$slope))
mod2.res$res.var<-rep("Com", length(mod2.res$slope)) 
mod2.res$res.sd<-apply(com.pm, 2,  sd)
mod2.res$resname<-c("SR", "MCD", "MND")
head(mod2.res); dim(mod2.res)
trmod.res<-rbind(mod1.res, mod2.res, com.mean.res) # tree model results 
trmod.res$species_group<-"Trees"

mod.res<-as.data.frame(rbind(mmod.res, trmod.res)) # result of tree and mammal model together 

# Add idicators signficant values
ngpos<-rep(0, length(mod.res$slope)) # negative or positive? 
ngpos[mod.res$slope>0]<-1 # where more than 0, +1
ngpos[mod.res$slope<0]<--1 # where less than 0, -1
overlap<-rep(0, length(mod.res$slope)) # overlap indicates signficance, when CI overlap, insignificant
overlap[mod.res$c1*mod.res$c2>0]<-1 # signficant 
sig<-ngpos*overlap # signficant neg and pos
mod.res$ngpos<-ngpos
mod.res$sig<-sig

postprob<-rep(0, length(mod.res$slope))
postprob[mod.res$prob.pos>0.94]<-1 # signficant 
postprob[mod.res$prob.neg>0.94]<- -1
mod.res$postprob<-postprob; save.image()

# Standardize and Scale coefficients 
raw.coef<-mod.res[c("slope", "sd", "c1", "c2")]
scaled.coef<-apply(raw.coef, 2, function (x) (x/mod.res$res.sd))
colnames(scaled.coef) <-c("std.slope", "std.sd",  "std.c1", "std.c2")
mod.res<-as.data.frame(cbind(mod.res, scaled.coef))
rownames(mod.res)<-NULL
head(mod.res) 
write.csv (mod.res, "modres.csv")

# Coefficient plots 
m=mod.res

# order and rename predictors
predictors=c('LDI',"LDI2", "LDI3",'ME.resid','HL.resid', "MAP")
m$predictor<-factor(m$predictor, levels=predictors)

# Coefficent plot of mammals 
mam<-m[m$model=="mod1"& m$species_group=="Mammals",]
mam$occu=mam.occu$occu1
mamof1 <- ggplot(mam, aes(x = reorder(resname, (resname))  , y =slope, ymin=c1, ymax =c2, 
                          colour=factor(postprob), shape=factor(postprob)))+ 
  scale_colour_manual(values = c('coral2','grey', "deepskyblue1"), 
                      labels=c("Negative", "Non-Signficant", "Positive")) +
  scale_shape_manual(values = c(19,1, 19), labels=c("Negative", "Non-Signficant", "Positive")) + 
  geom_point(aes(size =occu), alpha = 0.75) + scale_size_continuous(range = c(2, 6)) + 
  geom_linerange(size=0.5, alpha = 0.5 )+ coord_flip()+ 
  geom_hline(yintercept = 0,linetype = "dotted", size=0.5)+ 
  ylab("Parameter estimates") + xlab("Tree species") + 
  theme(legend.position="bottom", legend.title=element_blank(), 
        legend.text=element_text(colour="gray30", size = 10, family = "Calibri"), 
        axis.text=element_text(colour="gray50", size = rel(0.7), family = "Calibri", face="italic"), 
        axis.title.x = element_text(colour = "grey40", family = "Calibri", size=11), 
        axis.title.y = element_text(colour = "grey40", family = "Calibri", size=11),
        axis.ticks.y=element_blank(),
        strip.text.x = element_text(colour = "gray30", family = "Calibri", size=11),  panel.background = element_blank())+
  facet_grid( ~predictor, scales = "free", shrink=T, drop=TRUE) +
  guides(size=FALSE, colour = guide_legend(override.aes = list(size=4))); mamof1
ggsave(plot=mamof1, "mam_sps_occu.png",width = 7.4, height = 5) 

# Coefficient plot of trees
tr<-m[m$model=="mod1"& m$species_group=="Trees",]
tr$abund=tree.abund$abund
tcf1 <- ggplot(tr, aes(x = reorder(resname, (resname))  , y =slope, ymin=c1, ymax =c2, colour=factor(postprob), shape=factor(postprob)))+ scale_colour_manual(values = c('coral2','grey', "deepskyblue1"), labels=c("Negative", "Non-Signficant", "Positive")) +scale_shape_manual(values = c(19,1, 19), labels=c("Negative", "Non-Signficant", "Positive")) + 
  geom_point(aes(size =abund), alpha = 0.75) + scale_size_continuous(range = c(2, 6)) + 
  geom_linerange(size=0.5, alpha = 0.5 )+ coord_flip()+ 
  geom_hline(yintercept = 0,linetype = "dotted", size=0.5)+ 
  ylab("Parameter estimates") + xlab("Tree species") + 
  theme(legend.position="bottom", legend.title=element_blank(), 
        legend.text=element_text(colour="gray30", size = 10, family = "Calibri"), 
        axis.text=element_text(colour="gray50", size = rel(0.7), family = "Calibri", face="italic"), 
        axis.title.x = element_text(colour = "grey40", family = "Calibri", size=11), 
        axis.title.y = element_text(colour = "grey40", family = "Calibri", size=11),
        axis.ticks.y=element_blank(),
        strip.text.x = element_text(colour = "gray30", family = "Calibri", size=11),  panel.background = element_blank())+
  facet_grid( ~predictor, scales = "free", shrink=T, drop=TRUE) +
  guides(size=FALSE, colour = guide_legend(override.aes = list(size=4))); tcf1
ggsave(plot=tcf1, "trC.png",width = 7.4, height = 9.1) 

#Community level coefficient plots
# Community Occurence
com<-m[(m$model=="com.mean"),]

library(scales)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

com.occu <- ggplot(com, aes(x =species_group, y = std.slope, ymin=std.c1, ymax = std.c2, colour=factor(postprob), shape=factor(postprob)))+ 
  scale_colour_manual(values = c('coral2','grey40', "deepskyblue1"), labels=c("Negative", "Non-Signficant", "Positive")) +
  scale_shape_manual(values = c(19,1, 19), labels=c("Negative", "Non-Signficant", "Positive")) + 
  geom_point(size=4, alpha = 0.5) +    scale_y_continuous(breaks=number_ticks(2)) + 
  geom_linerange(size=0.75, alpha = 0.5 )+ coord_flip()+ 
  geom_hline(yintercept = 0,linetype = "dotted", size=0.5)+ 
  ylab("Parameter estimates") + xlab("") + 
  theme(legend.position="bottom", legend.title=element_blank(), 
        legend.text=element_text(colour="gray30", size = 11, family = "Calibri"), 
        axis.text=element_text(colour="gray50", size = 11, family = "Calibri"), 
        axis.title.x = element_text(colour = "grey40", family = "Calibri", size=12), 
        axis.title.y = element_text(colour = "grey40", family = "Calibri", size=12),
        axis.ticks.y=element_blank(),
        strip.text.x = element_text(colour = "gray30", family = "Calibri", size=12),  panel.background = element_blank()) + facet_grid( ~predictor,shrink=T,scales = "free", drop=TRUE) +   guides(size=FALSE, colour = guide_legend(override.aes = list(size=4)))
com.occu 
ggsave(plot=com.occu, "com.occu .png",width = 7.4, height = 3) 

# Community Diversity
div<-m[(m$model=="mod2"),]

# order and rename Diveristy variables
diversity.var=c("SR", "MCD", "MND")
diversity.var.rename<-c("Richness", "Turnover", "Nestedness")

div$resname<-factor(div$resname, levels=diversity.var)
levels(div$resname)<-diversity.var.rename

com.div <- ggplot(div, aes(x =resname, y = std.slope, ymin=std.c1, ymax = std.c2, colour=factor(postprob), shape=factor(postprob)))+ 
  scale_colour_manual(values = c('coral2','grey40', "deepskyblue1"), labels=c("Negative", "Non-Signficant", "Positive")) +
  scale_shape_manual(values = c(19,1, 19), labels=c("Negative", "Non-Signficant", "Positive")) + 
  geom_point(size=4, alpha = 0.5) +    scale_y_continuous(breaks=number_ticks(3)) + 
  geom_linerange(size=0.75, alpha = 0.5 )+ coord_flip()+ 
  geom_hline(yintercept = 0,linetype = "dotted", size=0.5)+ 
  ylab("Parameter estimates") + xlab("") + 
  theme(legend.position="bottom", legend.title=element_blank(), 
        legend.text=element_text(colour="gray30", size = 11, family = "Calibri"), 
        axis.text=element_text(colour="gray50", size = 11, family = "Calibri"), 
        axis.title.x = element_text(colour = "grey40", family = "Calibri", size=11), 
        axis.title.y = element_text(colour = "grey40", family = "Calibri", size=11),
        axis.ticks.y=element_blank(),
        strip.text.x = element_text(colour = "gray30", family = "Calibri", size=11),  panel.background = element_blank()) + facet_grid( species_group~predictor,shrink=T,scales = "free", drop=TRUE) +   guides(size=FALSE, colour = guide_legend(override.aes = list(size=4)))
com.div
ggsave(plot=com.div, "com.div.png",width = 7.5, height = 3) 

# Correlation between predictors
mam<-read.csv("mam.csv", header=TRUE); head(m) # mammal dataset 
tree<-read.csv("tr.csv", header=TRUE); head(tr) # tree abundance dataset
preds1<-as.data.frame(cbind(append(tree$ldi, mam$ldi), append(tree$loss_intensity, mam$loss_intensity), 
                            append(tree$miombo_extent, mam$miombo_extent), append(tree$MAP, mam$MAP)))
colnames(preds1)<-c("LDI", "HL", "ME", "MAP")
library(PerformanceAnalytics)
chart.Correlation (preds1)
save.image()


# mean loss of diveristy at different levels of fragmentation 

#Richness 
# at90% fragmentation 
sr90<-sr[sr$x>0.95,]
mean(sr90$y1); (mean(sr90$ymax1)-mean(sr90$ymin1)) / 3.92

# across all intermediate landscapes 
sr25<-sr[sr$x>0.25,]
mean(sr25$y1); (mean(sr$ymax1)-mean(sr25$ymin1)) / 3.92

nest[nest$y1==0,]


