# File to provide summary, PMF source apportionment
# results of  QCII speciation monitor in NYC
# 3/13/13
# ========================================================



#load libraries
library(xtable)
library(ggplot2)
library(scales)

#get data
dat <- read.csv("http://biostat.jhsph.edu/~jkrall/nycdat.csv", stringsAsFactors = F)
mdlsboth <- read.csv("http://biostat.jhsph.edu/~jkrall/mdls.csv", stringsAsFactors = F)
mdlmax <- as.numeric(mdlsboth[1, ])
names(mdlmax) <- colnames(mdlsboth)
mdlmin <- as.numeric(mdlsboth[2, ])
names(mdlmin) <- colnames(mdlsboth)
dat[, 1] <- as.Date(dat[, 1])


dat.dir <- "C:/Users/Jenna Krall/Dropbox/MDL_sourceapp/nyc_analysis_mdl"

##
# load data
load("dataNYC.RData")
load("apca_nyc.RData")


#mdl matrix
mdls <- matrix(rep(mdlmax, each = nrow(dat)), 
	ncol = length(mdlmax), byrow = F)
colnames(mdls) <- colnames(mdlsboth)


#get uncertainty as in Ito et al. 2004
unc <- 0.05 * dataNYC[[1]][,-1] + mdls


#fix snratio
dataNYC.snrat <- snratALL(dataNYC[[1]][, -1], unc, mdls)

#add back in PM2.5
datmcmc <- list()
for(i in 1 : dim(dataNYC[[4]])[3]) {
	datmcmc[[i]] <- data.frame(dataNYC[[1]][, 2], dataNYC[[4]][,,i])
}


dataNYCPMF <- list(dataNYC[[1]][, -1], dataNYC[[2]][, -1], 
	dataNYC.snrat[[1]], datmcmc)
unc.snrat <- dataNYC.snrat[[2]]
save(dataNYCPMF, unc.snrat, file = file.path(dat.dir, "dataNYC_pmf.RData"))







##########
# apply PMF


gwds <- "C:/Users/Jenna Krall/Dropbox/MDL_sourceapp/PMF/tryme"
setwd(gwds)

fp1 <- file.path(dat.dir, "pmf_nyc.RData")

# pmf.nyc <- pmf.all(dataNYCPMF, unc, unc.snrat, nf = 8, 
	# fp = fp1 )
# save(pmf.nyc, file =fp1)






########
# compare with APCA


load("pmf_nyc.RData")
cn <- colnames(dataNYC[[1]])[-c(1, 2)]

# #first match no missing between APCA and PMF
round(cor(pmf.nyc[[1]][[1]], source.conc), 3)
round(cor(pmf.nyc[[1]][[2]][-1, ], source.prof), 3)

dista.conc <- dist(t(cbind(pmf.nyc[[1]][[1]], source.conc)))
dista <- dist(t(cbind(pmf.nyc[[1]][[2]][-1, ], source.prof)))
ncol1 <- ncol(pmf.nyc[[1]][[2]])
ncol2 <- ncol(source.conc)
dista <- as.matrix(dista)[1 : ncol1, (ncol1 + 1) : (ncol1 + ncol2)]
dista.conc <- as.matrix(dista.conc)[1 : ncol1, (ncol1 + 1) : (ncol1 + ncol2)]

#1 is soil
#2 is sec sulfate
#3 is traffic
#4  is res oil



