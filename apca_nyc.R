# File to provide summary, APCA SA results of  QCIImonitor in NYC
# 8/20/13
# ========================================================




library(ncdf4)


#Obtain adjusted data
################################
################################
################################
################################
################################
################################


#summary of bdl values
bqcii <- bdlfun(dat, mdlmax)


# ####
# # Get adjusted data for NYC
dataNYC <- createDATs(dat, mdlmax)
set.seed(1597)
#create matrix of mdls
mdls <- matrix(rep(mdlmax, each = nrow(dat)), 
	ncol = length(mdlmax), byrow = F)
#run MCMC
datmcmc <- runGIBBS(dat[, -1], mdls, niter = 50000, 
      burnin = 25000, ndraws = 100, seed = 2361)
dataNYC[[4]] <- datmcmc
names(dataNYC)[4] <- "gibbs" 



# save(dataNYC, mdlmat, file = file.path(newcode.dir, "dataNYC_20aug13.RData"))
# #################################
# #################################
# #################################
# #################################
# #################################
# #################################


















datNYCold <- dataNYC
dataNYC <- datNYCold
################################
################################
################################
################################
################################
################################
####
# # APPLY APCA to each
set.seed(99466)
#only use first 10 for gibbs
# dataNYC[[4]] <- dataNYC[[4]][, , c(1 : 10)]
apca.nyc <- apca.all(dataNYC, tot = T)
# save(apca.nyc, file = "apca_nyc_16aug13.RData")



# #########
# #first ID sources in data
# #for ITO:
# #ec, OC, K for traffic
# #ammonium, S, for sec sulf
# #"aluminum", "calcium", "iron", "silicon","titanium" for soil
# #"nickel", "vanadium", "nitrate",  "chlorine", "lead", "zinc" res oil
# ####
# #A=soil, B=secsulf, C=traffic, D= resoil

# ######
# #match results to each other
# #fix reported
source.conc <- apca.nyc[[1]][[1]][, c(3, 2, 4, 1)]
source.prof <- apca.nyc[[1]][[4]][, c(3, 2, 4, 1)]
source.conc.all <- match.sources(apca.nyc, source.conc, "scores", holds = c(6, 7))
source.prof.all <- match.sources(apca.nyc, source.prof, "prof", holds = c(6, 7))


# ## DOUBLE CHECK
matchmat <- matrix(nrow = 13, ncol = 4)
pr.load(dataNYC[[1]][, -1], 8)
matchmat[1, ] <- c(3, 2, 4, 1)
pr.load(dataNYC[[2]][, -1], 8)
matchmat[2, ] <- c(2, 3, 4, 1)
pr.load(dataNYC[[3]][, -1], 8)
matchmat[3, ] <- c(3, 2, 4, 1)
for(j in 1 : 10) {
	print(j)
	print(pr.load(dataNYC[[4]][, -1, j], 8))
} #5, 7
matchmat[(4 : 13), ] <- matrix(rep(c(2, 3, 4, 1), 10), ncol = 4, byrow = T)
matchmat[3 + 5, ] <- c(2, 3, 4, 1)
matchmat[3 + 7, ] <- c(2, 3, 4, 1)


# #check 5 and 7
mismatch <- test.equal(matchmat, source.prof.all, source.conc.all)
matchmat[mismatch, ]
source.prof.all[mismatch, ]
source.conc.all[mismatch, ]
#just issues with soil in iterations 9 and 10 (6 and 7)
#agree with matchmat
####
#get rid of poor matching
# apca.nyc[[4]] <- apca.nyc[[4]][-mismatch]
source.prof.all[1, ] <- matchmat[1, ]







# #####
# # get scores/prof

scores.prof.apca.nyc <- reorder.fun.apca(apca.nyc, source.prof.all)
scores <- scores.prof.apca.nyc[[1]]
# # ci.fun(scores, trim = 0)
xtable(make.table(scores, trim = 0.1))

# save(scores.prof.apca.nyc, file = file.path(newcode.dir, "apcaNYC_27aug13.RData"))










#################################
#################################
#################################
#################################
#################################
#################################
#########
#Plot of residual oil and traffic


###Plot for paper

setwd(fig.dir)

dates <- read.csv(file.path(dat.dir, "QCII_thurston_date_26nov12.csv"))
datesc <- as.character(dates[,1])
#datesc2 <- datesc[which(substr(datesc, 1, 4) == "2001")]
seqs1 <- seq(1, 48)
datesc2 <- datesc[seqs1]
dateN <- as.Date(datesc2)
atseqs <- dateN[seqs1]


sours <- c("Soil", "Sec. SO2-4", "a. Traffic", "b. Residual Oil")
# sours <- c("Soil", "Sec. SO2-4", "Time series of Traffic in New York City", "Time series of Residual Oil in New York City")
seqs <- seq(1, length(datesc2), length = 10)
ylims <- rbind( c(-5, 20), c(-5, 30), c(-5, 25), c(-5, 10))
library(RColorBrewer)

#b and w
cols <- c(1, "grey60", "grey40", "grey90")
lwds <- c(1, 3, 2, 1)
ltys <- c(1, 1, 3, 1)
pchs <- c(20, 3,  1, 29)
fills <- c(rep("white", 3), "lightcoral")
# fills <- c(rep("white", 3), "thistle3")
# col4 <- "thistle4"


#color
# cols <- c(brewer.pal(4, "Dark2"))
# #lwds <- c(2, 1, 1, 2.5)
# #ltys <- c(1,1, 1, 1)
# fills <- c(rep("white", 3), cols[4])

par( mar = c(8, 4, 4, 2))
# pdf("NY_traffic_resoil_20aug13.pdf", height = 7, width = 12)
pdf("NY_traffic_resoil_color_thesis.pdf", height = 7, width = 12)

for(k in 3 : 4) {

	if(k == 3){
    
		par(mfrow = c(1,2), mar = c(10, 4.5, 4, 1.5))
	}else if (k == 4) {
		par(mar = c(10, 1.5, 4, 4.5))
		
		}


#set up plot
plot(dateN, scores[[1]][seqs1, k], 
     col = cols[1], type = "n",
      ylim =ylims[k, ], main = "", 
     xlab = "", ylab = "", 
     axes = T, lwd = lwds[1], lty = ltys[1], 
     las = 3, cex.axis = 1.3)
     box()
if(k == 3) { 
	mtext(expression(paste("Concentration (",
	mu, "g/m"^"3",")")), side = 2, line = 2.5,
	cex = 1.5)
	}


#likelihood
	# maxlhood <- apply(scores[[4]][seqs1, k, ], 1, max, na.rm = T)
	maxlhood <- apply(scores[[4]][seqs1, k, ], 1, quantile, probs = 0.75, na.rm = T)
	# minlhood <- apply(scores[[4]][seqs1, k, ], 1, min, na.rm = T)
	minlhood <- apply(scores[[4]][seqs1, k, ], 1, quantile, probs = 0.25, na.rm = T)
	ys <- c( minlhood, rev(maxlhood))
	xs <- c(dateN, rev(dateN))
	polygon(xs, ys, col = fills[4], border = NA)
  
  # #plot median
  # meds <- apply(scores[[4]][seqs1, k, ], 1, quantile, probs = 0.5, na.rm = T)
  	# points(dateN, meds, 
     # col = col4, lwd = 2, lty=ltys[1], type = "l")

  
  
  #plot reported
	points(dateN, scores[[1]][seqs1, k], 
     col = cols[1], lwd = lwds[1], lty=ltys[1], type = "l")
	points(dateN, scores[[1]][seqs1, k], 
     col = cols[1], lwd = lwds[1], pch = pchs[1])


	mtext(sours[k], line = 0.5, cex = 1.5)
  

  

#add 1/2MDL and exclude
	points(dateN, scores[[2]][seqs1, k], type = "l", col = cols[2], 
 		lwd = lwds[2], lty = ltys[2])
	points(dateN, scores[[2]][seqs1, k], col = cols[2], pch = pchs[2])
	
	points(dateN, scores[[3]][seqs1, k], type = "l", col = cols[3], 
     	lwd = lwds[3], lty = ltys[3])
	points(dateN, scores[[3]][seqs1,k],col = cols[3], pch = pchs[3])


  if(k != 3){
    swt <- c(4, 2, 1, 3)
    cols2 <- cols[swt]
    legend("topright", col =c("white", cols2[2:4]), 
    	legend = c("Likelihood", "1/2 MDL", "Reported", "Exclude"), 
    	lty = ltys[swt], lwd = lwds[swt]+1, cex = 1.5, pch = pchs[swt])
    d1 <- 30
    y1 <- 8.8
    polygon(as.Date(datesc2[c(d1, d1+7, d1+7, d1)], 
    	origin= "1970-01-01"), c(y1, y1, y1+1.2, y1+1.2), 
    	border = NA, col= fills[4])
graphics.off()
	}


}


















