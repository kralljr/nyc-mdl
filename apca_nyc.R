# File to provide summary, APCA source apportionment
# results of  QCII speciation monitor in NYC
# 8/20/13
# ========================================================



#load libraries
library(ncdf4)
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












#Obtain adjusted data by applying censoring adjustment methods
################################
################################
################################
################################
################################
################################


#get matrix of mdls
mdls <- matrix(rep(mdlmax, each = nrow(dat)), 
	ncol = length(mdlmax), byrow = F)
colnames(mdls) <- colnames(mdlsboth)	

#summary of bdl values
sort(round(bdlfun(dat, mdlmax)[[2]], 2))


# ####
# # Get adjusted data for NYC
# Reported
# 1/2 MDL
# Exclude
dataNYC <- createDATs(dat, mdlmax)


#run MCMC, get likelihood-based adjusted data
set.seed(1597)
#do not use PM2.5
datmcmc <- runGIBBS(dat[, -c(1, 2)], mdls[, -1], niter = 50000, 
      burnin = 25000, ndraws = 10, seed = 2361)
dataNYC[[4]] <- datmcmc
names(dataNYC)[4] <- "gibbs" 

# save(dataNYC, file = "dataNYC.RData")


# #################################
# #################################
# #################################
# #################################
# #################################
# #################################












# # APPLY APCA to each
################################
################################
################################
################################
################################
################################
####


#remove dates
dataNYCnodate <- dataNYC
for(i in 1 : 3) {
	dataNYCnodate[[i]] <- dataNYCnodate[[i]][, -1]
}
set.seed(99466)
apca.nyc <- apca.all(dataNYCnodate, tot = T)
# save(apca.nyc, file = "apca_nyc.RData")


#try without total
datT <- dataNYCnodate
for(i in 1 : 3) {
	datT[[i]] <- datT[[i]][, -1]
}
datT[[4]] <- datT[[4]][, -1, ]
set.seed(99466)
apca.nyc.nopm <- apca.all(datT, tot = F)



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
#new
source.conc <- apca.nyc[[1]][[1]][, c(2, 3, 8, 1)]
source.prof <- apca.nyc[[1]][[4]][, c(2, 3, 8, 1)]
source.conc.all <- match.sources(apca.nyc, source.conc, "scores")
source.prof.all <- match.sources(apca.nyc, source.prof, "prof")


#1/2MDL and exclude
#5321
#4321


# ## DOUBLE CHECK
matchmat <- matrix(nrow = 13, ncol = 4)
pr.load(dataNYC[[1]][, -c(1, 2)], 8)
matchmat[1, ] <- c(2, 3, 8, 1)
pr.load(dataNYC[[2]][, -c(1, 2)], 8)
matchmat[2, ] <- c(5, 3, 2, 1)
pr.load(dataNYC[[3]][, -c(1, 2)], 8)
matchmat[3, ] <- c(4, 3, 2, 1)

j <- 1
for(j in 1 : 10) {
	print(j)
	print(round(pr.load(dataNYC[[4]][, , j], 8), 2))
} #5, 7
matchmat[4, ] <- c(2, 4, 8, 1)
matchmat[5, ] <- c(2, 3, 4, 1)
matchmat[6, ] <- c(4, 2, 7, 1)
matchmat[7, ] <- c(5, 2, 3, 1)
matchmat[8, ] <- c(8, 2, 4, 1)
matchmat[9, ] <- c(4, 7, 3, 1)
matchmat[10, ] <- c(8, 2, 3, 1)
matchmat[11, ] <- c(5, 2, 4, 1)
matchmat[12, ] <- c(3, 2, 4, 1)
matchmat[13, ] <- c(1, 2, 4, 3)


# #check 5 and 7
mismatch <- test.equal(matchmat, source.prof.all, source.conc.all)
matchmat[mismatch, ]
source.prof.all[mismatch, ]
source.conc.all[mismatch, ]

source.prof.all <- matchmat


#######
#no PM
source.conc <- apca.nyc.nopm[[1]][[1]][, c(2, 3, 8, 1)]
source.prof <- apca.nyc.nopm[[1]][[4]][, c(2, 3, 8, 1)]
source.conc.all <- match.sources(apca.nyc.nopm, source.conc, "scores")
source.prof.all <- match.sources(apca.nyc.nopm, source.prof, "prof")







# #####
# # get scores/prof
scores.prof.apca.nyc <- reorder.fun.apca(apca.nyc, source.prof.all)
scores <- scores.prof.apca.nyc[[1]]
xtable(make.table(scores, trim = 0.1))









#################################
#################################
#################################
#################################
#################################
#################################
#########

###Plot for paper

dates <- dat[, 1]
datesc <- as.character(dates)
seqs1 <- seq(1, 48)
datesc2 <- datesc[seqs1]
dateN <- as.Date(datesc2)
atseqs <- dateN[seqs1]
n  <- length(dateN)

#create dataset
sc <- 0
t1 <- c("", "")
sc2 <- c(0, 0)
t2 <- c("", "")
datetot <- dateN[1]
datetot2 <- dateN[1]
types <- c("Reported", "1/2MDL", "Exclude", "Likelihood")
source <- c("Soil", "Sec. Sulfate", "Traffic", "Residual oil")
#for each source
for(i in 1 : 4) {
	for(j in 1 : 4) {
		print(c(i, j))

		if(i != 4) {
			sc <- c(sc, scores[[i]][seqs1, j])
		}else{
			sc4 <- apply(scores[[i]][seqs1, j, ], 1, 
				quantile, probs = 0.5, na.rm = T)
			sc <- c(sc, sc4)
		}
		others <- cbind(rep(types[i], n), rep(source[j], n))
		t1 <- rbind(t1, others)
		datetot <- c(datetot, dateN)
		if(i == 4) {
			top <- apply(scores[[i]][seqs1, j, ], 1, 
				quantile, probs = 0.75, na.rm = T)
			bottom <- apply(scores[[i]][seqs1, j, ], 1, 
				quantile, probs = 0.25, na.rm = T)	
			others <- cbind(bottom, top)
			sc2 <- rbind(sc2, others)
			
			others <- cbind( rep("Likelihood", n), rep(source[j], n))
			t2 <- rbind(t2, others)
			datetot2 <- c(datetot2, dateN)
			
		}
	}
	
	
}
t1 <- t1[-1, ]
sc <- sc[-1]
t2 <- t2[-1, ]
sc2 <- sc2[-1, ]
datetot <- datetot[-1]
datetot2 <- datetot2[-1]
dat1 <- data.frame(datetot, t1, sc)
colnames(dat1) <- c("Date", "Type", "Source",  "Conc" )
dat2 <- data.frame(datetot2, t2, sc2)
colnames(dat2) <- c("Date", "Type", "Source", "Bottom", "Top" )
dat1$Date <- as.Date(dat1$Date, origin = "1970-01-01")


lev1 <- c("Soil", "Sec. Sulfate", "Traffic", "Residual oil")
dat1$Source <- factor(dat1$Source, levels = lev1)
dat2$Source <- factor(dat2$Source, levels = lev1)
dat1$Type <- factor(dat1$Type, levels = c("Likelihood", "1/2MDL", "Exclude", "Reported"))

pd <- position_dodge(.4)
size1 <- 18
sizeline <- 0.8


s1 <- "Residual oil"
datU <- dat1[which(dat1$Source == s1), ]
datU <- dat1
p <- ggplot() + scale_colour_manual(name="",
                     breaks=c("Likelihood", "1/2MDL", "Exclude", "Reported"),
                     values = c("#CC6666", "grey20", "grey40", "grey70"))  
p <- p + geom_ribbon(data = dat2, aes(x = Date, 
	ymin = Bottom, ymax = Top, group = Type), fill = "red", alpha = 0.2)

p <- p +
	geom_line(data = datU, aes(x = Date, y = Conc, colour = Type, 
		group = Type, linetype = Type))  +
	geom_point(data = datU, aes(x = Date, y = Conc, colour = Type, 
		group = Type, shape = Type)) +
	scale_shape_manual(values=c(NA, 1,16, 3), name = "", 
		 breaks=c("Likelihood", "1/2MDL", "Exclude", "Reported")) + 
	scale_linetype_manual(values=c(1,2, 3, 1), name="", breaks=c("Likelihood", 
		"1/2MDL", "Exclude", "Reported"))
p <- p +    	# scale_colour_grey(name="",
                     # breaks=c("Likelihood", "1/2MDL", "Exclude", "Reported"),
                     # start = 0.1, end = 0.6)  +
      theme(axis.text.y=element_text(size=size1)) +
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size = size1)) +
      scale_x_date(labels = date_format("%m-%Y")) +
      	theme(panel.background = element_blank()) +
	 ylab(expression("Concentration (" * mu * 
		 "g/m"^"3"*")" )) + xlab("") +
		 theme(legend.text=element_text(size=size1))  
p <- p + facet_wrap(~Source, scales = "free_y", ncol = 1) + 
	theme(strip.text.x = element_text(size = size1), axis.title=element_text(size=size1))
# p


pdf("NYC_mdlsources.pdf", height = 10, width = 9)
p
graphics.off()






######## OLD PLOT



# # h + geom_ribbon(aes(ymin=level-1, ymax=level+1))

# sours <- c("Soil", "Sec. SO2-4", "a. Traffic", "b. Residual Oil")
# # sours <- c("Soil", "Sec. SO2-4", "Time series of Traffic in New York City", "Time series of Residual Oil in New York City")
# seqs <- seq(1, length(datesc2), length = 10)
# ylims <- rbind( c(-5, 20), c(-5, 30), c(-5, 25), c(-5, 10))
# library(RColorBrewer)

# #b and w
# cols <- c(1, "grey60", "grey40", "grey90")
# lwds <- c(1, 3, 2, 1)
# ltys <- c(1, 1, 3, 1)
# pchs <- c(20, 3,  1, 29)
# fills <- c(rep("white", 3), "lightcoral")
# # fills <- c(rep("white", 3), "thistle3")
# # col4 <- "thistle4"


# #color
# # cols <- c(brewer.pal(4, "Dark2"))
# # #lwds <- c(2, 1, 1, 2.5)
# # #ltys <- c(1,1, 1, 1)
# # fills <- c(rep("white", 3), cols[4])

# par( mar = c(8, 4, 4, 2))
# # pdf("NY_traffic_resoil_20aug13.pdf", height = 7, width = 12)
# pdf("NY_traffic_resoil_color_thesis.pdf", height = 7, width = 12)

# for(k in 3 : 4) {

	# if(k == 3){
    
		# par(mfrow = c(1,2), mar = c(10, 4.5, 4, 1.5))
	# }else if (k == 4) {
		# par(mar = c(10, 1.5, 4, 4.5))
		
		# }


# #set up plot
# plot(dateN, scores[[1]][seqs1, k], 
     # col = cols[1], type = "n",
      # ylim =ylims[k, ], main = "", 
     # xlab = "", ylab = "", 
     # axes = T, lwd = lwds[1], lty = ltys[1], 
     # las = 3, cex.axis = 1.3)
     # box()
# if(k == 3) { 
	# mtext(expression(paste("Concentration (",
	# mu, "g/m"^"3",")")), side = 2, line = 2.5,
	# cex = 1.5)
	# }


# #likelihood
	# # maxlhood <- apply(scores[[4]][seqs1, k, ], 1, max, na.rm = T)
	# maxlhood <- apply(scores[[4]][seqs1, k, ], 1, quantile, probs = 0.75, na.rm = T)
	# # minlhood <- apply(scores[[4]][seqs1, k, ], 1, min, na.rm = T)
	# minlhood <- apply(scores[[4]][seqs1, k, ], 1, quantile, probs = 0.25, na.rm = T)
	# ys <- c( minlhood, rev(maxlhood))
	# xs <- c(dateN, rev(dateN))
	# polygon(xs, ys, col = fills[4], border = NA)
  
  # # #plot median
  # # meds <- apply(scores[[4]][seqs1, k, ], 1, quantile, probs = 0.5, na.rm = T)
  	# # points(dateN, meds, 
     # # col = col4, lwd = 2, lty=ltys[1], type = "l")

  
  
  # #plot reported
	# points(dateN, scores[[1]][seqs1, k], 
     # col = cols[1], lwd = lwds[1], lty=ltys[1], type = "l")
	# points(dateN, scores[[1]][seqs1, k], 
     # col = cols[1], lwd = lwds[1], pch = pchs[1])


	# mtext(sours[k], line = 0.5, cex = 1.5)
  

  

# #add 1/2MDL and exclude
	# points(dateN, scores[[2]][seqs1, k], type = "l", col = cols[2], 
 		# lwd = lwds[2], lty = ltys[2])
	# points(dateN, scores[[2]][seqs1, k], col = cols[2], pch = pchs[2])
	
	# points(dateN, scores[[3]][seqs1, k], type = "l", col = cols[3], 
     	# lwd = lwds[3], lty = ltys[3])
	# points(dateN, scores[[3]][seqs1,k],col = cols[3], pch = pchs[3])


  # if(k != 3){
    # swt <- c(4, 2, 1, 3)
    # cols2 <- cols[swt]
    # legend("topright", col =c("white", cols2[2:4]), 
    	# legend = c("Likelihood", "1/2 MDL", "Reported", "Exclude"), 
    	# lty = ltys[swt], lwd = lwds[swt]+1, cex = 1.5, pch = pchs[swt])
    # d1 <- 30
    # y1 <- 8.8
    # polygon(as.Date(datesc2[c(d1, d1+7, d1+7, d1)], 
    	# origin= "1970-01-01"), c(y1, y1, y1+1.2, y1+1.2), 
    	# border = NA, col= fills[4])
# graphics.off()
	# }


# }


















