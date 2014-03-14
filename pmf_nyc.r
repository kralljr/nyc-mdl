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
dat.dir <- "/Users/jennakrall/Dropbox/MDL_sourceapp/nyc_analysis_mdl"

setwd(dat.dir)
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
load("dataNYC_pmf.RData")
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


concpmf <- pmf.nyc[[1]][[1]][, c(1, 2, 3, 4)]
profpmf <- pmf.nyc[[1]][[2]][-1, c(1, 2, 3, 4)]
rownames(profpmf) <- colnames(dataNYCPMF[[1]])[-1]

source.conc.all <- match.sources(pmf.nyc, concpmf, "scores", typeSA = "PMF")
source.prof.all <- match.sources(pmf.nyc, profpmf, 
	 "prof", typeSA = "PMF", cn = cn)

for(i in 2 : 3) {
	temp <- pmf.nyc[[i]][[2]][-1, ]
	rownames(temp) <- colnames(dataNYCPMF[[i]])[-1]
	print(round(temp, 3))
}

cn2 <- c("aluminum", "calcium", "titanium", "silicon", "ammonium_ion",
	"sulfur", "OC", "elemental_carbon", "potassium", "chlorine", "lead",
	"nickel", "nitrate", "vanadium", "zinc")
for(i in 1 : 10) {
	print(i)
	temp <- pmf.nyc[[4]][[i]][[2]][-1, ]
	rownames(temp) <- colnames(dataNYCPMF[[4]][[1]])[-1]
	print(round(temp[cn2, ], 3))
}


matchmat <- matrix(nrow = 13, ncol = 4)
matchmat[1, ] <- c(1, 2, 3, 4)
matchmat[2, ] <- c(6, 2, 3, 4)
matchmat[3, ] <- c(7, 4, 6, 3)
matchmat[4, ] <- c(6, 4, 5, 1)
matchmat[5, ] <- c(1, 2, 3, 8)
matchmat[6, ] <- c(6, 4, 3, 8)
matchmat[7, ] <- c(1, 4, 3, 8)
matchmat[8, ] <- c(1, 2, 4, 5)
matchmat[9, ] <- c(6, 2, 4, 8)
matchmat[10, ] <- c(6, 4, 3, 8)
matchmat[11, ] <- c(6, 4, 3, 5)
matchmat[12, ] <- c(3, 6, 5, 8)
matchmat[13, ] <- c(6, 4, 3, 8)



scores.prof.pmf.nyc <- reorder.fun.pmf(pmf.nyc, matchmat)
scores <- scores.prof.pmf.nyc[[1]]
xtable(make.table(scores, trim = 0.1))











########
# Plot

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


lev1 <- c("Traffic", "Residual oil", "Soil", "Sec. Sulfate")
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


pdf("NYC_mdlsources_pmf.pdf", height = 10, width = 9)
p
graphics.off()
