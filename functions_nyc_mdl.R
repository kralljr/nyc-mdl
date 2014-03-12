
###########
#### File for different functions for MDL analysis of NYC data
# 8/16/13











#########################
#Function to give summary of BDL amounts
# dat is matrix of data with rows as 
# observation days and columns as constituents
# mdls is vector of mdls for each constituent
# date is whether data contains column of dates
bdlfun <- function(dat, mdls, date = T) {


  #remove date column
  if(date) {
  	dat <- dat[, -1]
  }
  
  #find which values BDL
  bdls <- 1 * (sweep(dat, 2, mdls, "<"))

  #% BDL for each constituent
  cm <- colMeans(bdls, na.rm = T)
  print(summary(cm))
  names(cm) <- colnames(dat)
  print(table(cut(cm, breaks = 5)))
  
  #Print different data summaries
  cat("Number have BDL with <5% = ", length(which(cm < 0.05)), 
  	"/", ncol(bdls), "\n")
  cat("Num BDL =", 
        length(which(colSums(bdls, na.rm = T)> 0)), 
        "/", ncol(bdls), "\n")
  cat("Num BDL >25% =", length(which(cm > 0.25)), 
  	"/", ncol(bdls), "\n")  
  cat("Num BDL >50% =", length(which(cm > 0.50)), 
  	"/", ncol(bdls), "\n")  	
  cat("Name with BDL >50%:", colnames(dat)[which(cm>.5)] ,"\n")
  
  list(bdls, cm)
  }




##############
# Function to create different MDL datasets	
# dat is matrix of data with rows as 
# observation days and columns as constituents
# mdls is vector of mdls for each constituent
# date is whether data contains column of dates
createDATs <- function(dat, mdls, date = T) {	
	
	#whether to drop column of dates
  if(date) {
  	dates <- dat[, 1]
  	dat <- dat[, -1]
  	
  }
  
  #find which values BDL
  bdls <- 1 * (sweep(dat, 2, mdls, "<"))
  #summary of BDL values
  cm <- colMeans(bdls, na.rm = T)
  
  #1/2mdl
  first <- sweep(bdls == 1, 2, 1/2 * mdls, "*")
  dat.5mdl <- first + dat * (bdls == 0)
  
  #drop cons
  w.25 <- which(cm > .25)  
  if ( length( w.25) > 0 ) {
    dat.drop <- dat.5mdl[, -w.25]
  }else{
  	dat.drop <- dat.5mdl
  	}
  	
  #add back in dates
  dat <- data.frame(dates, dat)
  dat.5mdl 	<- data.frame(dates, dat.5mdl)
  dat.drop <- data.frame(dates, dat.drop)

	list(dat = dat, dat5 = dat.5mdl, datdrop = dat.drop)
}





##########
# Function to apply APCA to each of 4 adjustments
apca.all <- function(datlist, tot, nf = 8) {
	
	apcares <- vector(mode = "list", length = 4)
	for(i in 1:3) {
		if(tot == T) {
			dat1 <- datlist[[i]][, -1]
			tots <- datlist[[i]][, 1]
		}else{
			dat1 <- datlist[[i]]
			tots <- rowSums(dat1)
			}
		apcares[[i]] <- abspca(dat = dat1, tot= tots, nfactors = nf)
	}
	
	
	apcares[[4]] <- vector(mode = "list", length = dim(datlist[[4]])[3])
	for(i in 1 : dim(datlist[[4]])[3]) {
		if(tot == T) {
			dat1 <- datlist[[4]][, -1, i]
			tots <- datlist[[1]][, 1]
		}else{
			dat1 <- datlist[[4]][, , i]
			tots <- rowSums(dat1)
			}
		apcares[[4]][[i]] <- abspca(dat = dat1, 
			tot = tots, nfactors = nf)
		}
	
	names(apcares) <- c("dat", "dat5", "datdrop", "gibbs")
	apcares
}
  
  
  
  
  
#######
# Function to run gibbs imputation
runGIBBS <- function(dat, mdls, niter = 100, 
  burnin = 50, ndraws = 3, seed, 
  homedir = "/Users/jennakrall/Dropbox/MDL_sourceapp",
  wdraw = NULL, cov = NULL) {

  gw <- getwd()
  setwd(file.path(homedir, "git-gibbs-truncnorm"))
  
  write.csv(dat, file = "dat.csv", row.names = F)
  write.csv(mdls, file = "mdls.csv", row.names= F)	
  
  outfilename <- "test.nc"
	comm <- paste("./gibbs -n", niter, "-b", burnin,
		"-d", ndraws, "-r", seed, "dat.csv", "mdls.csv", outfilename)
		
  system(comm)	 
	
  data <- nc_open(outfilename)
  data2 <- ncvar_get(data, "data")
  data2 <- aperm(data2, c(2,1,3))
  
    dimnames(data2) <- list(seq(1, nrow(dat)), 
                          colnames(dat),
                          seq(1, ndraws))
  
  out <- data2
  
  if (is.null(cov) != TRUE) {
    cov <- ncvar_get(data, "covariance")
    mean <- ncvar_get(data, "mean")
    
    datcov <- findtln(dat, mdls, cov, mean, nreps = ndraws)
    colnames(datcov) <- colnames(dat)
    out <- datcov
    } 
  
  nc_close(data)

  
  setwd(gw)
  out
}





####
make.list <- function(list, type) {
	
	if(type == "scores") {
		index <- 1
	}else if(type == "prof"){
		index  <- 4
	}
	
	out <- list()
	for(i in 1 : 3) {
		out[[i]] <- list[[i]][[index]]
	}
	
	#fix gibbs
	out[[4]] <- list()
	for(j in 1 : length(list[[4]])) {
		out[[4]][[j]] <- list[[4]][[j]][[index]]
	}
	
	
	out
}




####
make.list.pmf <- function(list, type, cn) {
	
	if(type == "scores") {
		index <- 1
	}else if(type == "prof"){
		index  <- 2
	}
	
	out <- list()
	for(i in 1 : 3) {
		cn1 <- cn
		if(i == 3) {
			cn1 <- colnames(dataNYCPMF[[3]])[-1]
		}
		out[[i]] <- list[[i]][[index]]
		if(index == 2) {
			rownames(out[[i]]) <- cn1
		}
		
	}
	
	#fix gibbs
	out[[4]] <- list()
	for(j in 1 : length(list[[4]])) {
		out[[4]][[j]] <- list[[4]][[j]][[index]]
		if(index == 2) {
			rownames(out[[4]][[j]]) <- cn
		}
	}
	
	
	out
}



####
find.top.mtch <- function(sares1, sares2, type, holds = 0) {
	
	match <- vector()
	
	if(type == "scores") {
		cormat <- cor(sares1, sares2)
		for(i in 1 : nrow(cormat)) {
			match[i] <- which.max(cormat[i, ])
		}
		compare <- cormat
		
	}else if(type == "prof"){
		#if cons are dropped
		if(nrow(sares1) != nrow(sares2)) {
			sares1 <- sares1[which(rownames(sares1) %in% rownames(sares2)) ,]
		}
		
		dista <- dist(t(cbind(sares1, sares2)))
		ncol1 <- ncol(sares1)
		ncol2 <- ncol(sares2)
		dista <- as.matrix(dista)[1 : ncol1, (ncol1 + 1) : (ncol1 + ncol2)]
		for(i in 1 : nrow(dista)) {
			match[i] <- which.min(dista[i, ])
		}
		compare <- dista
	}	
	
	if(length(unique(match)) != ncol(sares1)) {
		# browser()
	}
	
	
	if(holds != 0) {
		print(compare)
		print(match)
	}
	print(match)
	match
}



####
# match sources
match.sources <- function(saresults, report, type, ncol1 = NULL,
	holds = 0, typeSA = "APCA", cn = NULL) {
	if(typeSA == "APCA") {
		scores <- make.list(saresults, type)
	}else{
		scores <- make.list.pmf(saresults, type, cn)
		}
	
	if(is.null(ncol1)) {
		ncol1 <- ncol(scores[[2]])
	}
	nreps <- length(scores[[4]])
	matchmat <- matrix(nrow = nreps + 3, ncol = ncol(report))
	
	for(i in 2 : 3) {

		ncol1 <- ncol(scores[[i]])
		mtch <- find.top.mtch(report, scores[[i]][, 1 : ncol1], type = type)
		matchmat[i, ] <- mtch
	}
	
	curr <- 4
	for(j in 1 : nreps) {
		if(j %in% holds) {
			holds1 <- j
		}else{holds1 <- 0}
		
		ncol1 <- ncol(scores[[4]][[j]])
		
		print(j)
		mtch <- find.top.mtch(report, scores[[4]][[j]][, 1 : ncol1], type = type,
			holds = holds1)
		matchmat[curr, ] <- mtch
		
		curr <- curr + 1
	}
	
	
	matchmat
}



###
# principal loadings
pr.load <- function(data, nfactors) {
	datasc <- stdize1(data, sd = 1)

	pc1 <- principal(datasc, nfactors, scores = TRUE, rotate = "varimax")
	pc1
}



####
# 
test.equal <- function(match, prof, conc) {
	mismatch <- 0
	
	for(i in 2 : nrow(match)) {
		t1 <- all.equal(match[i, ], prof[i, ])
		t2 <- all.equal(match[i, ], conc[i, ])
		t3 <- all.equal(conc[i, ], prof[i, ])
		
		
		if( !(t1 == TRUE & t2 == TRUE & t3 == TRUE)) {
			mismatch <- c(mismatch, i)
			
		}
	}
	mismatch[-1]
	
}
  
  
  
  
####
#reorder
reorder.fun.apca <- function(sares, ordermat) {
	
	scores <- list()
	prof <- list()
	for(i in 1 : 3) {
		scores[[i]] <- sares[[i]][[1]][, ordermat[i, ]]
		prof[[i]] <- sares[[i]][[4]][, ordermat[i, ]]
	}
	
	curr <- 4
	
	nreps <- length(sares[[4]])
	scores[[4]] <- array(dim = c(nrow(scores[[1]]), ncol(scores[[1]]), nreps))
	prof[[4]] <- array(dim = c(nrow(prof[[1]]), ncol(prof[[1]]), nreps))
	for(j in 1 : nreps) {
		scores[[4]][ , , j] <- sares[[4]][[j]][[1]][, ordermat[curr, ]]
		prof[[4]][, , j] <- sares[[4]][[j]][[4]][, ordermat[curr, ]]
		curr <- curr + 1
		
	}
	
	list(scores, prof)
	
}  



reorder.fun.pmf <- function(sares, ordermat) {
	
	scores <- list()
	prof <- list()
	for(i in 1 : 3) {
		scores[[i]] <- sares[[i]][[1]][, ordermat[i, ]]
		prof[[i]] <- sares[[i]][[2]][, ordermat[i, ]]
	}
	
	curr <- 4
	
	nreps <- length(sares[[4]])
	scores[[4]] <- array(dim = c(nrow(scores[[1]]), ncol(scores[[1]]), nreps))
	prof[[4]] <- array(dim = c(nrow(prof[[1]]), ncol(prof[[1]]), nreps))
	for(j in 1 : nreps) {

		temp1 <- as.matrix(sares[[4]][[j]][[1]])
		temp2 <- as.matrix(sares[[4]][[j]][[2]])
		scores[[4]][ , , j] <- temp1[, ordermat[curr, ]]
		prof[[4]][, , j] <- temp2[, ordermat[curr, ]]
		curr <- curr + 1
		
	}
	
	list(scores, prof)
	
}  


####
# 
make.table <- function(scores, trim = 0) {
	nsources <- ncol(scores[[1]])
	means <- matrix(nrow = 4, ncol = nsources)
	sds <- means
	
	lhood.mean <- apply(scores[[4]], c(1, 2), mean, na.rm = T, trim = trim)
	
	
	for(i in 1 : 4) {
		if(i < 4) {
			dat <- scores[[i]]
		}else{
			dat <- lhood.mean
			}
		
		means[i, ] <- colMeans(dat, na.rm = T)
		sds[i, ] <- apply(dat, 2, sd, na.rm = T)
	}
	
	#NO
	# mn <- matrix(nrow = dim(scores[[4]])[3], ncol = nsources)
	# sd1 <- mn
	# for(j in 1 : dim(scores[[4]])[3]) {
		# mn[j, ] <- colMeans(scores[[4]][,, j], na.rm = T, trim = trim)
		# sd1[j, ] <- apply(scores[[4]][,, j], 2, sd, na.rm = T)
	# }
	# ss <- sample(seq(1, nrow(mn)), 1)
	# means[4, ] <- apply(mn, 2, meanmd, trim = 0.4, na.rm = T, ss = ss)
	# sds[4, ] <- apply(sd1, 2, meanmd, trim = 0.2, na.rm = T, ss = ss)
	
	

	

	
	outtab <- means
	for(i in 1 : nrow(means)) {
		for(j in 1 : ncol(means)) {
			outtab[i, j] <- paste0(round(means[i, j], 2), 
				" (", round(sds[i, j], 2), ")")
			
		}
	}
	rownames(outtab) <- c("Reported", "1/2 MDL", "Exclude", "Likelihood")
	
	outtab[c(1, 4, 2, 3), ]
}



meanmd <- function(vec, na.rm = T, trim = 0, ss = 0) {
	# mean(vec, trim = trim, na.rm = na.rm)
	# median(vec, na.rm = na.rm)
	vec[ss]
	
}


find.trim <- function(scores, trim = 0) {
	qL <- apply(scores, c(2, 3), quantile, probs = trim, na.rm = T)
	qH <- apply(scores, c(2, 3), quantile, probs = 1 - trim, na.rm = T)
	
	
	means <- matrix(nrow = dim(scores)[2], ncol = dim(scores)[3])
	vars <- means
	for(i in 1 : dim(scores)[2]) {
		for(j in 1 : dim(scores)[3]) {
			wh.new <- which(scores[, i, j] >= qL[i, j] & 
				scores[, i, j] <= qH[i, j])
			
			means[i, j] <- mean(scores[wh.new, i, j])
			usemean <- mean(scores[, i, j], 
				trim = trim, na.rm = T)
			if(all.equal(means[i, j], usemean) == F) {browser()}
			vars[i, j] <- var(scores[wh.new, i, j])/length(wh.new)
		}
	}
	list(means, vars)
}


gibbs.se <- function(scores, trim = 0) {
	out <- find.trim(scores, trim = trim)
	means <- out[[1]]
	sesq <- out[[2]]
	
	m <- ncol(means)
	qbar <- apply(means, 1, mean, na.rm = T)
	sesqbar <- apply(sesq, 1, mean, na.rm = T)
	b <- 1/(m - 1) * rowSums(sweep(means, 1, qbar, "-") ^ 2)
	t <- sesqbar + (1 + 1/m) * b
	# df <- (m - 1) * (1 + m * sesqbar / ((m + 1) * b))
	# tstat <- qt(0.975, df)
	# lb <- round(qbar - tstat * sqrt(t), 2)
	# ub <- round(qbar + tstat * sqrt(t), 2)
	paste0(round(qbar, 2), " (", round(sqrt(t), 2), ")")
}


ci.fun <- function(scores, trim = 0) {
	ci.out <- matrix(nrow = 4, ncol = 4)
	
	for(i in 1 : 3) {
		
		means <- colMeans(scores[[i]], na.rm = T)
		ses <- apply(scores[[i]], 2, sd) / sqrt(nrow(scores[[i]]))
		# lb <- means - qt(0.975, (nrow(scores[[i]]) - 1)) * ses
		# ub <- means + qt(0.975, (nrow(scores[[i]]) - 1)) * ses
		# cis <- paste0(means, "(", lb, ", ", ub, ")")
		cis <- paste0(round(means, 2), " (", round(ses, 2), ")")
		ci.out[i, ] <- cis
	}
	
	ci.out[4, ] <- gibbs.se(scores[[4]], trim = trim)
	rownames(ci.out) <- c("Reported", "1/2 MDL", "Exclude", "Likelihood")
	
	
	ci.out[c(1, 4, 2, 3), ]
	
}




