#########
# File to estimate sources in NYC Queens College monitor
# 3/5/14
#########


#On 11/26/12 Data was downloaded from: 
# "http://www.dec.ny.gov/chemical/8888.html"
# File is "speciationdb.zip" with excel file "SpeciationDB_2009_01_06wBlk.xls"


### DO NOT USE
#On 3/5/14, New data was downloaded from: 
# "ftp://ftp.dec.state.ny.us/dar/library/speciationthru4q11.zip"
# File is excel file "PM25 Speciation thru 4q11.xlsx"
# Does not have EC/OC data?



home.dir <- "/Users/jennakrall/Dropbox/MDL_sourceapp/nyc_analysis_mdl"
datdir2 <- "/Users/jennakrall/Dropbox/SpatialFA/data/"
setwd(home.dir)

	
	
	
	
####
# Clean NYC speciation data
####	
#keep rows for QC
dat <- read.csv(file.path(home.dir, "SpeciationDB_2009_01_06wBlk_csv.csv"), fill = T, 
	stringsAsFactors = F)
	
	
WC <- which(dat[, 2] == "QCII")
dat <- dat[c(1, 2, WC), ]

#fix column names
name1 <- as.matrix(dat[1, ])
name2 <- as.matrix(dat[2, ])
cn <- paste(name1, name2)
colnames(dat) <- cn
dat <- dat[-c(1, 2), ]


#correct dates
dat$" Date" <- as.Date(dat$" Date")
yr <- as.numeric(substr(dat$" Date", 1, 4))
dat <- dat[which(yr <= 2002 & yr >= 2001), ]
#get rid of March date
dat <- dat[-1, ]


#airscode for MDL
airs <- as.character(unique(dat$AIRS))

#keep columns
saps1 <- sapply(strsplit(colnames(dat), " "), function(x) x[1])
saps2 <- sapply(strsplit(colnames(dat), " "), function(x) x[2])
#identify columns to keep
#note ammonium is ammonium ion
cols <- c("Barium", "Cadmium", "Chromium", "Cobalt", "Magnesium",
	"Molybdenum", "Sulfur", "Ammonium")
#potassium/sodium not potassium ion/sodium ion	
cols2 <- c("Aluminum", "Arsenic", "Bromine", "Calcium", "Chlorine",
	"Copper",  "Iron", "Potassium","Manganese", "Sodium",
	"Nickel", "Nitrate", "Phosphorus", "Lead", "Selenium",
	"Silicon", "Strontium", "Titanium", "Vanadium", 
	"Zinc", "Elemental",  "Organic" )
other <- c("Date", "PM25")
dat <- dat[, which(saps1 %in% c(cols, cols2) | saps2 %in% other)]

colnames(dat)[1] <- "Date"
saps1 <- sapply(strsplit(colnames(dat), " "), function(x) x[1])
saps2 <- sapply(strsplit(colnames(dat), " "), function(x) x[2])
dat <- dat[, -which(saps2 %in% c("Date", "Analysis"))]
dat <- dat[, -which(colnames(dat) %in% c("Potassium ( K+ )", "Sodium ( Na+ )"))]
dat <- dat[, -which(saps1 %in% c("Collocated", "Official"))]


#concentrations are numeric
for(i in 2 : ncol(dat)) {
	dat[, i] <- as.numeric(dat[, i])
}

#get complete cases 
dat <- dat[complete.cases(dat), ]


#fix column names and reorder
colnames(dat) <- tolower(sapply(strsplit(colnames(dat), " "), function(x) x[1]))
colnames(dat)[which(colnames(dat) == "organic")] <- "OC"
colnames(dat)[which(colnames(dat) == "elemental")] <- "elemental_carbon"
colnames(dat)[which(colnames(dat) == "ammonium")] <- "ammonium_ion"
dat <- data.frame(dat[, c(1, 2)], dat[,  sort(colnames(dat)[-c(1, 2)])])
colnames(dat)[1:2] <- c("Date", "PM25")





write.csv(dat, file = "nycdat.csv", row.names = F)










#####
# get CSN/MDL data
#####

spec <- readRDS(file.path(datdir2, "speciation_monitors.rds"))

airs1 <- paste0(substr(airs, 1, 5), ".", substr(airs, 6, 9))
spec <- spec[[which(names(spec) == airs1)]]

#restrict dates
yr <- as.numeric(substr(spec$Date, 1, 4))
spec1 <- spec[which(yr >= 2001 & yr <= 2002), ]
mon <- as.numeric(substr(spec$Date, 6, 7))
spec1 <- spec1[which((mon >= 4 & yr == 2001) | yr == 2002), ]


#separate MDL matrix from data matrix
saps2 <- sapply(strsplit(colnames(spec1), "\\."), function(x) x[2])
mdl <- spec1[, c(1, which(saps2 == "MDL"))]
colnames(mdl) <- sapply(strsplit(colnames(mdl), "\\."), function(x) x[1])
colnames(mdl)[1] <- "Date"
specdat <- spec1[, which(is.na(saps2))]


#fix column names
colsspec <- c(tolower(c(cols, cols2)), "PM25_SPEC")
colsspec[which(colsspec == "ammonium")] <- "ammonium_ion"
colsspec[which(colsspec == "elemental")] <- "elemental_carbon"
colsspec[which(colsspec == "organic")] <- "OC" 
# colsspec[which(colsspec == "organic")] <- "OCM_K14"  #no non-missing MDLs


#fix order of columns, fix column names
specdat <- specdat[, c(1, which(colnames(specdat) %in% colsspec))]
pm25 <- specdat$PM25_SPEC
specdat <- specdat[, -which(colnames(specdat) == "PM25_SPEC")]
specdat <- data.frame(specdat[, 1], pm25, specdat[, sort(colnames(specdat)[-1])])
colnames(specdat)[1:2] <- c("Date", "PM25")


#fix MDLS
mdl <- mdl[, c(1, which(colnames(mdl) %in% colsspec))]
specdat <- specdat[complete.cases(specdat), ]
mdl <- mdl[complete.cases(mdl), ]
mdl <- mdl[, -1]
pm25 <- mdl$PM25_SPEC
mdl <- mdl[, -which(colnames(mdl) == "PM25_SPEC")]
mdl <- data.frame(pm25, mdl[, sort(colnames(mdl))])
colnames(mdl)[1] <- "PM25"
mdlmax <- apply(mdl, 2, max)
mdlmin <- apply(mdl, 2, min)


mdls <- rbind(mdlmax, mdlmin)

write.csv(mdls, file = "mdls.csv", row.names = F)






#####
# get avg %mdl for each source
bdls <- 1 * (sweep(dat[, -1], 2, mdlmax, "<"))
pbdls <- round(apply(bdls, 2, mean), 2)

sour1 <- list()
sour1[[1]] <- pbdls[c("aluminum", "calcium", "iron", "silicon", "titanium")]
sour1[[2]]  <- pbdls[c("sulfur", "ammonium_ion")]
sour1[[3]]  <- pbdls[c("OC", "elemental_carbon", "potassium")]
sour1[[4]]  <- pbdls[c("nickel", "vanadium", "nitrate",	
	"chlorine", "lead", "zinc")]
	
meanmed <- function(x) {median(x, na.rm = T)}	
x <- sapply(sour1, meanmed)
names(x) <- c("soil", "sec sulf", "traff", "resoil")
x





#####
# check matching between datasets

#SPECIATION DATA missing a ton for sulfur, why differential missing
dates <- c(dat[, 1], specdat[, 1])
dates <- unique(dates[which(duplicated(dates))])

datmatch <- dat[which(dat[, 1] %in% dates), ]
specmatch <- specdat[which(specdat[, 1] %in% dates), ]

all.equal(sort(colnames(datmatch)), sort(colnames(specmatch)))

datmatch <- datmatch[, -1]
specmatch <- specmatch[, -1]
all.equal(trunc(datmatch, 2), trunc(specmatch, 2))


#organic carbon, EC, similar, but not same
 cor(datmatch$elemental_carbon, specmatch$elemental_carbon)
# [1] 0.9852339
cor(datmatch$OC, specmatch$OC)
# [1] 0.9990565