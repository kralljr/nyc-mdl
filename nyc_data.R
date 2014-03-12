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


dat <- read.csv(file.path(home.dir, "SpeciationDB_2009_01_06wBlk_csv.csv"), fill = T, 
	stringsAsFactors = F)
spec <- readRDS("/Users/jennakrall/Dropbox/SpatialFA/data/speciation_monitors.rds")
	
	
	
	
	
####
# Clean NYC speciation data
####	
#keep rows for QC
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
other <- "Date"
dat <- dat[, which(saps1 %in% c(cols, cols2) | saps2 %in% other)]

colnames(dat)[1] <- "Date"
saps2 <- sapply(strsplit(colnames(dat), " "), function(x) x[2])
dat <- dat[, -which(saps2 %in% c("Date", "Analysis"))]
dat <- dat[, -which(colnames(dat) %in% c("Potassium ( K+ )", "Sodium ( Na+ )"))]


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
dat <- data.frame(dat[, 1], dat[,  sort(colnames(dat)[-1])])
colnames(dat)[1] <- "Date"



write.csv(dat, file = "nycdat.csv", row.names = F)










#####
# get CSN/MDL data
#####


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

colsspec <- tolower(c(cols, cols2))
colsspec[which(colsspec == "ammonium")] <- "ammonium_ion"
colsspec[which(colsspec == "elemental")] <- "elemental_carbon"
colsspec[which(colsspec == "organic")] <- "OC" 
# colsspec[which(colsspec == "organic")] <- "OCM_K14"  #no non-missing MDLs


specdat <- specdat[, c(1, which(colnames(specdat) %in% colsspec))]
mdl <- mdl[, c(1, which(colnames(mdl) %in% colsspec))]
specdat <- specdat[complete.cases(specdat), ]
mdl <- mdl[complete.cases(mdl), ]
mdl <- mdl[, -1]
mdl <- mdl[, sort(colnames(mdl))]
mdlmax <- apply(mdl, 2, max)
mdlmin <- apply(mdl, 2, min)


write.csv(mdlmax, file = "mdlmax.csv", row.names = F)












#####
# check matching between datasets

#SPECIATION DATA missing a ton for sulfur, why differential missing
dates <- c(dat[, 1], specdat[, 1])
dates <- unique(dates[which(duplicated(dates))])

datmatch <- dat[which(dat[, 1] %in% dates), ]
colnames(datmatch) <- tolower(sapply(strsplit(colnames(datmatch), " "), function(x) x[1]))
colnames(datmatch)[which(colnames(datmatch) == "organic")] <- "OC"
colnames(datmatch)[which(colnames(datmatch) == "elemental")] <- "elemental_carbon"
colnames(datmatch)[which(colnames(datmatch) == "ammonium")] <- "ammonium_ion"
colnames(datmatch)[1] <- "Date"
specmatch <- specdat[which(specdat[, 1] %in% dates), ]

all.equal(sort(colnames(datmatch)), sort(colnames(specmatch)))

datmatch <- datmatch[, -1]
specmatch <- specmatch[, -1]
datmatch <- datmatch[, sort(colnames(datmatch))]
specmatch <- specmatch[, sort(colnames(specmatch))]

all.equal(trunc(datmatch, 2), trunc(specmatch, 2))


#organic carbon, EC, similar, but not same
 cor(datmatch[, 20], specmatch[, 20])
# [1] 0.9990474
cor(datmatch[, 12], specmatch[, 12])
# [1] 0.9855286