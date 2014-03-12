#########
# File to estimate sources in NYC Queens College monitor
# 3/5/14
#########

#On 3/5/14, Data was downloaded from: 
# "ftp://ftp.dec.state.ny.us/dar/library/speciationthru4q11.zip"
# File is excel file "PM25 Speciation thru 4q11.xlsx"

library(xlsx)
home.dir <- "/Users/jennakrall/Dropbox/MDL_sourceapp/nyc_analysis_mdl"

# dat <- read.xlsx(file.path(home.dir, "PM25 Speciation thru 4q11.xlsx"), 
	# sheetName = "Queens College")
dat <- read.csv(file.path(home.dir, "PM25 Speciation thru 4q11_queenscsv.csv"), fill = T)
dat <- dat[-c(1, 3), ]
dat <- dat[-c(2 : (which(dat[, 1] == "4/7/01") - 1)), ]

#fix date
dat[, 1] <- as.Date(dat[, 1], origin = "1970-01-01", format = "%m/%d/%y")


d1 <- as.matrix(dat[1, ])
colnames(dat) <- as.character(d1)
dat <- dat[-1, ]
colnames(dat)[1] <- "Date"

#remove NA colnames
whNA <- which(is.na(colnames(dat)))
dat <- dat[, -whNA]


#fix column names
saps1 <- sapply(strsplit(colnames(dat), " "), function(x) x[2])
dat <- dat[, c(1, which(saps1 == "Pm2.5" ),
	 which(saps1 == "Ion"), which(saps1 == "Nitrate"), 
	 which(saps1 == "denum"))]



#identify columns to keep
#note ammonium is ammonium ion
cols <- c("Barium", "Cadmium", "Chromium", "Cobalt", "Magnesium",
	"Molyb-", "Sulfur", "Ammonium")
#potassium/sodium not potassium ion/sodium ion	
cols2 <- c("Aluminum", "Arsenic", "Bromine", "Calcium", "Chlorine",
	"Copper",  "Iron", "Potassium","Manganese", "Sodium",
	"Nickel", "Nitrate", "Phosphorus", "Lead", "Selenium",
	"Silicon", "Strontium", "Titanium", "Vanadium", 
	"Zinc", "EC",  "OC" )





