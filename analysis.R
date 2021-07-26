library(RMySQL)
library(dplyr)
library(gt23)

# Define trials and cell types of interest.
trials <- c('CART19_CLL', 'CART19_ALL', 'CART19_AML')
cellTypes <- c('PBL', 'PBMC', 'Whole blood', 'T-cells')


# Read in sample data & subset to the subjects in this report group.
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 'select * from gtsp')


# Create a list of all GTSPs that passed through the INSPIIRED pipeline.
dbConn  <- dbConnect(MySQL(), group='intsites_miseq')
intSitesamples <- unname(unlist(dbGetQuery(dbConn, 'select sampleName from samples where sampleName like "%GTSP%"')))
intSitesamples <- unique(gsub('\\-\\d+$', '', intSitesamples))


# Standardize cell type namimng.
samples[grepl('whole blood', samples$CellType, ignore.case = TRUE),]$CellType <- 'Whole blood'
samples[grepl('t[\\-\\s]cells', samples$CellType, ignore.case = TRUE, perl = TRUE),]$CellType <- 'T-cells'


# Subset samples to only include with intSite data and those found in the defined trial and cell types.
samples <- subset(samples, CellType %in% cellTypes)
samples <- subset(samples, SpecimenAccNum %in% intSitesamples)
samples <- subset(samples, Trial %in% trials)

if(! file.exists('intSites.rds')){
  intSites <- getDBgenomicFragments(samples$SpecimenAccNum, 'specimen_management', 'intsites_miseq') %>%
              stdIntSiteFragments() %>%
              collapseReplicatesCalcAbunds() %>%
              annotateIntSites()
} else {
  intSites <- readRDS('intSites.rds')
}
  
# Update intSite metadata.
intSites$cellType2 <- samples[match(intSites$GTSP, samples$SpecimenAccNum),]$CellType
intSites$trial <- samples[match(intSites$GTSP, samples$SpecimenAccNum),]$Trial

saveRDS(intSites, 'intSites.rds')

r <- dplyr::group_by(data.frame(intSites), GTSP) %>%
     dplyr::mutate(sitesInSample = n_distinct(posid)) %>%
     dplyr::ungroup() %>%
     dplyr::filter(trial == 'CART19_CLL' & nearestFeature == 'STAT3') %>%
     dplyr::select(strand, posid, estAbund, timePoint, patient, cellType2, nearestFeature, nearestFeatureDist, sitesInSample) %>%
     dplyr::arrange(patient, cellType2, posid)

table(stringr::str_extract(unique(r$posid), '[\\-\\+]'))
