source("preprocessing_scripts/utils.R")

# Load morphological dataset
morphologicalDF <- read.table( paste0(config::get("data"), "morphologicalDF.tsv"), header=T, row.names=1)
stimuliDF <- read.table( paste0(config::get("data"), "stimuliDF.tsv"), header=T)

stimEnd <- 900 # Recommendation, last stim is about at 890 and then 300s later for spontaneous recovery.
stimStart <- min(stimuliDF$time) # When did the first tap start?
burnStart <- 300 # Recommendation
burnEnd <- stimStart # Equivalent.

# Identify timepoints of interest
TIME_SAFE_2ND = 605
TIME_SAFE_30TH = (max(stimuliDF$time) - 100)


intersection <- intersect(names(morphologicalDF), names(stimuliDF))
innerDF_SRC <- merge(morphologicalDF, stimuliDF, by=intersection, all.x=F, all.y=T)
burnIn <- morphologicalDF[morphologicalDF$time >= burnStart & morphologicalDF$time < burnEnd, ]

# Add ratios to innerDF
innerDF <- computeProbabilities(innerDF_SRC) # Important: Run this before dropRedundant ---v
innerDF <- dropRedundant(innerDF) # Important: Run this after computeRatios  ---^

# How many genes do we have
genes <- unique(innerDF$gene)
lines <- unique(innerDF$line)
plates <- unique(innerDF$plate)

# 1) Rescale ALL the data together
# innerDF.sc <- rescale(innerDF)
# burnIn.sc  <- rescale(burnIn) 

# 2) Create subsets
# after.sc <- innerDF.sc[innerDF.sc$time >= stimStart & innerDF.sc$time < stimEnd, ]
# spontaneous.sc <- innerDF.sc[innerDF.sc$time >= stimEnd, ]
# initial.sc <- getFirstN(after.sc, 1)
# second.sc <- getFirstN(after.sc[after.sc$time > TIME_SAFE_2ND,], 1)
# final3.sc <- getLastN(after.sc[after.sc$time <  TIME_SAFE_30TH,], 3)

after <- innerDF[innerDF$time >= stimStart & innerDF$time < stimEnd, ]
spontaneous <- innerDF[innerDF$time >= stimEnd, ]
initial <- getFirstN(after, 1)
second <- getFirstN(after[after$time > TIME_SAFE_2ND,], 1)
final3 <- getLastN(after[after$time <  TIME_SAFE_30TH,], 3)

# 3) a) Reduce each PLATE to a single value

# FIXME; This is useless now
after.cm <- applyColMeansByPlate( after  )
burnIn.cm <- applyColMeansByPlate( burnIn  )
initial.cm <- applyColMeansByPlate( initial  )
second.cm <- applyColMeansByPlate( second  )
final3.cm <- applyColMeansByPlate( final3 )
spontaneous.cm <- applyColMeansByPlate( spontaneous )

after.cm.mtrx <-  collapseToMatrix(after.cm)
burnIn.cm.mtrx <- collapseToMatrix(burnIn.cm)
spontaneous.cm.mtrx <- collapseToMatrix(spontaneous.cm)
initial.cm.mtrx <- collapseToMatrix(initial.cm)
final3.cm.mtrx <- collapseToMatrix(final3.cm)
second.cm.mtrx <- collapseToMatrix(second.cm)

after.cs.mtrx <-  collapseToMatrix(after.cm, FUN=sd_rm.na)
burnIn.cs.mtrx <- collapseToMatrix(burnIn.cm, FUN=sd_rm.na)
spontaneous.cs.mtrx <- collapseToMatrix(spontaneous.cm, FUN=sd_rm.na)
initial.cs.mtrx <- collapseToMatrix(initial.cm, FUN=sd_rm.na)
final3.cs.mtrx <- collapseToMatrix(final3.cm, FUN=sd_rm.na)
second.cs.mtrx <- collapseToMatrix(second.cm, FUN=sd_rm.na)

# b) Do the same but keep plate separate
after.pcm.mtrx <-  collapsePlatesToMatrix(after.cm)
burnIn.pcm.mtrx <- collapsePlatesToMatrix(burnIn.cm)
spontaneous.pcm.mtrx <- collapsePlatesToMatrix(spontaneous.cm)
initial.pcm.mtrx <- collapsePlatesToMatrix(initial.cm)
final3.pcm.mtrx <- collapsePlatesToMatrix(final3.cm)
second.pcm.mtrx <- collapsePlatesToMatrix(second.cm)

after.pcs.mtrx <-  collapsePlatesToMatrix(after.cm, FUN=sd_rm.na)
burnIn.pcs.mtrx <- collapsePlatesToMatrix(burnIn.cm, FUN=sd_rm.na)
spontaneous.pcs.mtrx <- collapsePlatesToMatrix(spontaneous.cm, FUN=sd_rm.na)
initial.pcs.mtrx <- collapsePlatesToMatrix(initial.cm, FUN=sd_rm.na)
final3.pcs.mtrx <- collapsePlatesToMatrix(final3.cm, FUN=sd_rm.na)
second.pcs.mtrx <- collapsePlatesToMatrix(second.cm, FUN=sd_rm.na)


dim(after.pcm.mtrx)
dim(after.pcs.mtrx)
dim(burnIn.pcm.mtrx)
dim(burnIn.pcs.mtrx)
dim(spontaneous.pcs.mtrx)
dim(spontaneous.pcm.mtrx)


# Fixing the missing row
temprow <- matrix(c(rep.int(NA,length(data))),nrow=1,ncol=ncol(spontaneous.pcm.mtrx))
oddOneOut.row <-data.frame(temprow)
oddOneOut.name <- setdiff(rownames(burnIn.pcm.mtrx), rownames(spontaneous.pcm.mtrx))
colnames(oddOneOut.row) <- colnames(spontaneous.pcm.mtrx)
rownames(oddOneOut.row) <- c(oddOneOut.name)
spontaneous.pcm.mtrx <- rbind(spontaneous.pcm.mtrx, oddOneOut.row)[rownames(burnIn.pcm.mtrx), ]
nrow(spontaneous.pcm.mtrx)
table(rownames(spontaneous.pcm.mtrx) == rownames(burnIn.pcm.mtrx))

nrow(initial.pcm.mtrx)
nrow(final3.pcm.mtrx)
nrow(second.pcm.mtrx)


############################################################
#  Means 

reversals <- 
  data.frame(
    initRespPb = initial.cm.mtrx[, c("probSR") ] ,
    initRespDuration = initial.cm.mtrx[, c("revDuration") ] ,
    initRespSpeed = initial.cm.mtrx[, c("speed") ] ,
    initRespDistance = initial.cm.mtrx[, c("revDistance") ] ,
    
    habituationPb = initial.cm.mtrx[, c("probSR") ] - final3.cm.mtrx[, c("probSR") ] ,
    habituationDuration = initial.cm.mtrx[, c("revDuration") ] - final3.cm.mtrx[, c("revDuration") ] ,
    habituationSpeed = initial.cm.mtrx[, c("speed") ]  - final3.cm.mtrx[, c("speed") ] ,
    habituationDistance = initial.cm.mtrx[, c("revDistance") ] - final3.cm.mtrx[, c("revDistance") ] ,
    
    spontPb = spontaneous.cm.mtrx[, c("probSR") ],
    spontDuration = spontaneous.cm.mtrx[, c("revDuration") ] ,
    spontSpeed = spontaneous.cm.mtrx[, c("speed") ] ,
    spontDistance = spontaneous.cm.mtrx[, c("revDistance") ] 
  )


reversals.plate <- 
  data.frame(
    initRespPb = initial.pcm.mtrx[, c("probSR") ] ,
    initRespDuration = initial.pcm.mtrx[, c("revDuration") ] ,
    initRespSpeed = initial.pcm.mtrx[, c("speed") ] ,
    initRespDistance = initial.pcm.mtrx[, c("revDistance") ] ,
    
    habituationPb = initial.pcm.mtrx[, c("probSR") ] - final3.pcm.mtrx[, c("probSR") ] ,
    habituationDuration = initial.pcm.mtrx[, c("revDuration") ] - final3.pcm.mtrx[, c("revDuration") ] ,
    habituationSpeed = initial.pcm.mtrx[, c("speed") ]  - final3.pcm.mtrx[, c("speed") ] ,
    habituationDistance = initial.pcm.mtrx[, c("revDistance") ] - final3.pcm.mtrx[, c("revDistance") ] ,
    
    spontPb = spontaneous.pcm.mtrx[, c("probSR") ],
    spontDuration = spontaneous.pcm.mtrx[, c("revDuration") ] ,
    spontSpeed = spontaneous.pcm.mtrx[, c("speed") ] ,
    spontDistance = spontaneous.pcm.mtrx[, c("revDistance") ] 
  )

############################################################
# Standard deviation

reversals.sd <- 
  data.frame(
    initRespPb = initial.cs.mtrx[, c("probSR") ] ,
    initRespDuration = initial.cs.mtrx[, c("revDuration") ] ,
    initRespSpeed = initial.cs.mtrx[, c("speed") ] ,
    initRespDistance = initial.cs.mtrx[, c("revDistance") ] ,
    
    habituationPb = sdOfTwoSlices( a=initial.cm, b=final3.cm, VAR = "probSR" ),
    habituationDuration = sdOfTwoSlices( a=initial.cm, b=final3.cm, VAR = "revDuration" ),
    habituationSpeed = sdOfTwoSlices( a=initial.cm, b=final3.cm, VAR = "speed" ),
    habituationDistance = sdOfTwoSlices( a=initial.cm, b=final3.cm, VAR = "revDistance" ),
    
    spontPb = spontaneous.cs.mtrx[, c("probSR") ],
    spontDuration = spontaneous.cs.mtrx[, c("revDuration") ] ,
    spontSpeed = spontaneous.cs.mtrx[, c("speed") ] ,
    spontDistance = spontaneous.cs.mtrx[, c("revDistance") ] 
  )

# reversals.plate.sd <- 
#   data.frame(
#     initRespPb = initial.pcs.mtrx[, c("probSR") ] ,
#     habituationPb = collapseToMatrix(initial.pcm[, c("probSR") ] - final3.pcm[, c("probSR") ], FUN=sd_rm.na) ,
#     initRespDuration = initial.pcs.mtrx[, c("revDuration") ] ,
#     habituationDuration = collapseToMatrix(initial.pcm[, c("revDuration") ] - final3.pcm[, c("revDuration") ] , FUN=sd_rm.na) ,
#     initRespSpeed = initial.pcs.mtrx[, c("speed") ] ,
#     habituationSpeed = collapseToMatrix(initial.pcm[, c("speed") ]  - final3.pcm[, c("speed") ]  , FUN=sd_rm.na) , 
#     initRespDistance = initial.pcs.mtrx[, c("revDistance") ] ,
#     habituationDistance = collapseToMatrix(initial.pcm[, c("revDistance") ] - final3.pcm[, c("revDistance") ]  , FUN=sd_rm.na),
#     spontDuration = spontaneous.pcs.mtrx[, c("revDuration") ] ,
#     spontSpeed = spontaneous.pcs.mtrx[, c("speed") ] ,
#     spontDistance = spontaneous.pcs.mtrx[, c("revDistance") ] 
#   )
# # reversals.plate.sd <- 
# #   data.frame(
# #     initRespPb = initial.pcs.mtrx[, c("probSR") ] ,
# #     habituationPb = initial.pcs.mtrx[, c("probSR") ] - final3.pcs.mtrx[, c("probSR") ] ,
# #     initRespDuration = initial.pcs.mtrx[, c("revDuration") ] ,
# #     habituationDuration = initial.pcs.mtrx[, c("revDuration") ] - final3.pcs.mtrx[, c("revDuration") ] ,
# #     initRespSpeed = initial.pcs.mtrx[, c("speed") ] ,
# #     habituationSpeed = initial.pcs.mtrx[, c("speed") ]  - final3.pcs.mtrx[, c("speed") ] ,
# #     initRespDistance = initial.pcs.mtrx[, c("revDistance") ] ,
# #     habituationDistance = initial.pcs.mtrx[, c("revDistance") ] - final3.pcs.mtrx[, c("revDistance") ] ,
# #     spontSR = spontaneous.pcs.mtrx[, c("probSR") ],
# #     spontDuration = spontaneous.pcs.mtrx[, c("revDuration") ] ,
# #     spontSpeed = spontaneous.pcs.mtrx[, c("speed") ] ,
# #     spontDistance = spontaneous.pcs.mtrx[, c("revDistance") ] 
# #   )
#############################################################

usableFeatures <- c(
  c("pathlen",
    "bias",
    "dir",
    "angular",
    "aspect",
    "width",
    "length",
    "kink",
    "crab",
    "speed",
    "morphwidth",
    "curve",
    "area",
    "midline"
  ),
  colnames(reversals)
)

# Burn-in average features versus plus reversals
# Per variant
wormData <- cbind(burnIn.cm.mtrx, reversals)
wormData.sc <- scaleCols(wormData)
wormData.plate <- cbind(burnIn.pcm.mtrx, reversals.plate)
wormData.plate.sc <- scaleCols(wormData.plate)

## Per plate
wormData.mtrx <- as.matrix(wormData)
wormData.sc.mtrx <- as.matrix(wormData.sc)
wormData.plate.mtrx <- as.matrix(wormData.plate)
wormData.plate.sc.mtrx <- as.matrix(wormData.plate.sc)

# .df include gene/plate name as column
wormGenes <- rownames(burnIn.cm.mtrx)
wormData.sc.df <- cbind(wormGenes, wormData.sc.mtrx)
wormData.feature.df <- wormData[, usableFeatures]
wormData.feature.sc.df <- wormData.sc.df[, usableFeatures]

wormPlates <- rownames(burnIn.pcm.mtrx)
wormData.plate.sc.df <- cbind(wormPlates, wormData.plate.sc.mtrx)
wormData.plate.feature.df <- wormData.plate[, usableFeatures]
wormData.plate.feature.sc.df <- wormData.plate.sc.df[, usableFeatures]


setdiff(colnames(wormData.sc.mtrx), usableFeatures) # Shows dropped features
setdiff(usableFeatures, colnames(wormData.sc.mtrx)) # Should be empty

wormData.feature.mtrx <- wormData.mtrx[, usableFeatures]
wormData.feature.sc.mtrx <- wormData.sc.mtrx[, usableFeatures]
colnames(wormData.feature.mtrx)

wormData.plate.feature.mtrx <- wormData.plate.mtrx[, usableFeatures]
wormData.plate.feature.sc.mtrx <- wormData.plate.sc.mtrx[, usableFeatures]
colnames(wormData.plate.feature.mtrx)


#===========================
writeObject(after.cm)
writeObject(burnIn.cm)
writeObject(spontaneous.cm)
writeObject(initial.cm)
writeObject(final3.cm)
writeObject(second.cm)

# writeObject(burnIn.cm.cpmtrx)
# burnIn.feature.sc.cm.cpmtrx <- burnIn.sc.cm.cpmtrx[, usableFeatures[usableFeatures %in% colnames(burnIn.sc.cm.cpmtrx)] ]
# writeObject(burnIn.feature.sc.cm.cpmtrx)
#===========================
writeObject(after.cm.mtrx)
writeObject(burnIn.cm.mtrx)
writeObject(spontaneous.cm.mtrx)
writeObject(initial.cm.mtrx)
writeObject(final3.cm.mtrx)
writeObject(second.cm.mtrx)
writeObject(reversals)
writeObject(reversals.sd)
#===========================
writeObject(after.pcm.mtrx)
writeObject(burnIn.pcm.mtrx)
writeObject(spontaneous.pcm.mtrx)
writeObject(initial.pcm.mtrx)
writeObject(final3.pcm.mtrx)
writeObject(second.pcm.mtrx)
writeObject(reversals.plate)
#===========================

writeObject(innerDF)


writeObject(wormData)
writeObject(wormData.feature.df)
writeObject(wormData.sc.df)

writeObject(wormData.mtrx)
writeObject(wormData.sc.mtrx)

writeObject(wormData.feature.mtrx)
writeObject(wormData.feature.sc.mtrx)


writeObject(wormData.plate)
writeObject(wormData.plate.feature.df)
writeObject(wormData.plate.sc.df)

writeObject(wormData.plate.mtrx)
writeObject(wormData.plate.sc.mtrx)

writeObject(wormData.plate.feature.mtrx)
writeObject(wormData.plate.feature.sc.mtrx)