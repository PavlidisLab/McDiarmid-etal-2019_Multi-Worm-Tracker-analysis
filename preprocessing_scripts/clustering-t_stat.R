source("preprocessing_scripts/utils.R")

## Set clustering parameters
NCPU=config::get("clustering")$ncpu # Use NCPU parallel local compute nodes
NBOOT=config::get("clustering")$nboot # Bootstrap rounds
SEED=config::get("clustering")$seed # Seed for randomization
CONFIDENCE=config::get("clustering")$alpha # alpha cutoff for clustering
#CONFIDENCE=c(0.90,0.95, 0.99) # It's also possible to specify multiple thresholds.

MAX.ONLY=T # Only show maximum clusters

if (exists("COMPUTE")) {
  parallel::stopCluster(COMPUTE)
}
COMPUTE <- parallel::makePSOCKcluster(NCPU)

reversals.src <- readObjectAsMatrix("reversals")[, ]
t.stat <- readObjectAsMatrix("t.stat")
isWT <- grepl("^N2", rownames(t.stat))
t.stat <- t.stat[!isWT,] # Don't cluster on the WTs

wormData.feature.mtrx <- t.stat
wormData.feature.mtrx.reversals <-t.stat[,colnames(reversals.src)]
wormData.feature.mtrx.morpho <- t.stat[,setdiff(colnames(t.stat), 
                                                colnames(reversals.src))]

wormData.feature.mtrx.learning <-t.stat[,LEARNING_FEATURES] 

wormData.feature.mtrx.plate <- readObjectAsMatrix("wormData.plate.feature.mtrx")
days.melted <- getDateInformationByPlate(wormData.feature.mtrx.plate)

# Set up clusters
pv.genes.all <- pvclust(data = t(wormData.feature.mtrx), 
                        nboot = NBOOT,
                        iseed = SEED, 
                        parallel = COMPUTE
)

pv.genes.morpho <- pvclust(data = t(wormData.feature.mtrx.morpho), 
                           nboot = NBOOT,
                           iseed = SEED, 
                           parallel = COMPUTE
)

pv.genes.reversals <- pvclust(data = t(wormData.feature.mtrx.reversals), 
                              nboot = NBOOT,
                              iseed = SEED, 
                              parallel = COMPUTE
)

pv.genes.learning <- pvclust(data = t(wormData.feature.mtrx.learning), # TODO/FIXME: Is this correct? 
                             nboot = NBOOT,
                             iseed = SEED, 
                             parallel = COMPUTE 
)

# Produce plots
pvplot(pv.genes.all, breaks=c(CONFIDENCE), filename = "t.test.pvclust_all.png", HIGHLIGHT = days.melted) 
pvplot(pv.genes.morpho, breaks=c(CONFIDENCE), filename = "t.test.pvclust_morpho.png", HIGHLIGHT = days.melted) 
pvplot(pv.genes.reversals, breaks = c(CONFIDENCE), filename = "t.test.pvclust_reversals.png", HIGHLIGHT = days.melted) 
pvplot(pv.genes.learning, breaks = c(CONFIDENCE), filename = "t.test.pvclust_learning.png", HIGHLIGHT = days.melted) 

for (file in paste0(config::get("clustering")$plots, c("t.test.pvclust_all.png" , "t.test.pvclust_morpho.png" , "t.test.pvclust_reversals.png", "t.test.pvclust_learning.png")) ){
  system( paste("convert -rotate 90 ", file, file))
}

# seplot(pv.genes.all)
# seplot(pv.genes.morpho)
# seplot(pv.genes.reversals)
# seplot(pv.genes.learning)

reversal.clusters <- pvpick(x = pv.genes.reversals, 
                            alpha = CONFIDENCE, 
                            pv="au", 
                            type="geq", 
                            max.only = MAX.ONLY)$clusters

morpho.clusters <- pvpick(x = pv.genes.morpho, 
                          alpha = CONFIDENCE, 
                          pv="au", 
                          type="geq", 
                          max.only = MAX.ONLY)$clusters

all.clusters <- pvpick(x = pv.genes.all, 
                       alpha = CONFIDENCE, 
                       pv="au", 
                       type="geq", 
                       max.only = MAX.ONLY)$clusters

learning.clusters <- pvpick(x = pv.genes.learning, 
                            alpha = CONFIDENCE, 
                            pv="au", 
                            type="geq", 
                            max.only = MAX.ONLY)$clusters


reversal.clusters.df <- meltClusters(reversal.clusters)
morpho.clusters.df <- meltClusters(morpho.clusters)
all.clusters.df <- meltClusters(all.clusters)
learning.clusters.df <- meltClusters(learning.clusters)

reversal.clusters.df$group <- paste0("T.REV", reversal.clusters.df$group)
morpho.clusters.df$group <- paste0("T.MORPH", morpho.clusters.df$group)
all.clusters.df$group <- paste0("T.ALL", all.clusters.df$group)
learning.clusters.df$group <- paste0("T.LEARNING", learning.clusters.df$group)

writeObject(reversal.clusters.df,filename = "t.test.reversal.clusters.df")
writeObject(morpho.clusters.df, filename = "t.test.morpho.clusters.df" )
writeObject(all.clusters.df, filename = "t.test.all.clusters.df")
writeObject(learning.clusters.df, filename = "t.test.learning.clusters.df")

