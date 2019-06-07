source("preprocessing_scripts/utils.R")

wormData.feature.mtrx <- readObjectAsMatrix("wormData.feature.mtrx")
wormData.feature.mtrx.plate <- readObjectAsMatrix("wormData.plate.feature.mtrx")

days.melted <- getDateInformationByPlate(wormData.feature.mtrx.plate)
days.WT <- grepl("^N2", days.melted$plate)
days.with.controls = unique(days.melted$day[ days.WT ])

row.names(wormData.feature.mtrx) <- cleanAlleleNames(wormData.feature.mtrx)
row.names(wormData.feature.mtrx.plate) <- cleanAlleleNames(wormData.feature.mtrx.plate)

plate.days <- gsub(pattern = ".*:", replacement = "", gsub(pattern = "_.*", replacement = "", row.names(wormData.feature.mtrx.plate)))

alleles <- clipPlate(cleanAlleleNames(wormData.feature.mtrx.plate))
wormData.feature.mtrx <- wormData.feature.mtrx[unique(alleles),]

wormData.feature.mtrx.plate.squashed <- simpleSquash(wormData.feature.mtrx.plate, 
                                                     groupVector = alleles)


if (!all(cleanAlleleNames(wormData.feature.mtrx) == cleanAlleleNames(wormData.feature.mtrx.plate.squashed)))  { stop("[ERROR] Allele name mismatch!") } 

df.isWT <- grepl("^N2", rownames(wormData.feature.mtrx.plate))

t.stat <- wormData.feature.mtrx.plate.squashed
t.stat[t.stat <= Inf ]  <- NA

t.pval <- wormData.feature.mtrx.plate.squashed
t.pval[t.pval <= Inf ]  <- NA

t.mean <- wormData.feature.mtrx.plate.squashed
t.mean[t.mean <= Inf ]  <- NA

CI95.lo <- wormData.feature.mtrx.plate.squashed
CI95.lo[CI95.lo <= Inf ]  <- NA

CI95.hi <- wormData.feature.mtrx.plate.squashed
CI95.hi[CI95.hi <= Inf ]  <- NA

for (feature in colnames(wormData.feature.mtrx.plate.squashed)){
  for (allele in cleanAlleleNames(wormData.feature.mtrx.plate.squashed) ){
    # Plates with allele of interest
    alleles.filter <- clipPlate(cleanAlleleNames(wormData.feature.mtrx.plate)) %in% allele
    
    # Days where allele of interest was run
    days.allele <- unique(days.melted[gsub(pattern=":.*", replacement = "", x=gsub(pattern=" ", replacement = "", x=days.melted$plate)) == allele, "day"])
    matched.WT.filter = df.isWT 
    
    grp1 <-   wormData.feature.mtrx.plate[alleles.filter,feature]
    grp2 <-   wormData.feature.mtrx.plate[matched.WT.filter,feature]
    
    # Number of controls
    N.Ctrl = length(unique(gsub(pattern = ":.*", replacement = "", x = names(grp2))))
    
    if (sum(matched.WT.filter) == 0){
      print(paste("No WT for",allele))
      next;  
    } else if ( length(days.allele) > N.Ctrl) { # If there's more "days" than controls, there's a problem.
      print(paste("Missing control for ",allele, ". Expected", length(days.allele), " and got ", N.Ctrl))
      next;  
    }
    
    t.test.results<-t.test(x = grp1, y = grp2)
    
    t.stat[allele, feature] <- t.test.results$statistic
    t.pval[allele, feature] <- t.test.results$p.value
    t.mean[allele, feature] <-  t.test.results$estimate[1] - t.test.results$estimate[2]
    CI95.lo[allele, feature] <- min(t.test.results$conf.int)
    CI95.hi[allele, feature] <- max(t.test.results$conf.int)
    
   
  }
}

## summary(t.stat) # Sanity check
## iheatmapper( t.stat ) # Visual sanity check

## ================================
## Save objects
writeObject(t.stat)
writeObject(t.pval)
writeObject(t.mean)
writeObject(CI95.lo)
writeObject(CI95.hi)