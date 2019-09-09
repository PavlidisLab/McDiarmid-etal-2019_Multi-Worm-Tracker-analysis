t.stat <- readObjectAsMatrix("t.stat") 
t.pval <- readObjectAsMatrix("t.pval") 
t.qval <- readObjectAsMatrix("t.qval") 

t.direction.template <- t.stat
t.direction.template[T] <- 0
t.direction.template[t.stat > 0] <- 1
t.direction.template[t.stat < 0] <- -1

## Absolute
sum.abs.strain = rowSums(t.direction.template)
sum.abs.phenotype = colSums(t.direction.template)

## Up
t.direction.up = t.direction.template
t.direction.up[t.direction.up != 1] <- 0
sum.up.strain = rowSums(t.direction.up)
sum.up.phenotype = colSums(t.direction.up)

## Down
t.direction.down = t.direction.template
t.direction.down[t.direction.down != -1] <- 0
sum.down.strain = rowSums(abs(t.direction.down))
sum.down.phenotype = colSums(abs(t.direction.down))


## Change template to remove q-vaL > 0.1
## Absolute
t.direction.template[t.qval > 0.1] <- 0
sum.abs.strain.q01 = rowSums(t.direction.template)
sum.abs.phenotype.q01 = colSums(t.direction.template)

## Up
t.direction.up = t.direction.template
t.direction.up[t.direction.up != 1] <- 0
sum.up.strain.q01 = rowSums(t.direction.up)
sum.up.phenotype.q01 = colSums(t.direction.up)

## Down
t.direction.down = t.direction.template
t.direction.down[t.direction.down != -1] <- 0
sum.down.strain.q01 = rowSums(abs(t.direction.down))
sum.down.phenotype.q01 = colSums(abs(t.direction.down))

## Dataframes
t.direction.count.strain = data.frame(
  Delta = sum.abs.strain,
  Up = sum.up.strain,
  Down = sum.down.strain,
  `Delta_q_lt_0.1` = sum.abs.strain.q01,
  `Up_q_lt_0.1` = sum.up.strain.q01,
  `Down_q_lt_0.1` = sum.down.strain.q01
)


t.direction.count.phenotype =  data.frame(
  Delta = sum.abs.phenotype,
  Up = sum.up.phenotype,
  Down = sum.down.phenotype,
  `Delta_q_lt_0.1` = sum.abs.phenotype.q01,
  `Up_q_lt_0.1` = sum.up.phenotype.q01,
  `Down_q_lt_0.1` = sum.down.phenotype.q01
)

N2s = t.direction.count.strain[grepl(rownames(t.direction.count.strain), pattern = "^N2"),] 
t.direction.count.strain = t.direction.count.strain[!grepl(rownames(t.direction.count.strain), pattern = "^N2"),] 

rownames(t.direction.count.strain) <- translate(rownames(t.direction.count.strain), dictionary = "alleles")
rownames(t.direction.count.phenotype) <- translate(rownames(t.direction.count.phenotype), dictionary = "features")

t.direction.count.strain <- rbind(t.direction.count.strain, N2s) # Re-add N2s
# View(t.direction.count.strain)
# View(t.direction.count.phenotype)
writeObject(t.direction.count.strain)
writeObject(t.direction.count.phenotype)
