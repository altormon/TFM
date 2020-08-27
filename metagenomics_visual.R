## Pavian

# devtools::install_github('fbreitwieser/pavian') ## requires library devtools
pavian::runApp(port=5000)

## Rarefaction Vegan

# install.packages("vegan")
library(vegan); library(openxlsx); library(readxl); library(edgeR); library(ggplot2); library(NOISeq)
filename <- "data.xlsx"
bioenv <- openxlsx::read.xlsx(filename, rowNames = F, colNames = F, detectDates = F)
rnames <- bioenv[2:nrow(bioenv), 1]
cnames <- as.integer(bioenv[1, 2:ncol(bioenv)])
dates <- as.Date(cnames, origin = "1899-12-30")
cnames <- dates
bioenv <- bioenv[2:nrow(bioenv),2:ncol(bioenv)]
rownames(bioenv) <- rnames
colnames(bioenv) <- cnames
bioenv <- t(bioenv)

# rarefaction (relation between number of sequences in a sample and number of taxa expected based on the abundances detected)
S <- specnumber(bioenv) # observed number of species
(raremax <- min(rowSums(bioenv)))
Srare <- rarefy(bioenv, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(bioenv, step = 20, sample = raremax, col = "blue", cex = 0.6, xlab = "Cumulative DNA Sequences", ylab = "Cumulative Taxa Discovered")

# alpha diversity (estimation of richness (number) and distribution of taxa in a sample)
alfa <- fisher.alpha(bioenv, base = exp(1))
barplot(alfa, ylab = "Shannon Index (log)")
