#############################################################
# Megan Brudz-Rodríguez
# GSE53798 - CSC 450
# Summary from GEO:
# We used global gene expression profiles from human B-cell cell
# lines to generate gene expression signatures for prediction
# of response to the drugs cyclophosphamide, doxorubicin or
# vincristine. The signatures were validated in two publicly
# available clinical cohorts.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse53798
#############################################################

# 1 libraries
library(GEOquery)
library(limma)


# 2 download process data
GSE53798 = getGEO("GSE53798")
GSE53798.expr = exprs(GSE53798[[1]])
GSE53798.p = pData(GSE53798[[1]])


# 3 boxplot
boxplot(GSE53798.expr, main = "processed data", xlab = "Sample",
        ylab = "Expression")

# no need to normalize data

# 4 number of samples and probes
# number of samples
ncol(GSE53798.expr)

# number of probes
nrow(GSE53798.expr)


# 5 comparing 
cellfrac = as.character(GSE53798.p$characteristics_ch1.6)

# separate "doxorubicin" from Resistant and Sensitive found via "View(GSE53798.p)"
sample.names = as.character(cellfrac)
s = strsplit(sample.names, ": ")
get.last <- function(x) {return(x[[2]])}  
cellfrac = sapply(s, get.last)

# amount in control and treatment group
sum(cellfrac == "Resistant")
sum(cellfrac == "Sensitive")


# 6 differentially expressed probes
groups = GSE53798.p$characteristics_ch1.6
levels(groups) = c("Intermediate", "Resistant", "Sensitive")

# remove "Intermediate"
index = groups != "Intermediate"
groups = groups[index]
GSE53798.expr = GSE53798.expr[,index]

groups = as.character(groups)

design = model.matrix(~0+groups)
head(design)

colnames(design) = c("resistant", "sensitive")
fit = lmFit(GSE53798.expr, design)
head(fit$coefficients)
head(fit$sigma)

contrast.matrix <- makeContrasts(resistant - sensitive,levels=design)

fit = contrasts.fit(fit, contrast.matrix)
head(fit$coefficients)
head(fit$sigma)

fit = eBayes(fit)

tt = topTable(fit,sort.by = "p")

tt.100 = topTable(fit,sort.by = "p", number = 100)
nrow(tt.100)

# 7 top probes
tt.50 = topTable(fit, number = 50)
set2 = tt.50
head(tt.50) # top probes
# top 3 probes w/ logFC & adj.p-val
tt.50 = topTable(fit, number = 3)
tt.50 = head(tt.50) # top 50 probes
cols = c(1,5)
tt.50[1:3, cols]
adj = tt.50$adj.P.Val
adj[3] # value of FDR using final adj.P.Value

# 8 probe with lowest adjusted p-val: 
min(tt.50$adj.P.Val)
logFC = tt.50[1,]$logFC
FC = 2**logFC
boxplot(tt.50[1,], main = "Fold Change = 0.1361205", xlab = "Probe 203349_s_at",
        ylab = "Expression")

# 9 find probe that corresponds to specific gene, get GPL
annotation(GSE53798[[1]])
# get platform for this dataset (the 1st dataset in the list)
platform = annotation(GSE53798[[1]])

# 10 download GPL
GPL = getGEO(platform)

# 11 top genes corresponding to all probes
genes.db = Table(GPL)
m = match(rownames(tt.100), genes.db$ID)
genes = as.character(genes.db$'Gene Symbol'[m])[]
keep = genes!=""
genes = genes[keep]
genes = strsplit(genes, " /// ")
get.first <-function(x) {
  x[1]
}
genes = sapply(genes, get.first)
genes = unique(genes)
write.table(genes[2:4], row.names = FALSE, quote = FALSE)

# 12 top 2 genes according to gene cards
# HEG1 and LGALS8


# 13 using DAVID for Gene Ontology (GO) terms of KEGG
# HEG1 is called HEG homolog 1 (zebrafish.  It has a role
# in regulating concentric heart growth in zebrafish.  It
# is also a protein coding gene.

# LGALS8 is called lectin, galactoside-binding, soluble, 8.
# It has been used in  development, differentiation, cell-cell
# adhesion, cell-matrix interaction, growth regulation,
# apoptosis, and RNA splicing.  It is widely expressed in
# tumoral tissues and seems to be involved in integrin-like
# cell interactions. Alternatively spliced transcript variants
# encoding different isoforms have been identified.


# 14 heatmap
# ten probes are from tt
probes2 = rownames(tt)


col.heat = colorRampPalette(c("yellow", "blue"))(200)
m = match(probes2, rownames(GSE53798.expr))

PL = getGEO("GPL570")
PL = Table(PL)

m.genes = match(probes2, PL$ID)
genes = PL$`Gene Symbol`[m.genes]
genes = paste0(genes, " (", probes2, ")")

rownames(GSE53798.expr)[m] = genes

col.response = as.integer(factor(groups))
col.response = c("darkred","darkgreen")[col.response]
table(col.response)

heatmap(GSE53798.expr[m,], col=col.heat, ColSideColors = col.response)

# venn diagram
area2 = probes2

