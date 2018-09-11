# Rosenberg et al. 2018 
# "Single-cell profiling of the developing mouse brain and spinal cord with split-pool barcoding"
# [DOI: 10.1126/science.aam8999]

# supplementary table 5 contains the TPM + 1 values across all cell populations
# URL http://science.sciencemag.org/highwire/filestream/707022/field_highwire_adjunct_files/4/aam8999_TableS5.xlsx

# load libraries
# have to be installed with install.packages() or BiocInstaller::biocLite()
library(reshape2)
library(data.table)
library(ggplot2)
library(readxl)

# path to Supp5 table
# read XLSX and convert to data.table
Supp5 <- "Data/Rosenberger_et_al_2018/aam8999_TableS5.xlsx"
TPMs <- readxl::read_excel(Supp5)
colnames(TPMs) <- TPMs[1,]
TPMs <- TPMs[-1,]
colnames(TPMs)[1] <- "GeneSymbol"
TPMs <- data.table::as.data.table(TPMs)

# get list of histone genes
histoneGenes <- data.table::fread("Data/histoneGenes.txt")

# get histone specific expression values [TPM +1 = transcript per million + 1 offset to avoid problems when log2 transforming]
# first try to match gene names from histone table to gene names in TPM table
# unfortunately only about half of the genes are matched
# have to use better mapping at later stage
histoneTPMs <- TPMs[which(toupper(TPMs$GeneSymbol) %in% histoneGenes$ApprovedSymbol)]
histoneTPMs <- rehshape2::melt(histoneTPMs, id.vars = 1, variable.name = "CellType")
histoneTPMs <- data.table::data.table(histoneTPMs)
data.table::setkey(histoneTPMs, GeneSymbol)
histoneTPMs$value <- as.double(histoneTPMs$value)

# reordering the factor levels according to the integer in the cell type name
histoneTPMs$CellType <- factor(as.character(histoneTPMs$CellType),
       levels(histoneTPMs$CellType)[order(as.integer(unlist(lapply(strsplit(as.character(histoneTPMs$CellType[1:73]), " "), function(x) x[1]))))])

# preparing colours for manual color assignment to use discriminative scheme
# generated using http://tools.medialab.sciences-po.fr/iwanthue/
n <- length(unique(histoneTPMs$GeneSymbol))
palette <- c("#5eb9d9", "#e45228", "#55e243", "#a33ad4", "#9eed34", "#5b41d0", "#e5e82f", "#412487", "#84eb65", "#e04acf", "#46b439",
             "#a42d90", "#b9e34e", "#586dd9", "#d6d44b", "#af79e4", "#83b537", "#7f479e", "#5ce491", "#d83879", "#4b9f52", "#dd66b6",
             "#c0e37f", "#322155", "#eab937", "#5768ab", "#ac9f2f", "#68a0e4", "#d4872d", "#63dcc2", "#dc3d4b", "#95db97", "#87295d",
             "#7e963c", "#df9ad9", "#3a6019", "#b1a5dd", "#a94421", "#a1dbde", "#81272b", "#cde3b5", "#3f1b23", "#e2cf7c", "#5e3e5b",
             "#d3a76f", "#243f49", "#d7805f", "#3c7f6a", "#d76a7d", "#3a703f", "#a36d9c", "#857023", "#5e7493", "#7b4c28", "#c8bed5",
             "#32391e", "#de9da0", "#609a9c", "#925b5e", "#8fab7f", "#998474", "#dac3ac", "#6c6941")
names(palette) <- unique(histoneTPMs$GeneSymbol)
palette <- data.table(GeneSymbol = names(palette), color = palette)

# actual plotting
# all mapped histones
ggplot(histoneTPMs, aes(CellType, log2(value), group = GeneSymbol)) +
  geom_line(aes(color = GeneSymbol)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_color_manual(values=palette$color)

# two candidates only
ggplot(histoneTPMs[c("H2afb3", "H2afz")], aes(CellType, log2(value), group = GeneSymbol)) +
  geom_line(aes(color = GeneSymbol)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_color_manual(values=palette$color)
       
