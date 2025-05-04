library(phyloseq)
library(tidyr)
library(dplyr)
library(ggplot2)


#Make the phyloseq object
map_file <- read.csv("~/Desktop/Xesto2025/xesto_reads/xesto454_metadata_new.csv", row.names=NULL, stringsAsFactors=TRUE)
# Specify which column contains the sample names (assuming it's named 'SampleID' in the map file)
rownames(map_file) <- map_file$SampleID  # Replace 'SampleID' with the actual column name in your map file
#map_file$SampleID <- NULL  # Optionally remove the 'SampleID' column if it's no longer needed
sample_data_ps <- sample_data(map_file)  # This now uses the 'SampleID' column as the row names



library(readr)
asv <- read_csv("asv.csv")


data <- read_csv("data.csv")
data <- as.data.frame(data)
asv <- as.data.frame(asv)
asv[] <- lapply(asv, as.numeric)
asv[is.na(asv)] <- 0
asv_ps <- otu_table(asv, taxa_are_rows = TRUE)

tax_ps <- tax_table(as.matrix(data))

ps <- phyloseq(tax_ps,asv_ps)
ps
ps1 = merge_phyloseq(ps,sample_data_ps)
ps1


#CSS normalization
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("metagenomeSeq")
a
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

library(metagenomeSeq)

sum(taxa_sums(ps1) == 0) #Calculate current number of zeros in data
ps1 <- filter_taxa(ps1, function(x) sum(x) != 0, TRUE) #Filters zeros from data
sum(taxa_sums(ps1) == 0) #Recalculates number of zeros






sort(sample_sums(ps1))
ps1s = subset_samples(ps1, SampleID != "5_11_Soil") #Filtering out to fix counts
ps1s = subset_samples(ps1s, SampleID != "5_29_Soil") #Filtering out to fix counts
ps1s = subset_samples(ps1s, SampleID != "5_11_Water") #Filtering out to fix counts
ps1s
sort(sample_sums(ps1s))
MX <- phyloseq_to_metagenomeSeq(ps1s)
MX

p <- cumNormStatFast(MX) #Normalization (Cumulative Sums)
p
MX<-cumNorm(MX,p=p)

normFactors(MX) #Returns normalization factors
normmybiom<-MRcounts(MX,norm=T)

exportMat(normmybiom, file = "Xest_CSS_norm.txt", sep = "\t")
b <- MRexperiment2biom(MX, norm = T)

library(biomformat)
write_biom(b, biom_file = "Xesto_CSS_norm.biom")
library(RCurl)
import_biom2_script <- getURL("https://gist.githubusercontent.com/jnpaulson/324ac1fa3eab1bc7f845/raw/2ef62334d4e9bc5446a5ee6dd198f52484097dae/import_biom2.R", ssl.verifypeer = FALSE)
eval(parse(text = import_biom2_script))

biom <- read_biom("Xesto_CSS_norm.biom")
biom2 <- import_biom2(biom)
sort(sample_sums(biom2))

sort(sample_sums(ps1s))
norm_ps<-merge_phyloseq(biom2,sample_data_ps,tax_ps)
norm_ps

#Here, we will calculate the gap statistic (need for eventual k-means evaluation). Below is an add-on script
library("cluster")
theme_set(theme_bw())
exord = ordinate(ps1, method="MDS", distance="jsd")
pam1 = function(x, k){list(cluster = pam(x,k, cluster.only=TRUE))}
x = phyloseq:::scores.pcoa(exord, display="sites")
gskmn = clusGap(x[, 1:2], FUN=kmeans, nstart=20, K.max = 6, B = 500)
gskmn

gap_statistic_ordination = function(exord, FUNcluster, type="DiseaseStates", K.max=6, axes=c(1:2), B=500, verbose=interactive(), ...){
  require("cluster")
  #   If "pam1" was chosen, use this internally defined call to pam
  if(FUNcluster == "pam1"){
    FUNcluster = function(x,k) list(cluster = pam(x, k, cluster.only=TRUE))     
  }
  # Use the scores function to get the ordination coordinates
  x = phyloseq:::scores.pcoa(exord, display=type)
  #   If axes not explicitly defined (NULL), then use all of them
  if(is.null(axes)){axes = 1:ncol(x)}
  #   Finally, perform, and return, the gap statistic calculation using cluster::clusGap  
  clusGap(x[, axes], FUN=FUNcluster, K.max=K.max, B=B, verbose=verbose, ...)
}

#Now we can plot it.
plot_clusgap = function(gap_statistic_ordination, title="Gap Statistic calculation results"){
  require("ggplot2")
  gstab = data.frame(gap_statistic_ordination$Tab, k=1:nrow(gap_statistic_ordination$Tab))
  p = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size=5)
  p = p + geom_errorbar(aes(ymax=gap+SE.sim, ymin=gap-SE.sim))
  p = p + ggtitle(title)
  return(p)
}

gs = gap_statistic_ordination(exord, "pam1", B=50, verbose=FALSE)



phylum_colors <- c( "#E69F00", "#56B4E9",  "#0072B2", "#D55E00","#009E73","#F0E442","#999999", "#000000","#CC79A7","#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                    "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
                    "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")

colors <-c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", 
  "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", 
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FF0000", "#F1F1F1", "#C2B0D0", 
  "#7D92A2", "#AB5C5B", "#9BCE61", "#F8A1C0", "#FEFF00", "#1F78B4", "#33A02C", "#FF00FF", 
  "#FF4C3B", "#4E79A7", "#F4A82A", "#D8BFD8", "#E6B3B1", "#A9D08E", "#FF6F91", "#D2D2D2", 
  "#ADFF2F", "#BFD8D2", "#F5A623", "#8B0000", "#A52A2A", "#0D98BA", "#8C9E91", "#D50000", 
  "#4A90E2", "#11C0FF", "#0A7D8C", "#B9E5FF", "#48C9B0", "#79C7B2", "#A3E4D7", "#7F5C81", 
  "#2E7D32", "#003C73", "#F1C6E7", "#1D96B0", "#C1E6FC","#599861", "#CC5A53", "#E5E1E6", "#5A3C5E", 
  "#FAE5D3", "#5D4080", "#003366", "#CA3E8F", "#A0E0A4", "#F2BB77", "#3F51B5", "#D6DAFF")
ps1_gg <- ps1_rel %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at phylum level
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Genus)                                      # Sort data frame alphabetically by phylum

barplot<-ggplot(ps1_gg, aes(x ="DiseaseState", y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(angle=270, hjust=1,size=18))+
  theme(axis.text.y = element_text(size=24))+
  theme(legend.text=element_text(size=18),legend.title=element_text(size=24))+
  theme(strip.text.x = element_text(size = 26))+
  theme(axis.text=element_text(size=24),axis.title=element_text(size=24))+
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Abundance \n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



