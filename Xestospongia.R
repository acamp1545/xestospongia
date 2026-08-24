if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

library(phyloseq)
library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)
install.packages("compositions")
library(compositions)
library(devtools)
#Install qiime2r with github source code tarfile.
library(qiime2R)

setwd("/Users/alexandracampbell/Desktop/Xesto2025/xesto_reads")

list.files()
ASVs<-read_qza("denoise/table.qza")
names(ASVs)
ASVs$data[1:5,1:5]
ASVs$type
ASVs<-ASVs$data
colnames(ASVs) <- gsub("\\.fastq.*$", "", colnames(ASVs))

metadata <- read.csv("xesto454_metadata_new.csv", row.names = "SampleID", stringsAsFactors = FALSE)
head(row.names(metadata))

taxonomy<-read_qza("xesto_454_taxonomy.qza") 
head(taxonomy$data)
taxonomy<-parse_taxonomy(taxonomy$data)
head(taxonomy)

root_tree<-read_qza("rooted-tree.qza")
tree<- root_tree$data
plot(tree)

OTU <- otu_table(as.matrix(ASVs), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(taxonomy))
MAP <- sample_data(metadata)

sample_names(OTU)[1:5]
sample_names(MAP)[1:5]



physeq_object <- phyloseq(OTU, TAX, MAP,tree)

# Handle zeros before CLR transformation (e.g., add a pseudocount)
# A common approach is to add a small pseudocount (e.g., 1) to all counts
physeq_obj_clr <- transform_sample_counts(physeq_object, function(x) x + 1)

# Perform CLR transformation
# Note: clr() function from 'compositions' package expects samples in rows
otu_clr <- t(otu_table(physeq_obj_clr)) # Transpose to have samples in rows
otu_clr_transformed <- clr(otu_clr)

# Convert back to phyloseq otu_table format (taxa in rows)
OTU_clr <- otu_table(t(otu_clr_transformed), taxa_are_rows = TRUE)
physeq_obj_clr_final <- phyloseq(OTU_clr, TAX, MAP)



dist_matrix <- phyloseq::distance(physeq_obj_clr_final, method = "euclidean")

# --- 3. Constrained Analysis of Principal Coordinates (CAP) Ordination ---

# Define the constraining variable from your map file (e.g., 'Group')
# Replace 'Group' with the actual column name in your map file
constraining_variable <- "DiseaseState" 

cap_ordination <- capscale(dist_matrix ~ DiseaseState,data = metadata)

#plot setup
site_scores <- scores(cap_ordination, display = "sites")
site_df <- as.data.frame(site_scores)
site_df$SampleID <- rownames(site_df)
site_df <- cbind(site_df, metadata[rownames(site_df), , drop = FALSE])

# Plot
ggplot(site_df, aes(x = CAP1, y = CAP2, color = DiseaseState)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "CAP Ordination", x = "CAP1", y = "CAP2") +
  theme(plot.title = element_text(hjust = 0.5))

#SIMPER Analysis
otu_simper <- t(otu_table(physeq_obj_clr_final))
grouping_factor <- sample_data(physeq_obj_clr_final)[[constraining_variable]]

simper_result <- simper(otu_simper, group = grouping_factor)

# View summary of SIMPER results for specific comparisons (e.g., Group1 vs Group2)
# Replace 'Group1' and 'Group2' with your actual group names
# summary(simper_result, ordered = TRUE) # Shows all pairwise comparisons
summary(simper_result, ordered = TRUE, digits = 3, rev.ord = TRUE) # Example for specific comparison
names(simper_result)

###FUNCTION###
# A function to extract and prepare data for plotting for a single pair:
prepare_simper_data <- function(simper_pair, pair_name) {
  # The species names are stored in the row names
  data_frame <- as.data.frame(simper_pair)
  data_frame$species <- rownames(data_frame)
  # Calculate cumulative percentage explicitly if needed for sorting in plot
  data_frame <- data_frame %>%
    arrange(desc(average)) %>% # Sort by contribution
    mutate(cumul_contrib = cumsum(average / sum(average) * 100))
  data_frame$pair <- pair_name
  return(data_frame)
}



# Create directory for SIMPER plots (no warnings if it already exists)
dir.create("simper_plots", showWarnings = FALSE)

# Extract the list of comparisons from your SIMPER result
comparisons <- names(simper_result)


# Loop through each comparison
# Create directory for SIMPER plots
dir.create("simper_plots", showWarnings = FALSE)

comparisons <- names(simper_result)

for (comp in comparisons) {
  comp_data <- summary(simper_result)[[comp]]
  
  # Filter to species contributing to first 90% cumulative contribution
  comp_data_filtered <- comp_data[comp_data$cumsum < 0.90, ]
  
  # Skip if no species meet criteria
  if(nrow(comp_data_filtered) == 0){
    message(paste("Skipping comparison", comp, "- no species with cumsum < 0.90"))
    next
  }
  
  # Contribution column
  contrib_col <- "average"
  
  # Convert contribution to numeric (just in case)
  comp_data_filtered[[contrib_col]] <- as.numeric(as.character(comp_data_filtered[[contrib_col]]))
  
  # Reorder by decreasing contribution
  ord <- order(comp_data_filtered[[contrib_col]], decreasing = TRUE)
  comp_data_filtered <- comp_data_filtered[ord, ]
  
  species_names <- rownames(comp_data_filtered)
  filename <- paste0("simper_plots/plot_", comp, ".png")
  
  # Open PNG device
  png(filename, width = 800, height = 600)
  par(mar = c(10, 4, 4, 2) + 0.1)
  
  # Create bar plot
  barplot(comp_data_filtered[[contrib_col]],
          names.arg = species_names,
          las = 2,
          ylab = "Contribution (average)",
          main = paste("Species Contribution for Comparison:", comp),
          sub = "Species ordered by decreasing contribution")
  
  # Optional mean line
  abline(h = mean(comp_data_filtered[[contrib_col]]), col = "red", lty = 2)
  
  # Close device
  dev.off()
}

print("Plots saved in the 'simper_plots' directory")


comp_data <- summary(simper_result)[[comparisons[1]]]
str(comp_data)
head(comp_data)

#Richness plots of target taxa
install.packages("ggpubr")
library(ggpubr)

physeq_rel  = transform_sample_counts(physeq_object, function(x) x / sum(x) )
physeq_rel_df<-psmelt(physeq_rel)

blasto<-subset_taxa(physeq_rel, Genus=="Blastopirellula")
rhodo<-subset_taxa(physeq_rel, Order=="Rhodobacterales")
alpha<-subset_taxa(physeq_rel, Class=="Alphaproteobacteria")
planc<-subset_taxa(physeq_rel, Class=="Planctomycetes")
pseudoalt<-subset_taxa(physeq_rel, Genus=="Pseudoalteromonas")
gamma<-subset_taxa(physeq_rel, Class=="Gammaproteobacteria")

custom.colors <- function(n) {
  palette <- c("dodgerblue1", "skyblue4", "chocolate1", "seagreen4",
               "bisque3", "red4", "purple4", "mediumpurple3",
               "maroon", "dodgerblue4", "skyblue2", "darkcyan",
               "darkslategray3", "lightgreen", "bisque",
               "palevioletred1", "black", "gray79", "lightsalmon4",
               "darkgoldenrod1")
  if (n > length(palette))
    warning('palette has duplicated colours')
  rep(palette, length.out=n)
}

options(ggplot2.discrete.colour = c("red", "#af01ef"))
options(ggplot2.discrete.colour= list(c("red", "#af01ef"), custom.colors(99)))

# Source - https://stackoverflow.com/a
# Posted by jan-glx, modified by community. See post 'Timeline' for change history
# Retrieved 2025-11-11, License - CC BY-SA 4.0

# `custom.colors` defined as in question
scale_custom <- function(aesthetics, scale_name= "custom", ..., palette = custom.colors) discrete_scale(aesthetics = aesthetics, scale_name, ..., palette = palette)
scale_colour_custom <- function(...) scale_custom("colour", ...)
options(ggplot2.discrete.colour = scale_colour_custom)



title="Relative Abundance of Blastopirellula Across Disease States"
bar_blasto<-plot_bar(blasto, "Description","Abundance", title=title)
bar_blasto
ggsave(bar_blasto, file="bar_blasto.png")

title2="Relative Abundance of Rhodobacterales Across Disease States"
bar_rhodo<-plot_bar(rhodo, "Description","Abundance", "Genus", title=title2) + scale_fill_viridis_d()
bar_rhodo
ggsave(bar_rhodo, file="bar_rhodo.png")

title3="Relative Abundance of Alphaproteobacteria Across Disease States"
bar_alpha<-plot_bar(alpha, "Description","Abundance", "Family", title=title3) + scale_fill_viridis_d(direction=-1)
bar_alpha
ggsave(bar_alpha, file="alphaprot.png")

title4="Relative Abundance of Planctomycetaceae Across Disease States"
bar_planc<-plot_bar(planc, "Description","Abundance", "Genus", title=title4) + scale_fill_viridis_d()
bar_planc
ggsave(bar_planc, file="planctomycet.png")

title5="Relative Abundance of Pseudoalternomonas Across Disease States"
bar_pseudoalt<-plot_bar(pseudoalt, "Description","Abundance", "Genus", title=title5) + scale_fill_viridis_d()
bar_pseudoalt
ggsave(bar_pseudoalt, file="pseudoalt.png")

title6="Relative Abundance of Gammaproteobacteria Across Disease States"
bar_gamma<-plot_bar(gamma, "Description","Abundance", "Genus", title=title6) + scale_fill_viridis_d()
bar_gamma
ggsave(bar_gamma, file="gamma.png")

kruskal_results<-list()

Blastopirellula <-c("Blastopirellula")

for (taxon in Blastopirellula) {
  formulat_str <-paste(taxon, "~DiseaseState")
  kruskal_test_result <- kruskal.test(as.formula(formula_str), data = physeq_rel_df)
  kruskal_results[[taxon]] <-kruskal_test_result$p.value
}
  
}