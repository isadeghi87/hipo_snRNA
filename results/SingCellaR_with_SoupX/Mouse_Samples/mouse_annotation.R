# Set working directory
setwd("/Users/i439h/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/hipo_temp/results/SingCellaR_with_SoupX/Mouse_Samples/")

# Load the processed data
mouse <- readRDS("./all_mouse_samples.SoupX.SingCellaR.Harmony.v1.rds")


# Define cell type mapping from clusters
cell_type_mapping <- c(
  "cl1"  = "OPC",
  "cl2"  = "G2",
  "cl3"  = "MES2_like",
  "cl4"  = "OPC_like",
  "cl5"  = "Neuron",
  "cl6"  = "OPC",
  "cl7"  = "OPC",
  "cl8"  = "Microglia",
  "cl9"  = "OC",
  "cl10" = "NPC",
  "cl11" = "NPC1_like",
  "cl12" = "NFOC",
  "cl13" = "Astrocyte",
  "cl14" = "Macrophage",
  "cl15" = "NPC",
  "cl16" = "Microglia",
  "cl17" = "Neuron",
  "cl18" = "OC",
  "cl19" = "Astrocyte",
  "cl20" = "Neuron",
  "cl21" = "NPC",
  "cl22" = "Pericytes",
  "cl23" = "Neuron",
  "cl24" = "Neuron",
  "cl25" = "Neuron",
  "cl26" = "Endothelial",
  "cl27" = "Neuron"
)

# Assign cell types to `mouse@sc.clusters`
mouse@sc.clusters$Cell_Type <- cell_type_mapping[as.character(mouse@sc.clusters$louvain_cluster)]

# Check if cell types are assigned correctly
table(mouse@sc.clusters$Cell_Type)
head(mouse@sc.clusters)


# Define Tumor vs. Normal Conditions
condition_mapping <- c(
  "cl1"  = "Tumor", "cl2"  = "Tumor", "cl3"  = "Tumor", "cl4"  = "Tumor",
  "cl5"  = "Normal", "cl6"  = "Tumor", "cl7"  = "Tumor", "cl8"  = "Tumor",
  "cl9"  = "Tumor", "cl10" = "Tumor", "cl11" = "Tumor", "cl12" = "Tumor",
  "cl13" = "Tumor", "cl14" = "Tumor", "cl15" = "Tumor", "cl16" = "Tumor",
  "cl17" = "Normal", "cl18" = "Normal", "cl19" = "Normal", "cl20" = "Normal",
  "cl21" = "Tumor", "cl22" = "Normal", "cl23" = "Normal", "cl24" = "Normal",
  "cl25" = "Normal", "cl26" = "Normal", "cl27" = "Normal"
)

# Assign condition to `mouse@sc.clusters`
mouse@meta.data$Condition <- condition_mapping[as.character(mouse@sc.clusters$louvain_cluster)]

# Check the assignment
table(mouse@meta.data$Condition)
head(mouse@meta.data)

# Count cells per condition
condition_counts <- table(mouse@sc.clusters$Condition)
condition_proportions <- prop.table(condition_counts)

# Convert to DataFrame
condition_df <- as.data.frame(condition_proportions)
colnames(condition_df) <- c("Condition", "Proportion")

# Save results
write.csv(condition_df, "tumor_vs_normal_proportions.csv", row.names = FALSE)

# Print proportions
print(condition_df)

# Visualize Tumor vs. Normal Cell Distribution
#(a) Bar Plot of Cell Proportions
ggplot(condition_df, aes(x = Condition, y = Proportion, fill = Condition)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ggtitle("Proportion of Tumor vs Normal Cells") +
  xlab("Condition") + ylab("Proportion") +
  scale_fill_manual(values = c("Tumor" = "red", "Normal" = "blue"))

ggsave("tumor_vs_normal_distribution.pdf", width = 6, height = 4)

# Compute cell type proportions for Tumor vs Normal
cell_type_composition <- table(mouse@sc.clusters$Cell_Type, mouse@sc.clusters$Condition)
cell_type_proportions <- prop.table(cell_type_composition, margin = 2)

# Convert to DataFrame
cell_type_df <- as.data.frame(as.table(cell_type_proportions))
colnames(cell_type_df) <- c("Cell_Type", "Condition", "Proportion")
cell_type_df = subset(cell_type_df,Proportion != 0) 

# Save results
write.csv(cell_type_df, "cell_type_composition_tumor_vs_normal.csv", row.names = FALSE)

# Plot
ggplot(cell_type_df, aes(x = Condition, y = Proportion, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "stack",width = 0.5,color ='black',show.legend = F) +
  geom_text(aes(label = Cell_Type), 
            position = position_stack(vjust = 0.5),  # Center labels in bars
            size = 4, color = "white", fontface = "bold") +  # Adjust label size and color
  theme_minimal() +
  ggtitle("Cell Type Composition in Tumor vs Normal") +
  xlab("Condition") + ylab("Proportion")


ggsave("cell_type_composition_tumor_vs_normal.pdf", width = 10, height = 6)

saveRDS(mouse,"./all_mouse_samples.SoupX.SingCellaR.Harmony.annotated.Rds")
