# Ryan Gallagher
# 03.12.2025
#
# Check KO experiment efficiency
#
#

library(tidyverse)

norm.cts = read.csv("./data/normalized_counts.csv", row.names = 1)
raw.cts = read.csv("./data/counts.csv", row.names = 1)

meta = read.csv("./data/meta.csv")

get_abundances_for_gene_on_plate = function(gene.f, plate.f) {
  print(paste("Processing gene:", gene.f, "on plate:", plate.f))
  
  relevant_samples = meta %>% filter(gene == gene.f, plate == plate.f)
  treatment_samples = relevant_samples %>% filter(condition == "treatment" & startsWith(sample, gene))
  control_samples = meta %>% filter(condition != "treatment", plate == plate.f)
  
  if (nrow(treatment_samples) == 0 | nrow(control_samples) == 0) {
    print("Warning: No treatment or control samples found!")
    return(data.frame())  # Return an empty dataframe safely
  }
  
  treatment_abundance = raw.cts[row.names(raw.cts) == gene.f, treatment_samples$sample, drop = FALSE]
  control_abundance = raw.cts[row.names(raw.cts) == gene.f, control_samples$sample, drop = FALSE]
  
  if (ncol(treatment_abundance) < 3 | ncol(control_abundance) < 6) {
    print("Warning: Insufficient sample columns for renaming!")
    print(treatment_abundance)
    print(control_abundance)
    return(data.frame())  # Return an empty dataframe safely
  }
  
  treatment_abundance = treatment_abundance %>%
    dplyr::rename(trt1 = 1, trt2 = 2, trt3 = 3)
  control_abundance = control_abundance %>%
    dplyr::rename(trac1 = 1, trac2 = 2, trac3 = 3, wt1 = 4, wt2 = 5, wt3 = 6)
  
  # Compute averages
  treatment_avg = rowMeans(treatment_abundance, na.rm = TRUE)
  control_avg = rowMeans(control_abundance, na.rm = TRUE)
  
  # Add new columns
  result_df = cbind(treatment_abundance, control_abundance)
  result_df = result_df %>% 
    mutate(treatment_avg = treatment_avg, control_avg = control_avg)
  
  return(result_df)
}

baby.df = meta %>% select(gene, plate, condition) %>% distinct() %>% filter(gene != "BACH1", gene != "RAP1GDS1", gene != "SH3BP5") %>%
  filter(condition == "treatment") %>% select(-condition)

results_df = baby.df %>%
  pmap_dfr(~ get_abundances_for_gene_on_plate(..1, ..2) %>%
             mutate(gene = ..1, plate = ..2))

# Reshape data and compute % difference
plot_data = results_df %>%
  select(gene, treatment_avg, control_avg) %>%
  mutate(percent_diff = ((treatment_avg - control_avg) / control_avg) * 100) %>%  # Compute % Diff
  pivot_longer(cols = c(treatment_avg, control_avg), 
               names_to = "condition", 
               values_to = "average_abundance") %>%
  mutate(condition = recode(condition, 
                            "treatment_avg" = "Treatment", 
                            "control_avg" = "Control"))

# Compute max abundance per gene to determine sorting and grouping
gene_order = plot_data %>%
  group_by(gene) %>%
  summarise(max_abundance = max(average_abundance, na.rm = TRUE)) %>%
  arrange(desc(max_abundance)) %>%  # Sort from high to low
  mutate(rank = row_number())        # Assign rank for grouping

# Merge sorted rank into the main data
plot_data = plot_data %>%
  left_join(gene_order, by = "gene") %>%
  mutate(group = ntile(rank, 4))  # Split into 4 groups based on rank

# Create annotation text for % Diff (only for Treatment bars)
plot_data = plot_data %>%
  mutate(label = ifelse(condition == "Treatment", paste0(round(percent_diff, 1), "%"), NA))

# Store plots in a list
plot_list = list()

# Loop over each group and create separate sorted plots
for (g in unique(plot_data$group)) {
  plot_subset = plot_data %>% 
    filter(group == g) %>% 
    arrange(desc(max_abundance))  # Sort genes within each group
  
  p = ggplot(plot_subset, aes(x = reorder(gene, -max_abundance), 
                               y = average_abundance, fill = condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label), 
              position = position_dodge(width = 0.9), 
              vjust = -0.5, size = 4, na.rm = TRUE) +  # Annotate % Diff above bars
    labs(title = paste("Gene Expression Abundances - Group", g),
         x = "Gene",
         y = "Average Abundance") +
    scale_fill_manual(values = c("Treatment" = "blue", "Control" = "red")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Store plot in list
  plot_list[[paste0("Group_", g)]] = p
}

# Print individual plots
print(plot_list[["Group_1"]])  # View plot for Group 1
print(plot_list[["Group_2"]])  # View plot for Group 2
print(plot_list[["Group_3"]])  # View plot for Group 3
print(plot_list[["Group_4"]])  # View plot for Group 4


## Heatmap showing total data

# Compute % Difference and categorize into groups
heatmap_data = results_df %>%
  select(gene, treatment_avg, control_avg) %>%
  mutate(percent_diff = ((treatment_avg - control_avg) / control_avg) * 100) %>%
  mutate(group = ifelse(percent_diff >= 0, "Upregulated", "Downregulated"))

# Assign a sorting key before arranging
heatmap_data = heatmap_data %>%
  mutate(sort_key = ifelse(group == "Upregulated", -percent_diff, percent_diff)) %>%  # Upregulated (Descending), Downregulated (Ascending)
  arrange(group, sort_key) %>%  # Sort within each group correctly
  group_by(group) %>%
  mutate(row = (row_number() - 1) %/% 10,  # Adjust for layout
         col = (row_number() - 1) %% 10,
         label = paste0(gene, "\n", round(percent_diff, 1), "%")) %>%
  ungroup()


ggplot(heatmap_data, aes(x = col, y = -row, fill = percent_diff)) +  # Negative row for correct ordering
  geom_tile(color = "white", size = 0.5) +  # White borders for readability
  geom_text(aes(label = label), size = 3) +  # Gene name + % difference inside each square
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Blue = Down, Red = Up
  labs(title = "Heatmap of % Difference in Targeted Gene Expression",
       subtitle = "Negative % change are targeted gene knocked-down. \nPositive % change are the opposite.",
       fill = "% Difference") +
  theme_minimal() +
  facet_grid(rows = vars(group), scales = "free", space = "free_y") +  # Ensure proper spacing
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Titles for each group
    axis.title = element_blank(),  # Remove axis labels
    axis.text = element_blank(),   # Remove axis tick labels
    axis.ticks = element_blank(),  # Remove tick marks
    panel.grid = element_blank(),   # Remove background grid
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  )



#### THIS IS SCUFFED - NEED TO LOAD ENV FROM combineData...R file:

gene_plot_data = gene_counts %>%
  left_join(select(heatmap_data, gene, percent_diff), by = "gene")  # Merge percent_diff values

# Replace NA values with 0 for visualization (optional, depends on dataset)
gene_plot_data$percent_diff[is.na(gene_plot_data$percent_diff)] = 0

ggplot(gene_plot_data, aes(x = plate_row, y = plate_col, fill = percent_diff)) +
  geom_tile(color = "black") +  # Tile grid with black borders
  geom_text(aes(label = gene), color = "black", size = 3) +  # Label each tile with the gene name
  scale_fill_gradient2(
    name = "% Difference",
    low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "white"
  ) +
  facet_wrap(~ plate, scales = "free_x") +  # Group by plate
  theme_minimal() +
  theme(
    axis.text = element_blank(),  
    axis.ticks = element_blank(), 
    axis.title = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),  # Bold facet labels
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    legend.position = "right"  # Keeps legend for color scale
  ) +
  labs(
    title = "Plate Layout & Percent Difference in Knock-Down Expression Data",
    subtitle = "Most abundant genes from initial list on P6 descending to P1",
    fill = "% Change in Expression"
  )
