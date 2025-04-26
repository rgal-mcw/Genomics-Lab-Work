# Ryan Gallagher
# 03.12.2025
#
# This file will manipulate the EditCo Knock-Down experiment data

library(tidyverse)

# Read in files
files = list.files(path='data', pattern='combined_counts_Batch[1-3]\\.txt', full.names=TRUE)
data_list = lapply(files, read.table, header=TRUE, sep='\t')
result = do.call(cbind, data_list)

# Build metadata
meta = tibble(
  sample = names(result),  # Full sample name
  gene = sub("\\..*", "", names(result)),  # Everything before the first '.'
) %>%
  mutate(
    num_dots = str_count(sample, "\\."),  # Count number of dots in the sample name
    plate = case_when(
      num_dots == 2 ~ str_extract(sample, "(?<=\\.)[^._]+"),  # After first '.' and before first '_'
      num_dots == 1 ~ str_extract(sample, "P\\d+"),  # Extract P# directly
      TRUE ~ NA_character_
    ),
    well = case_when(
      num_dots == 2 ~ str_extract(sample, "(?<=_)[^_]+(?=_)"),  # After first '_' and before last '_'
      num_dots == 1 ~ str_extract(sample, "(?<=P\\d)[A-H]\\d{2}"),  # After P# and before first '_'
      TRUE ~ NA_character_
    )
  ) %>%
  select(-num_dots) %>%
  mutate(
    condition = if_else(gene %in% c("WT", "TRAC"), "control", "treatment")
  )


# Write our files
write.csv(result, './data/counts.csv', row.names=T)
write.csv(meta, './data/meta.csv', row.names=FALSE)


# Define grid size (adjust columns per row)
gene_counts = tibble(gene = sub("\\..*", "", names(result))) %>%
  count(gene, name = "count") %>%
  left_join(meta %>% select(gene, plate), by = "gene") %>%
  distinct()

num_cols = 10 
gene_counts <- tibble(gene = sub("\\..*", "", names(result))) %>%
  count(gene, name = "count") %>%
  left_join(meta %>% select(gene, plate), by = "gene") %>%
  group_by(plate, gene) %>%  # Ensure each gene appears only once per plate
  summarise(count = max(count), .groups = "drop") %>%  # Keep only max count per gene-plate
  arrange(plate, count == 18, gene) %>%
  group_by(plate) %>%  
  mutate(
    plate_row = (row_number() - 1) %/% num_cols + 1,  # Compute row index within plate
    plate_col = (row_number() - 1) %% num_cols + 1    # Compute column index within plate
  ) %>%
  ungroup()

# Plot as a grid
ggplot(gene_counts, aes(x = plate_row, y = plate_col, fill = factor(count))) +
  geom_tile(color = "black") +
  geom_text(aes(label = gene), color = "black", size = 3) +  
  scale_fill_manual(
    name = "Replicate Count",
    values = c("1" = "red", "2" = "yellow", "3" = "green", "18" = "gray"), 
    labels = c("1" = "1", "2" = "2", "3" = "3", "18" = "control"),  
    na.value = "white"
  ) +
  facet_wrap(~ plate, scales = "free_x") +  # Ensure each plate is grouped
  theme_minimal() +
  theme(axis.text = element_blank(),  
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=15, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9)) + # Make plate labels more visible
  labs(title = "Plate Layout & Replicate Counts in our Knock-Down Expression Data",
       subtitle = "Most abundant genes from initial list on P6 descending to P1")
