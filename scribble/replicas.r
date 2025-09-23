library(tidyverse)
library(GGally)

data.dir <- "/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/selections"
save.dir = file.path("/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/analysis/")
dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)

neg_ctrls <- c("AG24_1", "AG24_2", "AG24_3")
naive     <- c("AG24_19", "AG24_20", "AG24_21")
protein   <- c("AG24_10", "AG24_11", "AG24_12")

read_counts <- function(name, group) {
  df <- read.delim(file.path(data.dir, name, "counts.txt"), sep = "\t")
  df$name <- name
  df$group <- group
  df
}

# Read all data
data <- bind_rows(
  lapply(neg_ctrls, read_counts, group = "neg_ctrl"),
  lapply(naive,     read_counts, group = "naive"),
  lapply(protein,   read_counts, group = "protein")
)

# Quick check: total counts per sample
data %>%
  group_by(name) %>%
  summarise(total = sum(count))

# Plot replicate comparisons per group
for (grp in unique(data$group)) {
  df_group <- data %>%
    filter(group == grp) %>%
    select(code_1, code_2, name, count) %>%
    pivot_wider(names_from = name, values_from = count, values_fill = 0)
  
  g = df_group %>%
    select(-code_1, -code_2) %>%
    ggpairs(
      upper = list(continuous = wrap("smooth", alpha = 0.3, size = 0.2)),
      lower = list(continuous = wrap("cor", size = 3))
    ) +
    ggtitle(paste("Replicate comparisons for", grp))
  ggsave(file.path(save.dir, paste0("replicates_", grp, ".png")), g, width = 8, height = 6)
}
