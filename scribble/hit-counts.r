# %%
library(tidyverse)
library(edgeR)
library(data.table)
library(purrr)

# %%
cpm = FALSE
exp.name = 'test'
analysis.name = ifelse(cpm, 'cpm', 'counts')
data.dir <- file.path("/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/selections")
save.dir = file.path("/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/analysis/", exp.name, analysis.name)
dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)

samples = data.table(
  selection = c("AG24_1", "AG24_2", "AG24_3",
                "AG24_10", "AG24_11", "AG24_12",
                "AG24_19", "AG24_20", "AG24_21"),
  group = c(
    rep('no_protein', 3),
    rep('protein', 3),
    rep('naive', 3)
  )
)

# %%
read_counts <- function(name) {
  df <- read.delim(file.path(data.dir, name, "counts.txt"), sep = "\t")
  df$name <- name
  df
}

read_dataset <- function(samples) {
  data <- map2_dfr(samples$selection, samples$group, ~ {
    d <- read_counts(.x)
    d$group <- .y
    d
  })
  
  data$group <- factor(data$group)
  data
}

get_hits.counts = function(data, cpm=FALSE) {
  
  if(cpm){
    data.counts <- data %>%
      group_by(name) %>%
      mutate(count = count / sum(count) * 1e6) %>%
      ungroup()
  }else{
    data.counts <- data
  }
  
  data.avg <- data.counts %>%
    group_by(code_1, code_2, group) %>%
    summarise(mean = mean(count), .groups = "drop")
  
  stats = data.avg %>% 
    pivot_wider(names_from = group, values_from = mean, values_fill = 0) %>%
    mutate(enrichment = protein - no_protein, sticky = no_protein - naive)
  
  hits <- stats %>%
    arrange(desc(enrichment)) %>%
    head(100)
  
  sticky <- stats %>%
    arrange(desc(sticky)) %>%
    head(100)
  
  result = list(
    stats = stats,
    hits = hits,
    sticky = sticky
  )
  
}

data = read_dataset(selections, grps)
result.counts = get_hits.counts(data=data, cpm=cpm)

save.path = file.path(save.dir, 'stats.csv')
write_csv(result.counts$stats, file=save.path)

save.path = file.path(save.dir, 'hits.csv')
write_csv(result.counts$hits, file=save.path)

save.path = file.path(save.dir, 'sticky.csv')
write_csv(result.counts$sticky, file=save.path)

uniq_grps = unique(groups)
for(grp in uniq_grps){
  fname = paste0(grp, '.csv')
  save.path = file.path(save.dir, fname)
  result.counts$stats %>% 
    select(code_1, code_2, all_of(grp)) %>%
    mutate(count = .data[[grp]]) %>% 
    select(-all_of(grp)) %>%
    write_csv(file=save.path)
}



