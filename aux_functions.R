# aux function for EukBank analysis

library(Biostrings)
library(tidyverse)
library(speedyseq)
library(patchwork)
library(ggtext)
library(ComplexHeatmap)
library(viridis)
library(ggtree)
library(treeio)
library(phytools)

theme_set(theme_minimal(base_size = 12))

envplot_colors <-
  set_names(
    c(
      "royalblue3", "#00A04B",
      "#FFB400", "#FF5A5F", "#7B0051"
    ),
    nm = c(
      "marine_water", "marine_sediment",
      "freshwater", "freshwater_sediment", "soil"
    )
  )

envplot_plotter <- function(df, facet, x = "envplot", y = ".abundance", scales_ = "free") {
  plot <-
    df %>%
    mutate(envplot = factor(envplot, levels = c(
      "marine_water", "marine_sediment", "marine_ice",
      "soil", "land_water", "freshwater", "freshwater_sediment"
    ))) %>%
    select(
      "envplot" = all_of(x),
      "abundance" = all_of(y),
      "facet" = all_of(facet),
      everything()
    ) %>%
    ggplot(aes(
      x = envplot,
      y = abundance
    )) +
    geom_jitter(aes(color = envplot),
      width = 0.2,
      height = 0,
      alpha = 0.1
    ) +
    geom_violin(aes(fill = envplot),
      scale = "width",
      alpha = 0.5
    ) +
    scale_color_manual(values = envplot_colors) +
    scale_fill_manual(values = envplot_colors) +
    facet_wrap(~facet, scales = scales_) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      strip.text = element_text(size = 12)
    ) +
    labs(
      x = "",
      y = "Relative read abundance (%)"
    )

  return(plot)
}

read_blast <- function(path) {
  read_tsv(path, col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
}

df_to_fasta <- function(df, out_file, header_col = "amplicon", seq_col = "sequence") {
  df %>%
    select(all_of(c(header_col, seq_col))) %>%
    deframe() %>%
    DNAStringSet() %>%
    writeXStringSet(out_file)
}

factor_envplot <- function(df) {
  df %>%
    mutate(envplot = factor(envplot, levels = c(
      "marine_water", "marine_sediment",
      "freshwater", "freshwater_sediment", "soil"
    )))
}

remove_empty_taxa <- function(physeq) {
  require(phyloseq)
  prune_taxa(x = physeq, taxa = taxa_sums(physeq) > 0)
}

remove_empty_samples <- function(physeq) {
  prune_samples(x = physeq, samples = sample_sums(physeq) > 0)
}

counts_to_percent <- function(physeq, total_col) {
  # transform counts to percentage but using a total defined in sample_data

  otu_tab <-
    as.matrix(otu_table(physeq))

  totals <- sample_data(physeq)[[total_col]]

  if (identical(rownames(sample_data(physeq)), colnames(otu_tab))) {
    new_otu_tab <-
      sweep(100 * otu_tab, 2, totals, `/`)

    new_physeq <- physeq

    otu_table(new_physeq) <- otu_table(new_otu_tab, taxa_are_rows = T)

    return(new_physeq)
  } else {
    message("Sample names order does not match betwen OTU and metadata tables")
  }
}

tax_glom_custom <- function(physeq, tax.col) { # adapted from https://github.com/joey711/phyloseq/issues/941#issuecomment-391018770

  require(tidyverse)

  reorderedtax <- data.frame(tax_table(physeq)) %>%
    # Puts the selected column first with everything else aftwds
    select(one_of(tax.col), everything()) %>%
    as.matrix() %>%
    tax_table()

  tax_table(physeq) <- reorderedtax

  return(tax_glom(physeq, taxrank = tax.col))
}

shared_asvs_plotter <- function(taxa) {
  df <- # fix axis
    asvs_diff_envplots %>%
    mutate(name = if_else(name == "mean_abun", "Mean relative\nabundance (%)", "Occurrence (%)")) %>%
    filter(tax == taxa) %>%
    bind_rows(tibble(
      envplot = c("marine_water", "marine_sediment", "soil", "freshwater", "freshwater_sediment"),
      name = c(rep("Mean relative\nabundance (%)", 3), rep("Occurrence (%)", 2)),
      taxogroup3 = rep(taxa, 5)
    ) %>%
      factor_envplot())
  
  p_abun <-
    ggplot() +
    geom_jitter(
      data = df %>% filter(type != "shared" | is.na(type), name == "Mean relative\nabundance (%)"),
      aes(
        color = color_points,
        x = envplot, y = value
      ),
      width = 0.2, height = 0
    ) +
    scale_x_discrete(
      drop = FALSE,
      labels = ~ str_to_sentence(str_replace(.x, "_", " "))
    ) +
    geom_jitter(
      data = df %>% filter(type == "shared", name == "Mean relative\nabundance (%)"), aes(x = envplot, y = value),
      color = "black", shape = 21, fill = "white", width = 0.2, height = 0
    ) +
    geom_boxplot(data = df %>% filter(name == "Mean relative\nabundance (%)"), aes(fill = envplot, x = envplot, y = value), alpha = 0.5, outlier.shape = NA) +
    scale_color_manual(values = envplot_colors) +
    scale_fill_manual(values = envplot_colors) +
    scale_y_log10(
      labels = c("10<sup>-8</sup>", "10<sup>-6</sup>", "10<sup>-4</sup>", 0.01, 0.1, 1),
      limits = c(
        min(asvs_diff_envplots %>% filter(name == "mean_abun") %>% .$value),
        max(asvs_diff_envplots %>% filter(name == "mean_abun") %>% .$value)
      ),
      breaks = c(10^c(-8, -6, -4, -2, -1, 0))
    ) +
    guides(fill = "none", color = "none") +
    theme_bw() +
    labs(x = "", y = "Mean relative\nabundance (%)", title = taxa) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_markdown()
    )
  
  p_occ <-
    ggplot() +
    geom_jitter(
      data = df %>% filter(type != "shared" | is.na(type), name == "Occurrence (%)"),
      aes(
        color = color_points,
        x = envplot, y = value
      ),
      width = 0.2, height = 0
    ) +
    scale_x_discrete(
      drop = FALSE,
      labels = ~ str_to_sentence(str_replace(.x, "_", " "))
    ) +
    geom_jitter(
      data = df %>% filter(type == "shared", name == "Occurrence (%)"), aes(x = envplot, y = value),
      color = "black", shape = 21, fill = "white", width = 0.2, height = 0
    ) +
    geom_boxplot(data = df %>% filter(name == "Occurrence (%)"), aes(fill = envplot, x = envplot, y = value), alpha = 0.5, outlier.shape = NA) +
    scale_color_manual(values = envplot_colors) +
    scale_fill_manual(values = envplot_colors) +
    scale_y_log10(
      labels = c(0.01, 0.1, 1, 10, 100),
      limits = c(
        min(asvs_diff_envplots %>% filter(name == "occ_perc") %>% .$value),
        max(asvs_diff_envplots %>% filter(name == "occ_perc") %>% .$value)
      ),
      breaks = c(0.01, 0.1, 1, 10, 100)
    ) +
    guides(fill = "none", color = "none") +
    theme_bw() +
    labs(x = "", y = "Occurrence (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  
  p_abun / p_occ
}

plot_mast_tree <- function(group, offset = 0.1){
  
  tree <- 
    trees[[group]]
  
  tree_names <-  
    tibble(name = tree$tip.label) |> 
    mutate(clade = if_else(str_detect(name, '-outgroup'), 'Outgroup', str_remove(name, '.*_')),
           clade = if_else(str_detect(clade, '[0-9]$'), NA_character_, clade)) |> # remove references without subclade
    filter(!is.na(clade))
  
  clade_node_numbers <- 
    tree_names |> 
    group_by(clade) |> 
    summarise(node_num = MRCA(tree, name))
  
  outgroup_num <- 
    clade_node_numbers |> 
    filter(clade == 'Outgroup') |> 
    pull(node_num)
  
  len_root <- 
    tree$edge.length[tree$edge[, 2] == outgroup_num]/2
  
  rooted_tree <- 
    reroot(tree, node.number = outgroup_num, position = len_root)
  
  node_num_all <- # first node number of tree without outgroup
    MRCA(rooted_tree, tree_names |> filter(clade != 'Outgroup') |> pull(name))
  
  clade_node_numbers_rooted <- 
    tree_names |> 
    group_by(clade) |> 
    summarise(node_num = MRCA(rooted_tree, name)) |> 
    filter(node_num != node_num_all)
  
  bootstrap_values_filtered <-
    as_tibble(rooted_tree) |>
    filter(!label %in% rooted_tree$tip.label) |>
    mutate(label = if_else(label == 'Root', '', label),
           shalrt = as.numeric(str_remove(label, '\\/.*')),
           ufboot = as.numeric(str_remove(label, '.*\\/')),
           # label = case_when(node %in% clade_node_numbers_rooted$node_num ~ label,
           #                   shalrt >= 80 & ufboot >= 95 ~ label,
           #                   TRUE ~ '')
    ) 
  
  rooted_tree$node.label <- bootstrap_values_filtered$label
  
  p_tree <- 
    ggtree(rooted_tree) + 
    # geom_text(aes(label=node), hjust=-.3) +
    geom_tiplab(size = 2.5) +
    geom_nodelab(size = 3) +
    geom_treescale()
  
  for (cla in clade_node_numbers_rooted |> filter(clade != 'Outgroup') |> pull(clade)){
    
    num_node <- 
      clade_node_numbers_rooted |> 
      filter(clade == cla) |> 
      pull(node_num)
    
    p_tree <- 
      p_tree +
      geom_cladelab(node = num_node, label = cla, align = T, offset = offset, size = 4)
    
  }
  
  return(p_tree)
}

vegan_formatter <- function(df, row_names, sample_col = 'Sample', abun_col = 'Abundance', fill = 0){
  
  df_wide <- 
    df %>% 
    dplyr::select(sample_col, abun_col, row_names) %>% 
    pivot_wider(names_from = sample_col,
                values_from = abun_col, 
                values_fill = fill) %>% 
    column_to_rownames(row_names) %>% 
    t()
  
  return(df_wide)
  
}