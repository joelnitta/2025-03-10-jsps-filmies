# Rescale a tree
rescale_tree_len <- function(tree, scale) {
  tree$edge.length <- tree$edge.length /
    max(phytools::nodeHeights(tree)) *
    scale
  tree
}

# Create a bootstrap tibble for plotting tree
make_bs_tbl <- function(tree, cutoff = 0) {
  # Make a dataframe (tibble) with node IDs (integers) and their corresponding bootstrap support values
  # The tibble has two columns: one called "node", the other can be named as we like (here, "bootstrap")
  bs_tibble <- tibble(
    # hard-code node ID: internal nodes start after tip nodes,
    # and phy$node.label is in the same order as internal nodes
    node = 1:Nnode(tree) + Ntip(tree),
    bootstrap = parse_number(tree$node.label)
  )
  if (cutoff > 0) {
    bs_tibble <-
      bs_tibble |>
      mutate(bootstrap = if_else(bootstrap < cutoff, "", bootstrap))
  }
  bs_tibble
}

# Plot full sampling nuclear tree
plot_full_nuclear <- function(hybpiper_dna_tree, samples_with_subgen) {
  hybpiper_dna_tree_plotting <- hybpiper_dna_tree

  hybpiper_dna_tree_data <- tibble(
    tip = hybpiper_dna_tree$tip.label
  ) |>
    left_join(
      samples_with_subgen,
      relationship = "one-to-one",
      by = join_by(tip)
    )

  tips_hym <- hybpiper_dna_tree_data |>
    filter(genus == "Hymenophyllum") |>
    pull(tip)

  hybpiper_dna_tree_plotting <-
    hybpiper_dna_tree_plotting |>
    root(tips_hym) |>
    rescale_tree_len(100)

  # Prepare bootstrap tbl
  bs_tibble <- make_bs_tbl(hybpiper_dna_tree_plotting)

  dna_tree_plot <- ggtree(hybpiper_dna_tree_plotting) %<+%
    bs_tibble %<+%
    hybpiper_dna_tree_data +
    # add extra horizontal space for species names
    xlim(0, 120)

  dna_tree_plot +
    geom_tiplab(aes(label = species), size = 3) +
    geom_nodelab(aes(label = bootstrap), hjust = -0.1, size = 3) +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = glue::glue("GoFLAG マーカー（{n_markers_nuc}個）"),
      size = 6,
      hjust = 0.1,
      vjust = 1,
      family = "HiraKakuPro-W3"
    ) +
    annotate(
      "text",
      x = Inf,
      y = 0,
      label = "iqtreeによるML法\n1000回ブートストラップ",
      size = 6,
      hjust = 1,
      vjust = -0.5,
      family = "HiraKakuPro-W3"
    )
}

prep_clade_labels <- function(tree, samples_with_subgen, tax_level) {
  clade_data <- tibble(
    tip = tree$tip.label
  ) |>
    left_join(samples_with_subgen, relationship = "one-to-one", by = "tip") |>
    filter(!is.na({{tax_level}}))

  singleltons <- clade_data |>
    add_count({{tax_level}}) |>
    filter(n == 1) |>
    select(tip, {{tax_level}}) |>
    mutate(mrca = match(tip, tree$tip.label)) |>
    select(-tip)

  clade_data |>
    add_count({{tax_level}}) |>
    filter(n > 1) |>
    group_by({{tax_level}}) |>
    summarize(mrca = getMRCA(tree, tip)) |>
    select(mrca, {{tax_level}}) |>
    bind_rows(singleltons)
}
