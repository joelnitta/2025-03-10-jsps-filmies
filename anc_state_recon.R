suppressPackageStartupMessages(library(tidyverse))
library(targets)
library(ape)
library(ggtree)
library(ggimage)
library(phytools)

source("functions.R")

# Define location of targets cache from goflag filmies analysis
filmy_store <- "data/_targets_2025-03-04"

tar_load(c(anc_recon, holo_epi_tips, samples_with_subgen), store = filmy_store)

samples_with_subgen <- mutate(
  samples_with_subgen,
  subgenus = str_replace_all(subgenus, "Davalliopsis \\?", "Pachychaetum")
) |>
  mutate(
    subgenus = case_when(
      species == "Didymoglossum fulgens" ~ "Microgonium",
      species == "Didymoglossum sublimbatum" ~ "Microgonium",
      .default = subgenus
    )
  )

# Extract tibble of ancestral state percentages
anc_state_tbl <- anc_recon$ace |>
  as.data.frame() |>
  rownames_to_column("node") |>
  as_tibble() |>
  mutate(node = parse_number(node)) |>
  pivot_longer(names_to = "trait", values_to = "prob", -node) |>
  mutate(n_states = str_count(trait, "\\+") + 1) |>
  separate_rows(trait, sep = "\\+") |>
  mutate(prob = prob / n_states) |>
  summarize(prob = sum(prob), .by = c(node, trait)) |>
  pivot_wider(names_from = trait, values_from = prob, id_cols = node)

# Set trait colors
# Do not change the order of elements in this vector!
# they are hard-coded for pie cols to match tip cols
trait_colors <- c(
  terrestrial = "#D55E00",
  epi_tree_fern = "#F0E442",
  epiphytic = "#009E73",
  epipetric = "#0072B2"
)

# To preview colors
# scales::show_col(palette.colors()) # nolint

# Set trait labels
trait_labels <- c(
  terrestrial = "Terrestrial",
  epiphytic = "Epiphytic",
  epi_tree_fern = "Epi. on tree fern",
  epipetric = "Epipetric"
)

# Rescale tree
tree <- attributes(anc_recon)$tree |>
  rescale_tree_len(100)

# Make pie charts
# - node pie charts
node_pies <-
  anc_state_tbl |>
  select(node, all_of(names(trait_colors))) |>
  ggtree::nodepie(
    cols = 2:4,
    color = trait_colors,
    outline.color = "black",
    outline.size = 0.2
  )

# - tip pie charts
curr_state_tbl <-
  holo_epi_tips |>
  select(tip, trait = growth_habit) |>
  mutate(n_states = str_count(trait, "\\+") + 1) |>
  separate_rows(trait, sep = "\\+") |>
  mutate(prob = 1 / n_states) |>
  pivot_wider(
    names_from = trait,
    values_from = prob,
    id_cols = tip,
    values_fill = 0
  )

tip_pies <- tibble(tip = tree$tip.label) %>%
  mutate(node = 1:nrow(.)) |>
  left_join(curr_state_tbl, by = "tip") |>
  select(node, all_of(names(trait_colors))) |>
  ggtree::nodepie(
    cols = 2:5,
    color = trait_colors,
    outline.color = "black",
    outline.size = 0.2
  )
# - combine
pies <- c(node_pies, tip_pies)

# Make simplified growth habit tip data for legend
growth_habit_simple <- holo_epi_tips |>
  select(tip, species, growth_habit) |>
  separate(growth_habit, "growth_habit", sep = "\\+", extra = "drop") |>
  assertr::verify(n_distinct(growth_habit) == 4)

# Prepare clade labels
subgenus_labs <- prep_clade_labels(tree, samples_with_subgen, subgenus) |>
  # Drop Hymenophyllum subgen
  filter(subgenus != "Mecodium")

adj <- 50

genus_labs <- prep_clade_labels(tree, samples_with_subgen, genus) |>
  # Set offset by whether the genus has a subgenus or not
  left_join(
    unique(select(samples_with_subgen, genus, subgenus)),
    by = "genus"
  ) |>
  summarize(
    n_subgenus = n_distinct(subgenus, na.rm = TRUE),
    .by = c(mrca, genus)
  ) |>
  mutate(
    offset_dist = case_when(
      # Special case for outgroup Hymenophyllum: don't show subgenus
      genus == "Hymenophyllum" ~ 3 + adj,
      n_subgenus > 0 ~ 28 + adj,
      n_subgenus == 0 ~ 3 + adj
    )
  )

# Make base plot with pie charts
dna_tree_plot <- ggtree::ggtree(tree) %<+%
  growth_habit_simple

base_plot <- inset(dna_tree_plot, pies, width = 0.04, height = 0.04)

# Add rest of features
plot <- base_plot +
  # Need extra horizontal space for labels
  xlim(0, 200) +
  geom_tiplab(aes(label = species), size = 3, offset = 3) +
  geom_cladelab(
    data = genus_labs,
    mapping = aes(node = mrca, label = genus, offset = offset_dist),
    hjust = -0.1,
    extend = 0.25,
    barsize = 1,
    fontsize = 4
  ) +
  geom_cladelab(
    data = subgenus_labs,
    mapping = aes(node = mrca, label = subgenus),
    offset = 45,
    hjust = -0.1,
    extend = 0.25,
    barsize = 0.5,
    fontsize = 4
  ) +
  geom_tippoint(
    aes(fill = growth_habit),
    size = 0,
    shape = 22,
    color = "black",
    alpha = 0
  ) +
  scale_fill_manual(
    values = trait_colors,
    name = "Growth habit",
    breaks = names(trait_labels),
    labels = trait_labels
  ) +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.position = c(0.15, 0.95)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1)))

ggsave(
  plot = plot,
  file = "images/anc_recon_plot.png",
  height = 16,
  width = 10,
  units = "in"
)
