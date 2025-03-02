library(phytools)
library(ggimage)

tar_load(c(combined_sanger_tree, traits, samples_with_subgen), store = "data/_targets-2025-02-26")

holo_epi_traits <-
  traits |>
  mutate(holo_epi = case_when(
    ecology == "t" ~ "not_epiphyte",
    ecology == "he" ~ "not_epiphyte",
    growth_form == "i" ~ "individual_epiphyte",
    growth_form == "c" ~ "colonial_epiphyte",
    .default = NA_character_
  )) |>
  select(genus, specific_epithet, holo_epi)

holo_epi_tips <-
  tibble(tip = combined_sanger_tree$tip.label) |>
  left_join(samples_with_subgen, by = "tip", relationship = "one-to-one") |>
  select(tip, genus, specific_epithet, species) |>
  left_join(
    holo_epi_traits,
    by = c("genus", "specific_epithet"),
    relationship = "many-to-one"
  ) |>
  filter(!is.na(holo_epi)) |>
  mutate(
    keep = case_when(
      genus == "Hymenophyllum" & tip != "Hymenophyllum_apiculatum" ~ FALSE,
      .default = TRUE
    )
  ) |>
  filter(keep) |>
  select(-keep) |>
  # Keep only one sample per species
  slice_head(n = 1, by = species)

tips_keep <- intersect(holo_epi_tips$tip, combined_sanger_tree$tip.label)

tree <- ape::keep.tip(combined_sanger_tree, tips_keep) |>
  ape::root("Hymenophyllum_apiculatum") |>
  ladderize()

holo_epi_tips_use <- holo_epi_tips |>
  filter(tip %in% tips_keep)

holo_epi_tips_vec <-
  holo_epi_tips_use |>
  pull(holo_epi) |>
  set_names(holo_epi_tips_use$tip)

epi_mod_er <- fitMk(tree, holo_epi_tips_vec, model = "ER")
epi_mod_sym <- fitMk(tree, holo_epi_tips_vec, model = "SYM")
epi_mod_ard <- fitMk(tree, holo_epi_tips_vec, model = "ARD")

epi_mod_aov <- anova(epi_mod_er, epi_mod_sym, epi_mod_ard)

epi_mod_ancr <- ancr(epi_mod_ard)

plot(epi_mod_ancr)

anc_state_tbl <- epi_mod_ancr$ace |>
  as.data.frame() |>
  rownames_to_column("node") |>
  as_tibble() |>
  mutate(node = parse_number(node))

# Do not change the order of elements in this vector!
# hard-coded for pie cols to match tip cols

# To preview colors
# scales::show_col(palette.colors()) # nolint

trait_colors <- c(
  colonial_epiphyte = "#009E73", # grey in Dubuisson
  individual_epiphyte = "#0072B2", # black in Dubuission
  not_epiphyte = "#D55E00" # white in Dubuisson
)

trait_labels <- c(
  colonial_epiphyte = "Colonial epiphyte",
  individual_epiphyte = "Individual epiphyte",
  not_epiphyte = "Not epiphyte"
)

pies <- nodepie(
  anc_state_tbl,
  cols = 2:4,
  color = trait_colors,
  outline.color = "black",
  outline.size = 0.2
  )

tree <- attributes(epi_mod_ancr)$tree |>
  rescale_tree_len(100)

 ggtree(tree) + xlim(0, 120)

dna_tree_plot <- ggtree(tree) %<+%
  holo_epi_tips +
  # add extra horizontal space for species names
  xlim(0, 130) +
  geom_tiplab(aes(label = species), size = 5, offset = 1.5) +
  geom_tippoint(aes(fill = holo_epi), size = 5, shape = 22, color = "black") +
  scale_fill_manual(
    values = trait_colors, name = "Growth habit",
    breaks = names(trait_labels),
    labels = trait_labels
  ) +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.position.inside = c(1, 0)
  )

anc_res_plot <- inset(dna_tree_plot, pies, width = 0.04, height = 0.04)

ggsave(plot = anc_res_plot, file = "images/anc_res_plot.png", height = 11, width = 12)
