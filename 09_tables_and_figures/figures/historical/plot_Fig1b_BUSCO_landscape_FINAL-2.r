# Fig. 1b — final two-row publication layout
# Top row: Mucoromycotina, y-axis 55–100%
# Bottom row: Mortierellomycotina + Glomeromycotina, y-axis 35–100%
# Run all in RStudio, then select the local TSV/CSV/XLSX table.

pkgs <- c("ggplot2", "dplyr", "tidyr", "readr", "readxl", "tools", "patchwork")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(tools)
library(patchwork)

input_file <- file.choose()
ext <- tolower(file_ext(input_file))

if (ext %in% c("tsv", "txt")) {
  dat <- read_tsv(input_file, show_col_types = FALSE, na = c("", "NA", "NaN"))
} else if (ext == "csv") {
  dat <- read_csv(input_file, show_col_types = FALSE, na = c("", "NA", "NaN"))
} else if (ext %in% c("xlsx", "xls")) {
  dat <- read_excel(input_file, na = c("", "NA", "NaN"))
} else {
  stop("Please select a TSV, CSV, XLSX, or XLS file.")
}

required <- c(
  "subphylum", "species_no", "species_name_current",
  "genome_BUSCO_C", "public_protein_BUSCO_C", "braker_protein_BUSCO_C"
)
missing_cols <- setdiff(required, names(dat))
if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

busco_cols <- c("genome_BUSCO_C", "public_protein_BUSCO_C", "braker_protein_BUSCO_C")
sub_order <- c("Mucoromycotina", "Mortierellomycotina", "Glomeromycotina")

dat <- dat %>%
  mutate(
    across(all_of(busco_cols), as.numeric),
    species_name_current = trimws(as.character(species_name_current)),
    subphylum = factor(trimws(as.character(subphylum)), levels = sub_order)
  ) %>%
  arrange(subphylum, species_name_current, species_no) %>%
  group_by(subphylum) %>%
  mutate(x = row_number()) %>%
  ungroup()

if (nrow(dat) != 202) warning("Input contains ", nrow(dat), " rows rather than 202.")
if (anyDuplicated(dat$species_no)) warning("Duplicated species_no values detected.")
if (any(vapply(dat[busco_cols], function(z) any(z < 0 | z > 100, na.rm = TRUE), logical(1)))) {
  stop("At least one BUSCO value is outside 0–100.")
}

series_colours <- c(
  "Genome" = "#4D4D4D",
  "Public proteins" = "#0072B2",
  "BRAKER proteins" = "#D55E00"
)
series_shapes <- c(
  "Genome" = 16,
  "Public proteins" = 21,
  "BRAKER proteins" = 17
)
series_fills <- c(
  "Genome" = "#4D4D4D",
  "Public proteins" = "white",
  "BRAKER proteins" = "#D55E00"
)

make_panel <- function(panel_data, y_min, show_y = TRUE, show_species = FALSE) {
  panel_data <- panel_data %>%
    arrange(species_name_current, species_no) %>%
    mutate(x = row_number())

  panel_name <- as.character(panel_data$subphylum[1])
  panel_title <- paste0(panel_name, "  ·  n = ", nrow(panel_data))

  segment_data <- panel_data %>%
    rowwise() %>%
    mutate(
      ymin = min(c_across(all_of(busco_cols)), na.rm = TRUE),
      ymax = max(c_across(all_of(busco_cols)), na.rm = TRUE)
    ) %>%
    ungroup()

  long_data <- panel_data %>%
    select(species_no, species_name_current, x, all_of(busco_cols)) %>%
    pivot_longer(all_of(busco_cols), names_to = "series", values_to = "BUSCO_C") %>%
    filter(!is.na(BUSCO_C)) %>%
    mutate(series = factor(
      series,
      levels = busco_cols,
      labels = c("Genome", "Public proteins", "BRAKER proteins")
    ))

  if (y_min >= 55) {
    y_breaks <- c(60, 70, 80, 90, 100)
  } else if (y_min > 0) {
    y_breaks <- c(40, 50, 60, 70, 80, 90, 100)
  } else {
    y_breaks <- seq(0, 100, 20)
  }

  reference_lines <- c(80, 90, 95)
  reference_lines <- reference_lines[reference_lines >= y_min]

  p <- ggplot() +
    geom_hline(
      yintercept = reference_lines,
      colour = c("#D9D9D9", "#D9D9D9", "#E8E8E8")[seq_along(reference_lines)],
      linewidth = c(0.27, 0.27, 0.23)[seq_along(reference_lines)],
      linetype = c("22", "22", "dotted")[seq_along(reference_lines)]
    ) +
    geom_segment(
      data = segment_data,
      aes(x = x, xend = x, y = ymin, yend = ymax),
      colour = "#BDBDBD", linewidth = 0.23, alpha = 0.68
    ) +
    geom_point(
      data = long_data,
      aes(x = x, y = BUSCO_C, colour = series, shape = series, fill = series),
      size = 1.55, stroke = 0.68, alpha = 0.97
    ) +
    scale_colour_manual(values = series_colours, name = NULL, drop = FALSE) +
    scale_shape_manual(values = series_shapes, name = NULL, drop = FALSE) +
    scale_fill_manual(values = series_fills, name = NULL, drop = FALSE) +
    scale_y_continuous(
      breaks = y_breaks,
      minor_breaks = NULL,
      expand = expansion(mult = c(0.012, 0.025))
    ) +
    coord_cartesian(ylim = c(y_min, 100), clip = "off") +
    labs(
      title = panel_title,
      x = NULL,
      y = if (show_y) "Complete BUSCOs (%)" else NULL
    ) +
    guides(
      colour = guide_legend(nrow = 1, byrow = TRUE,
                            override.aes = list(size = 3, alpha = 1)),
      shape = guide_legend(nrow = 1, byrow = TRUE),
      fill = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    theme_classic(base_size = 9, base_family = "sans") +
    theme(
      plot.title = element_text(
        face = "bold", colour = "#222222", size = 9.5,
        hjust = 0.5, margin = margin(0, 0, 4, 0)
      ),
      panel.grid.major.y = element_line(colour = "#EFEFEF", linewidth = 0.27),
      panel.border = element_rect(colour = "#8C8C8C", fill = NA, linewidth = 0.36),
      axis.title = element_text(colour = "#202020", size = 9),
      axis.text.y = element_text(colour = "#333333", size = 7.5),
      axis.line = element_line(colour = "#333333", linewidth = 0.40),
      axis.ticks = element_line(colour = "#333333", linewidth = 0.32),
      legend.position = "top",
      legend.text = element_text(size = 8.5),
      legend.spacing.x = grid::unit(0.20, "cm"),
      legend.key.width = grid::unit(0.48, "cm"),
      legend.key.height = grid::unit(0.42, "cm"),
      plot.margin = margin(4, 6, 4, 5)
    )

  if (!show_y) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  if (show_species) {
    p +
      scale_x_continuous(
        breaks = panel_data$x,
        labels = panel_data$species_name_current,
        expand = expansion(mult = c(0.004, 0.004))
      ) +
      theme(
        axis.text.x = element_text(
          angle = 90, hjust = 1, vjust = 0.5,
          size = 4.2, face = "italic"
        )
      )
  } else {
    p +
      scale_x_continuous(
        breaks = NULL,
        expand = expansion(mult = c(0.004, 0.004))
      ) +
      theme(axis.ticks.x = element_blank())
  }
}

assemble_final <- function(show_species = FALSE) {
  p_muc <- make_panel(
    filter(dat, subphylum == "Mucoromycotina"),
    y_min = 55,
    show_y = TRUE,
    show_species = show_species
  )

  p_mort <- make_panel(
    filter(dat, subphylum == "Mortierellomycotina"),
    y_min = 35,
    show_y = TRUE,
    show_species = show_species
  )

  p_glom <- make_panel(
    filter(dat, subphylum == "Glomeromycotina"),
    y_min = 35,
    show_y = FALSE,
    show_species = show_species
  )

  lower_row <- p_mort + p_glom + plot_layout(widths = c(57, 22))

  (p_muc / lower_row) +
    plot_layout(heights = c(1, 1), guides = "collect") +
    plot_annotation(
      caption = paste0(
        "Genomes are ordered alphabetically within each subphylum. ",
        "The upper and lower rows use y-axis ranges of 55–100% and 35–100%, respectively."
      ),
      theme = theme(
        plot.caption = element_text(
          hjust = 0.5, size = 7.8, colour = "#444444",
          margin = margin(5, 0, 0, 0)
        ),
        plot.margin = margin(3, 4, 3, 4)
      )
    ) &
    theme(legend.position = "top", legend.justification = "center")
}

p_final <- assemble_final(show_species = FALSE)
p_labelled <- assemble_final(show_species = TRUE)

output_dir <- dirname(normalizePath(input_file))

ggsave(
  file.path(output_dir, "Fig1b_BUSCO_FINAL_two_rows.pdf"),
  p_final, width = 7.2, height = 6.6, units = "in"
)
ggsave(
  file.path(output_dir, "Fig1b_BUSCO_FINAL_two_rows.tiff"),
  p_final, width = 7.2, height = 6.6, units = "in",
  dpi = 600, compression = "lzw"
)
ggsave(
  file.path(output_dir, "Fig1b_BUSCO_FINAL_two_rows_large_preview.pdf"),
  p_final, width = 14.4, height = 10.2, units = "in"
)
ggsave(
  file.path(output_dir, "Fig1b_BUSCO_FINAL_two_rows_labelled.pdf"),
  p_labelled, width = 30, height = 13, units = "in",
  limitsize = FALSE
)

cat(
  "\nSaved beside the selected table:\n",
  "  Fig1b_BUSCO_FINAL_two_rows.pdf\n",
  "  Fig1b_BUSCO_FINAL_two_rows.tiff\n",
  "  Fig1b_BUSCO_FINAL_two_rows_large_preview.pdf\n",
  "  Fig1b_BUSCO_FINAL_two_rows_labelled.pdf\n",
  sep = ""
)

print(p_final)
