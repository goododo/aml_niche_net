# theme_lcc.R ----
# One figure theme for the whole patch, sourced by 13 / 14 / 15 so the 15 figures cannot drift apart.
#
# PRESENTATION RULES, set here once:
#   1. FIGURES CARRY DATA ONLY. Captions are suppressed at render time -- they stay in each script's
#      labs() call as provenance for whoever maintains the code, but never reach the image. Multi-line
#      subtitles are truncated to their first line. The explanatory text lives in FIGURE_LEGENDS.md.
#   2. TEXT IS BLACK AND LARGE. Base size 15; axis text, titles, strips and legends all use primary
#      ink. The earlier grey-on-off-white secondary ink was too faint once figures were scaled down
#      into slides.
#   3. Colours are unchanged -- the validated categorical order, sequential for magnitude, diverging
#      only for signed quantities.

PAL <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100",
         magenta = "#e87ba4", green = "#008300", violet = "#4a3aa7", red = "#e34948")
INK  <- c(primary = "#0b0b0b", secondary = "#0b0b0b", muted = "#3d3d3d")
SURF <- "#ffffff"
GRID <- "#dcdbd6"

theme_lcc <- function(base = 15) {
  ggplot2::theme_minimal(base_size = base) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = SURF, colour = NA),
      panel.background = ggplot2::element_rect(fill = SURF, colour = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = GRID, linewidth = 0.3),
      axis.text   = ggplot2::element_text(colour = INK[["primary"]], size = base - 1),
      axis.title  = ggplot2::element_text(colour = INK[["primary"]], size = base),
      plot.title  = ggplot2::element_text(colour = INK[["primary"]], face = "bold", size = base + 3),
      plot.subtitle = ggplot2::element_text(colour = INK[["primary"]], size = base - 1),
      plot.caption  = ggplot2::element_blank(),          # rule 1: never render captions
      strip.text  = ggplot2::element_text(colour = INK[["primary"]], face = "bold", size = base),
      legend.title = ggplot2::element_text(colour = INK[["primary"]], size = base - 1),
      legend.text  = ggplot2::element_text(colour = INK[["primary"]], size = base - 1),
      legend.position = "top")
}

# strip the caption and any explanatory tail of the subtitle before writing
.strip_text <- function(p) {
  # patchwork objects are NOT ggplots: an earlier version returned early here, so composed figures
  # (G1, P5, P6) kept their captions and long subtitles while every single-panel figure lost them.
  if (inherits(p, "patchwork"))
    return(p & ggplot2::theme(plot.caption = ggplot2::element_blank()))
  if (!inherits(p, "ggplot")) return(p)
  p <- p + ggplot2::labs(caption = NULL)
  lab <- tryCatch(p$labels$subtitle, error = function(e) NULL)
  if (is.null(lab)) lab <- tryCatch(p@labels$subtitle, error = function(e) NULL)
  if (!is.null(lab) && is.character(lab) && length(lab) == 1L) {
    lab <- sub("\\s*\n.*$", "", lab)                 # keep only the first line
    # anything still long is an explanation, not a label: it belongs in FIGURE_LEGENDS.md
    p <- p + ggplot2::labs(subtitle = if (nchar(lab) > 72) NULL else lab)
  }
  p
}

# Canvas scale. The base font went from 11 to 15 without the canvases growing, which pushed strip
# labels and tick labels into each other. Scaling every figure by 1.25 restores the room; text is
# still ~10% larger relative to the plot than before, which was the point of the size increase.
FIG_SCALE <- 1.25

save_fig <- function(p, name, w, h) {
  p <- .strip_text(p)
  w <- w * FIG_SCALE; h <- h * FIG_SCALE
  for (e in c("png", "pdf"))
    ggplot2::ggsave(file.path(LCC_FIG_DIR, paste0(name, ".", e)), p,
                    width = w, height = h, dpi = 200, bg = SURF)
  message("    wrote ", name)
}
