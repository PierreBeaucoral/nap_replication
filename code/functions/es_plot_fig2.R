# ==============================================================================
# code/functions/es_plot_fig2.R
# Single-curve event-study plot in the convention of Figure 2
# (03_main_results.R): simultaneous bands drawn as vertical bars, points
# without connecting lines, the base period e = -1 plotted at zero, a tick at
# every event time, no in-figure title. Used by stages 04 and 07.
# ==============================================================================

#' Event-study plot in the style of Figure 2
#'
#' @param es data.frame with event_time, ATT, Lower, Upper
#' @param colour colour of bars and points
#' @param x_lab x-axis title
#' @param y_lab y-axis title
#' @return a ggplot object
es_plot_fig2 <- function(es, colour = "#2E86C1",
                         x_lab = "Years relative to NAP adoption",
                         y_lab = "ATT estimate") {
  ggplot2::ggplot(es, ggplot2::aes(x = event_time, y = ATT)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = -0.5, colour = "grey30", linetype = "dotted") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = Lower, ymax = Upper), colour = colour,
                            linewidth = 0.6, alpha = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(colour = colour, shape = 16, size = 2.5, na.rm = TRUE) +
    ggplot2::scale_x_continuous(breaks = seq(min(es$event_time), max(es$event_time))) +
    ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL, x = x_lab, y = y_lab) +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(family = "serif", size = 12),
                   panel.grid.minor = ggplot2::element_blank())
}
