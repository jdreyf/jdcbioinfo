#' Wrapper for print out PNG
#'
#' Print out PNG in in HTML format. This is useful for R Markdown reports.
#'
#' @param path The vector of PNG file names (<=3 files)
#' @param title The vector of figure titles
#' @param caption The vector of figure captions
#' @param one_col_width The width of one-column figure
#' @return NULL
#' @export

show_figure <- function(path, title, caption, one_col_width = c("50%", "80%")) {
  stopifnot(length(path) == length(title), length(path) == length(caption))
  stopifnot(length(path) < 4, length(title) < 4, length(caption) < 4)
  stopifnot(one_col_width %in% c("50%", "80%"))

  path <- paste0(path, ".png")
  stopifnot((file.exists(path)))

  ncol = length(path)
  esc <- function(x) as.character(htmltools::htmlEscape(x, attribute = TRUE))

  if (ncol == 1) {
    one_col_width <- match.arg(one_col_width,
                               choices = c("50%", "80%"),
                               several.ok = FALSE)

    if (one_col_width == "50%") {
      cat(sprintf(paste0(
        '<div class="one-column-center-container">\n',
        '<div class="one-column-center plot-card">\n',
        '<a href="%s"><img src="%s" alt="%s" style="width:100%%;"></a>\n',
        '<div class="figure-caption"><strong>%s.</strong> %s</div>\n',
        '</div>\n</div>\n<div class="fix-clear"></div>\n\n',
        '<p>Download: <a href="%s">%s</a></p>\n\n'),
        esc(path), esc(path), esc(title), esc(title), esc(caption),
        esc(path), esc(basename(path))))
    } else if (one_col_width == "80%") {
      cat(sprintf(paste0(
        '<div class="one-column-center-container">\n',
        '<div class="one-column-center-80 plot-card">\n',
        '<a href="%s"><img src="%s" alt="%s" style="width:100%%;"></a>\n',
        '<div class="figure-caption"><strong>%s.</strong> %s</div>\n',
        '</div>\n</div>\n<div class="fix-clear"></div>\n\n',
        '<p>Download: <a href="%s">%s</a></p>\n\n'),
        esc(path), esc(path), esc(title), esc(title), esc(caption),
        esc(path), esc(basename(path))))
    } else {
      stop("Invalid colum width. Use either 50% or 80%")
    }

  } else if (ncol == 2) {
    cat(sprintf(paste0(
      '<div class="two-column-container">\n',
      '<div class="two-column-left plot-card">\n',
      '<a href="%s"><img src="%s" alt="%s" style="width:100%%;"></a>\n',
      '<div class="figure-caption"><strong>%s.</strong> %s</div>\n',
      '</div>\n',
      esc(path[1]), esc(path[1]), esc(title[1]), esc(title[1]), esc(caption[1]),
      '<div class="two-column-right plot-card">\n',
      '<a href="%s"><img src="%s" alt="%s" style="width:100%%;"></a>\n',
      '<div class="figure-caption"><strong>%s.</strong> %s</div>\n',
      '</div>\n',
      esc(path[2]), esc(path[2]), esc(title[2]), esc(title[2]), esc(caption[2]),
      '</div>\n</div>\n<div class="fix-clear"></div>\n\n',
      '<p>Download: <a href="%s">%s</a> & <a href="%s">%s</a></p>\n\n'),
      esc(path[1]), esc(basename(path[1])), esc(path[2]), esc(basename(path[2]))
      ))
  } else if (ncol == 3) {
    cat(sprintf(paste0(
      '<div class="three-column-container">\n',
      '<div class="three-column-left plot-card">\n',
      '<a href="%s"><img src="%s" alt="%s" style="width:100%%;"></a>\n',
      '<div class="figure-caption"><strong>%s.</strong> %s</div>\n',
      '</div>\n',
      esc(path[1]), esc(path[1]), esc(title[1]), esc(title[1]), esc(caption[1]),
      '<div class="three-column-center plot-card">\n',
      '<a href="%s"><img src="%s" alt="%s" style="width:100%%;"></a>\n',
      '<div class="figure-caption"><strong>%s.</strong> %s</div>\n',
      '</div>\n',
      esc(path[2]), esc(path[2]), esc(title[2]), esc(title[2]), esc(caption[2]),
      '</div>\n',
      '<div class="three-column-right plot-card">\n',
      '<a href="%s"><img src="%s" alt="%s" style="width:100%%;"></a>\n',
      '<div class="figure-caption"><strong>%s.</strong> %s</div>\n',
      '</div>\n',
      esc(path[3]), esc(path[3]), esc(title[3]), esc(title[3]), esc(caption[3]),
      '</div>\n</div>\n<div class="fix-clear"></div>\n\n',
      '<p>Download: <a href="%s">%s</a> & <a href="%s">%s</a> & <a href="%s">%s</a></p>\n\n'),
      esc(path[1]), esc(basename(path[1])), esc(path[2]), esc(basename(path[2])),
      esc(path[3]), esc(basename(path[3]))
    ))
  } else {
    stop("Too many figures. Use 1 to 3 figures only")
  }


}
