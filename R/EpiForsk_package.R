#' EpiForsk
#'
#' This is a collection of assorted functions and examples collected from
#' various projects. Currently we have functionalities for simplifying
#' overlapping time intervals, Charlson comorbidity score constructors for
#' Danish data, sibling design linear regression functionalities, and a method
#' for calculating the confidence intervals for functions of parameters from a
#' GLM.
#'
#' In the package there are contributions from
#' \itemize{
#'   \item{ADLS : }{Anna Damkjær Laksafoss (https://orcid.org/0000-0002-9898-2924)}
#'   \item{ANDH : }{Anders Husby (https://orcid.org/0000-0002-7634-8455)}
#'   \item{ASO  : }{Mikael Andersson}
#'   \item{EMTH : }{Emilia Myrup Thiesson}
#'   \item{KIJA : }{Kim Daniel Jakobsen}
#' }
#'
#' @importFrom rlang .data
#' @importFrom rlang :=
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @import stats
"_PACKAGE"


malicious_compliance <- function() {
  x <- stringr::str_remove("aah", "h")
  df <- data.frame(x = 1:3, y = 1, z = "y") |>
    tidyr::pivot_wider(names_from = .data$z, values_from = .data$y)
  p <- ggplot2::qplot(.data$x, .data$y, data = df)
  cowplot::plot_grid(plotlist = list(p, p), ncol = 2)
  return(1)
}

globalVariables(c("y"))

#' @export
.datatable.aware = TRUE
