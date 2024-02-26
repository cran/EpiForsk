#' Wrapper for the `odds_ratio_function()`to perform several similar analyses
#' in one go
#'
#' The function is intended to make it easy to get OR's for several similar
#' models in one go, where either the same analysis is performed except for one
#' variable or the same analysis is performed but by each variable (each level
#' of the variable is analysed separately).
#'
#' @param normaldata A data frame or data frame extension (e.g. a tibble).
#' @param outcomevar A character vector naming factor variables in normaldata
#'   to use as outcomes in separate models.
#' @param expvars A character vector naming exposure variables (either numeric
#'   or factors) to use in separate models.
#' @param adjustment_fixed A character vector naming adjustment variables to
#'   include in all models. NULL is the default resulting in no fixed
#'   adjustment.
#' @param by_var A character vector specifying a factor on which to run the
#'   analyses completely separate for all levels. It only works with one
#'   variable (default is NULL). NOTE: NA and "" levels will not be used but all
#'   other levels will have separate models.
#' @param number_decimals An integer giving the number of decimals to show in
#'   the standardized output (default is two decimals).
#' @param alpha A scalar, between 0 and 1 specifying the desired significance
#'   level of the confidence intervals (default is 0.05 which will yield the
#'   usual 95% confidence interval).
#' @param regtype A character string specifying the analysis method. Can either
#'   be "logistic" for logistic regression (the default) or "log-linear" for
#'   log-linear regression. Log-linear regression can only be used with
#'   binomial, unconditional analysis.
#' @param matchgroup Character string specifying a variable in normaldata to
#'   condition the analysis on. Can only be used in binomial logistic regression
#'   models (default is NULL).
#' @param matchtiemethod Character string specifying the method for ties when
#'   using a matched/conditional analysis. The default options is "exact",
#'   however this option does not take weights into account for the analysis, so
#'   if weights (other than 1) are used, another option should be selected.
#'   Other options are "approximate", "efron", and "breslow" - for further
#'   explanations, see documentation for \link[survival]{clogit}.
#' @param values_to_remove Character vector specifying values to remove from
#'   ALL variables used in the regression before the analysis (default is NULL).
#'   This is useful if some value(s) are used consistently to encode
#'   missing/irrelevant in the data (e.g. c("888", "987") - normal missing (NA)
#'   don't need to be specified as it will be removed automatically. Do NOT
#'   remove the reference values as this will lead to unexpected results!
#' @param weightvar A character string specifying a numeric variable in
#'   normaldata with pre-calculated weights for observations in the analysis.
#'   The default value NULL corresponds to weight 1 for all observations.
#' @param surveydata A Boolean specifying whether the data comes from a survey
#'   (default is FALSE).
#' @param textvar A character string with text (like a note) to be added to the
#'   output. The default value NULL corresponds to no added note.
#' @param model_object A Boolean. If TRUE, returns the raw output object from
#'   the analysis instead of the standard output. This might be useful to see
#'   information not included in the standardized output (default is FALSE).
#'
#' @details It's possible to have same variable in `expvars` and
#' `adjustment_fixed`.
#'
#' When a model results in an error, the function will not stop - it continues
#' with other models until done BUT in the output the error text can be seen.
#'
#' @return A standardized analysis object with results from multiple models.
#'
#' @author ASO
#'
#' @seealso [odds_ratio_function()] to perform a single logistic or log-linear
#' regression giving a standardized output table.
#'
#' @examples
#' # Data to use
#' data("infert", package = "datasets")
#' infert2 <- infert |>
#'   dplyr::mutate(
#'     Age_grp = relevel(as.factor(dplyr::case_when(
#'       age < 25 ~ "<25",
#'       25 <= age & age < 35 ~ "25-<35",
#'       age >= 35 ~ "35+"
#'     )), ref="25-<35"),
#'     Parity_grp = relevel(as.factor(dplyr::case_when(
#'       parity == 1 ~ "1",
#'       parity >= 2 & parity <= 3 ~ "2-3",
#'       parity > 3 ~ "4+"
#'     )), ref="2-3"),
#'     induced = relevel(as.factor(induced), ref="0"),
#'     case = relevel(as.factor(case), ref="0"),
#'     spontaneous = relevel(as.factor(spontaneous), ref="0")
#'   )
#'
#' # Two outcomes (Parity_grp, case) with their own set of models, three
#' # variables included in separate models (spontaneous,induced and education)
#' # and one variable that is included in all models (Age_grp)
#' test <- odds_ratio_function_repeated(
#'   normaldata = infert2,
#'   outcomevar = c("Parity_grp","case"),
#'   expvars = c("spontaneous","induced","education"),
#'   adjustment_fixed = c("Age_grp")
#' )
#'
#' # One outcome (case), two variables included in separate models
#' # (spontaneous and induced), one variable included in all models (Age_grp)
#' # and all analyses made for each level of another variable (Parity_grp)
#' test2 <- odds_ratio_function_repeated(
#'   normaldata = infert2,
#'   outcomevar = c("case"),
#'   expvars = c("spontaneous","induced"),
#'   adjustment_fixed = c("Age_grp"),
#'   by_var = "Parity_grp"
#' )
#'
#' @export

odds_ratio_function_repeated <- function(
    normaldata,
    outcomevar,
    expvars,
    adjustment_fixed = NULL,
    by_var = NULL,
    number_decimals = 2,
    alpha = 0.05,
    regtype = c("logistic", "log-linear"),
    matchgroup = NULL,
    matchtiemethod = c("exact", "approximate", "efron", "breslow"),
    values_to_remove = NULL,
    weightvar = NULL,
    surveydata = FALSE,
    textvar = NULL,
    model_object = FALSE
) {
  regtype <- match.arg(regtype)
  matchtiemethod <- match.arg(matchtiemethod)
  new_expvars_prp <- expvars
  new_expvars_prp2 <- unlist(strsplit(new_expvars_prp, ":"))
  new_expvars <- unlist(strsplit(new_expvars_prp2, "[*]"))
  func_var_names <- unique(
    c(outcomevar, new_expvars, adjustment_fixed, by_var, weightvar)
  )

  func_table1 <- normaldata |>
    dplyr::select(dplyr::all_of(func_var_names))

  if (is.null(by_var)) {
    by_var_level_count <- 1
  } else {
    by_var2 <- dplyr::pull(
      unique(dplyr::distinct(dplyr::select(func_table1, {{ by_var }})))
    )
    if (is.numeric(by_var2)) {
      by_var3 <- by_var2[!is.na(by_var2)]
    }
    if (is.factor(by_var2) | is.character(by_var2)) {
      by_var3 <- as.character(by_var2[!is.na(by_var2) & by_var2 != ""])
    }
    by_var_level_count <- length(by_var3)
  }

  outcome_var_count <- length(outcomevar)
  expvars_var_count <- length(expvars)

  for (k in seq_len(by_var_level_count)) {
    if (is.null(by_var)) {
      By_var_name <- c("None")
      func_table1_2 <- func_table1
    } else {
      by_var_level <- dplyr::nth(by_var3, n = k)
      By_var_name <- paste(by_var, "=", by_var_level)
      print(paste0("By_var: ", By_var_name))
      func_table1_2 <- dplyr::filter(
        func_table1,
        as.character(.data[[by_var]]) == by_var_level
      )
    }

    for (i in seq_len(outcome_var_count)) {
      Outcome_var_name <- dplyr::nth(outcomevar, n = i)
      print(paste0("Outcome: ", Outcome_var_name))

      for (j in seq_len(expvars_var_count)) {
        Expvar_var_name <- dplyr::nth(expvars, n = j)
        print(paste0("Expvar: ", Expvar_var_name))
        new_expvar <- unique(c(Expvar_var_name, adjustment_fixed))
        func_res1 <- try_catch_warnings(
          odds_ratio_function(
            normaldata = func_table1_2,
            outcomevar = Outcome_var_name,
            expvars = new_expvar,
            number_decimals = number_decimals,
            alpha = alpha,
            regtype = regtype,
            matchgroup = matchgroup,
            matchtiemethod = matchtiemethod,
            values_to_remove = values_to_remove,
            weightvar = weightvar,
            surveydata = surveydata,
            textvar = textvar,
            model_object = model_object
          ),
          character = TRUE
        )
        if (model_object == FALSE & func_res1$error == '') {
          func_res2_prp <- func_res1$value |>
            dplyr::mutate(
              By_name = By_var_name,
              Outcome_name = Outcome_var_name,
              Expvar_name = Expvar_var_name
            )
          if (func_res1$warning != '') {
            func_res2 <- func_res2_prp |>
              dplyr::mutate(
                Warning = dplyr::case_when(
                  .data$term == "(Intercept)" ~ paste0(func_res1$warning),
                  TRUE ~ ""
                )
              )
          } else {
            func_res2 <- func_res2_prp
          }
          if (k==1 & i==1 & j==1) {
            func_table2 <- func_res2
          } else {
            func_table2 <- dplyr::bind_rows(func_table2, func_res2)
          }
        } else if (model_object == TRUE & func_res1$error == '') {
          if (k==1 & i==1 & j==1) {
            func_table2 <- c(
              "By" = By_var_name,
              "Outcome" = Outcome_var_name,
              "Expvar" = Expvar_var_name,
              "Warning" = func_res1$warning,
              func_res1$value
            )
          } else {
            func_table2 <- c(
              func_table2,
              "By" = By_var_name,
              "Outcome" = Outcome_var_name,
              "Expvar" = Expvar_var_name,
              "Warning" = func_res1$warning,
              func_res1$value
            )
          }
        } else if (model_object == FALSE & func_res1$error != '') {
          func_res2 <- as.data.frame(func_res1$error) |>
            dplyr::rename(Error = 1) |>
            dplyr::mutate(
              By_name = By_var_name,
              Outcome_name = Outcome_var_name,
              Expvar_name = Expvar_var_name
            ) |>
            dplyr::relocate("Error", .after = dplyr::last_col())
          if (k == 1 & i == 1 & j == 1) {
            func_table2 <- func_res2
          } else {
            func_table2 <- dplyr::bind_rows(func_table2, func_res2)
          }
        } else if (model_object == TRUE & func_res1$error != '') {
          if (k == 1 & i == 1 & j == 1) {
            func_table2 <- c(
              "By" = By_var_name,
              "Outcome" = Outcome_var_name,
              "Expvar" = Expvar_var_name,
              "Error" = func_res1$error
            )
          } else {
            func_table2 <- c(
              func_table2,
              "By" = By_var_name,
              "Outcome" = Outcome_var_name,
              "Expvar" = Expvar_var_name,
              "Error" = func_res1$error
            )
          }
        }
      }
    }
  }

  if (model_object == FALSE) {
    func_table3 <- func_table2 |>
      dplyr::relocate("Expvar_name", .before = 1) |>
      dplyr::relocate("Outcome_name", .before = 1) |>
      dplyr::relocate("By_name", .before = 1)
  } else{
    func_table3 <- func_table2
  }
  return(func_table3)
}
