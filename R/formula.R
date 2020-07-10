#' Parse formula
#' 
#' Parse formula and return dataset containing selected columns. 
#' Interactions are supported for numerical columns only. 
#' An interaction column is the product of all interacting columns.
#'
#' @param formula Object of class \code{formula} or \code{character} describing the model to fit.
#' @param data Training data of class \code{data.frame}.
#' @param env The environment in which the left hand side of \code{formula} is evaluated.
#'
#' @return Dataset including selected columns and interactions.
parse.formula <- function(formula, data, env = parent.frame()) {
  f <- as.formula(formula)
  t <- terms(f, data = data)
  
  ## Get dependent var(s)
  if (is.matrix(data)) {
    response <- data.frame(eval(f[[2]], envir = data.frame(data[, all.vars(f[[2]]), drop = FALSE]), enclos = env))
  } else {
    response <- data.frame(eval(f[[2]], envir = data, enclos = env))
  }
  colnames(response) <- deparse(f[[2]])
  
  ## Get independent vars
  independent_vars <- attr(t, "term.labels")
  interaction_idx <- grepl(":", independent_vars)
  
  ## Error if illegal column name
  if (!all(make.names(independent_vars[!interaction_idx]) == independent_vars[!interaction_idx])) {
    stop("Error: Illegal column names in formula interface. Fix column names or use alternative interface in ranger.")
  }
  
  ## Shortcut if no interactions
  if (all(!interaction_idx)) {
    return(data.frame(response, data[, independent_vars, drop = FALSE], check.names = FALSE))
  }
  
  ## Get interaction columns
  if (any(interaction_idx)) {
    interaction_vars <- independent_vars[interaction_idx]
    dat_interaction <- sapply(strsplit(interaction_vars, ":"), function(x) {
      if (any(!sapply(data[, x, drop = FALSE], is.numeric))) {
        stop("Error: Only numeric columns allowed in interaction terms.")
      }
      with(data, eval(parse(text = paste(x, collapse = "*"))))
    })
    colnames(dat_interaction) <- interaction_vars
  }
  
  ## Get main effect columns
  if (any(!interaction_idx)) {
    main_vars <- independent_vars[!interaction_idx]
    dat_main <- data[, main_vars, drop = FALSE]
  }
  
  ## Return combined data frame
  if (any(!interaction_idx)) {
    data.frame(response, dat_main, dat_interaction, check.names = FALSE)
  } else {
    data.frame(response, dat_interaction, check.names = FALSE)
  }
}