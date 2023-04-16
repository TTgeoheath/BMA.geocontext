#' model.overall
#'
#' This function calculates the averaged coefficients (i.e.,beta, p-value, z-value, standard error, hazard ratio,
#' confidence interval) weighted by posterior possibility
#'
#' @param data A data frame containing the survival time, status,covariates,and exposures delineated by different geographic contexts
#' @param surv_time Te column name of survivial time
#' @param status The column name of status
#' @param covariates A character vector of column names for the fixed covariates
#' @param exposure_colnames A character vector of column names for the exposure variables, the column names of exposures should be 
#' the combination of prefixes and buffer sizes, and the same exposures should have the same prefixes, for example, the column names 
#' of greenspace for buffer sizes 100, 200, 300, 400 could be bufnd100, bufnd200,bufnd300,bufnd400; for noise could be bufno100,
#' bufno200, bufno300, and bufno400.Note, the prefixes can only be alphabets, e.g., bufnd,bufpm,bufno..
#' @return A dataframe includes the the coefficients, posterior possibility, and BIC for each exposure-specific cox model
#' @examples 
#' \dontrun{Take our provided fake mortality data (mortality_exposure.csv) as an example:
#' data <- read.csv("~path/mortality_exposure.csv")
#' head(data)
#' covariates <- colnames(data[,c(3:7)])
#' exposure_colnames <- colnames(data[,c(8:29)])
#' overall_model <- model.overall(data,"Survival_time","Status",covariates,exposure_colnames)
#' }
#' @import 
#' survival
#' MASS
#' dplyr
#' @importFrom stats BIC confint na.omit
#' @export
model.overall<- function(data,surv_time,status,covariates,exposure_colnames) {
  # data: a data frame containing the covariates, survival time, status, and exposur
  # surv_status:  a character vector of column names for survival time and status
  # covariates: a character vector of column names for the fixed covariates
  # exposure_colnames: a character vector of column names for the exposure variables
  # number_exposures: the number of different exposures included in the model
  if (!is.data.frame(data))
    stop("data must be a data.frame!")
   if (!is.character(surv_time))
    stop("'surv_time' must be a charactor!")
  if (!is.character(status))
    stop("'status' must be a charactor!")
  if (!is.character(covariates))
    stop("'covariates' must be a charactor vector!")
  if (!is.character(exposure_colnames))
    stop("'exposure_colnames' must be a charactor vector!")
  # Step 1: Create all possible combinations of models
  surv_time <-data[,c(surv_time)]
  status <- data[,c(status)]
  covariates <- data[,c(covariates)]
  ##split the exposure colnames based on the prefixes
  prefixes <- unique(sub("\\d+", "", exposure_colnames))
  split_by_prefix <- function(prefix) {
    return(grep(prefix, exposure_colnames, value = TRUE))
  }
  # Apply the function to each prefix and create a named list of vectors
  split_exposure_colnames_list <- lapply(prefixes, split_by_prefix)
  # Generate all possible combinations of buffer sizes
  combinations <- do.call(expand.grid, split_exposure_colnames_list)
  ## Transfer the factor datatype to characters
  combinations[] <- lapply(combinations, function(x) {
    if (is.factor(x)) {
      return(as.character(x))
    } else {
      return(x)
    }
  })
  
  ## transfer it to a list
  data_list <- lapply(1:nrow(combinations), function(i) as.list(combinations[i, ])) 
  ## Unnamed each elements in the list
  model_combinations <- lapply(model_combinations, function(x) {
    unname(x)
  })
  ##Transfer factors to characters
  model_combinations<- lapply(model_combinations, as.character) 
  # Step 3: Fit Cox regression models for each combination of variables and save the models
  cox_models <- lapply(model_combinations, function(x) {
    # Get the columns corresponding to the current model combination
    dataset <- data[,x]
    dataset <- cbind(dataset,surv_time,status, covariates)
    # Fit the Cox regression model
    cox_models <- coxph(Surv(surv_time,status) ~ ., data = dataset,na.action=na.omit)
    #cox_models <- coxph(Surv(Survival_time,Status) ~ ., data = cox_models,na.action=na.omit)
    # Return the fitted model
    return(cox_models)
  })
  
  # Step 4: Calculate the BIC for each model
  model_bics <- sapply(cox_models, function(x) {
    # Calculate the BIC for the model
    bic <- BIC(x)
    
    # Return the BIC value
    return(bic)
  })
  
  # Step 5: Calculate the posterior probabilities, hazard ratios, and confidence intervals for each model
  likelihoods <- exp(-0.5 * (model_bics - mean(model_bics)))
  prior <- 1 / length(model_combinations)
  posterior_probs <- likelihoods * prior / sum(likelihoods * prior)
  
  # Get the hazard ratios and confidence intervals for each model
  ci <- sapply(cox_models, function(x) {
    # Get theconfidence intervals for the exposure variables
    ci <-  exp(confint(x, level = 0.95))
    # Return the confidence intervals
    return(ci)
  })
  ##extract the higher ci and lower ci for exposures.
  lower_ci <- ci[1:length(prefixes),]
  lower_ci <- t(lower_ci)
  lowerci_weighted <- lower_ci*posterior_probs
  overall_lowerci <- colSums(lowerci_weighted)
  overall_lowerci <-as.data.frame(matrix(overall_lowerci, nrow = 1, ncol = length(overall_lowerci), byrow = TRUE))
  low_ci_names <- paste0(prefixes, c("_lowerci_overall"))
  count_variables <- ncol(covariates)+length(prefixes)
  higher_ci <- ci[(1+count_variables):(count_variables+length(prefixes)),]
  higher_ci <- t(higher_ci)
  higherci_weighted <- higher_ci*posterior_probs
  overall_higherci <- colSums(higherci_weighted)
  overall_higherci <-as.data.frame(matrix(overall_higherci, nrow = 1, ncol = length(overall_higherci), byrow = TRUE))
  high_ci_names <- paste0(prefixes, c("_higherci_overall"))
  ci_names <- c(low_ci_names,high_ci_names)
  ci_overall <- cbind(overall_lowerci,overall_higherci)
  colnames(ci_overall) <- ci_names
  ## extract coefficients(i.e., beta,hazard ratio, standard error, z, and p value) 
  coef<- sapply(cox_models, function(x) {
    # Get the coefficients
    coef <- coef(summary(x)) 
    return(coef)
  })
  ## twist the matrix
  coef<- t(coef)
  ## get the beta for exposures
  Beta <- coef[,1:length(prefixes)]
  ## get the Hazard Ratio for exposures
  HR <- coef[,(1+count_variables):(count_variables+length(prefixes))]
  ## get the standard error for exposures for each model
  SE <- coef[,(2*count_variables+1):(2*count_variables+length(prefixes))]
  ## get the z_value for exposures for each model
  z <- coef[,(3*count_variables+1):(3*count_variables+length(prefixes))]
  ## get the P_value for exposures for each model
  P_value <- coef[,(4*count_variables+1):(4*count_variables+length(prefixes))]
  coef_names <- c(paste0(prefixes, c("_Beta")),paste0(prefixes, c("_HR")),
                  paste0(prefixes, c("_SE")),paste0(prefixes, c("_z")),
                  paste0(prefixes, c("_P_value")))  
  coef <- cbind(Beta,HR,SE,z,P_value)
  coef_weighted <- coef*posterior_probs
  overall_coef <- colSums(coef_weighted)
  overall_coef <-as.data.frame(matrix(overall_coef, nrow = 1, ncol = length(overall_coef), byrow = TRUE))
  colnames(overall_coef) <- coef_names 
  # Step 6: Create a data frame with the model combinations, posterior probabilities, hazard ratios, and confidence intervals
  result <- t(cbind(overall_coef,ci_overall))
  
  return(result)
}
