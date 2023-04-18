#' lac.limit
#' 
#' The function calculates the precise point where the lacunarity curve flattened out.
#' 
#' @param data A dataframe includes the boxsizes of the gliding boxes and their corresponding lacuanrity value
#' @param box_sizes The increased box_sizes of gliding boxes in lacunarity analysis
#' @param lacunarity_values The corresponding lacunarity values of the increased box sizes given by lacunarity analysis
#' @param start The start value of the curvature simulation. This function used a power function a * x^b + c to simulate 
#' the lacunarity curvature. Through simulation, the flattened point can be obtained by computing the maximum curvature. 
#' A start = list(a = 1, b = -1, c = 1) was set as default. Users can adjusted it based on their data.
#' @return Maximum carvature point
#' @examples  
#' \dontrun{
#' lac_ndvi <- read.csv("~path/lac_ndvi.csv")
#' maxcurv <- lac.limit(lac_ndvi, "distance", "GBL")
#' }
#' @import 
#' soilphysics
#' pracma
#' minpack.lm
#' @importFrom stats coef
#' @export
lac.limit<- function(data, box_sizes, lacunarity_values,start = list(a = 1, b = -1, c = 1)) {
  if (!is.data.frame(data))
    stop("data must be a data.frame!")
  if (!is.character(box_sizes))
    stop("'box_sizes'must be a charactor!")
  if (!is.character(lacunarity_values))
    stop("'lacunarity_values' must be a charactor!")
  if (!is.list(start))
    stop("start must be a list!")
  # Step 1: extract the box_sizes and lacunarity from the dataset
  box_sizes <- data[,c(box_sizes)]
  lacunarity_values <- data[,c(lacunarity_values)]
  # Fit a power function to the data
  power_model <- nlsLM(lacunarity_values ~ a * box_sizes^b+c,data=data,start=start)
  # Create a function for the fitted power curve
  fitted_power_curve <- function(x, a, b, c) {
    return(a * x^b + c)
  }
  # Calculate the box size where the maximum curvature appears
  a <- coef(power_model)["a"]
  b <- coef(power_model)["b"]
  c <- coef(power_model)["c"]
  
  # Use do.call() to pass the parameters to the fitted_power_curve function
  result <- soilphysics::maxcurv(x.range = range(box_sizes), fun = function(x) do.call(fitted_power_curve, c(list(x), coef(power_model))), method = "pd")
  
  return(result)
}
