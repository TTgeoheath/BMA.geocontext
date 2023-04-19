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
lac.limit <- function(data,box_widths, lacunarity_values, start = list(a = 1, b = -1, c = 1)) {
  if (!is.data.frame(data))
   stop("data must be a data.frame!")
  if (!is.numeric(box_widths))
  stop("'box_widths' must be a numeric vector!")
  if (!is.numeric(lacunarity_values))
    stop("'lacunarity_values' must be numeric vector!")
  if (!is.list(start))
    stop("start must be a list!")
  
  # Step 1
  ## Fit the data to a power function
  power_model <- nlsLM(lacunarity_values ~ a * box_widths^b + c, data = data, start = start)
  
  ## Extract the coefficients
  a <- coef(power_model)["a"]
  b <- coef(power_model)["b"]
  c <- coef(power_model)["c"]
  # Calculate the box size where the maximum curvature appears
  max_curvature <- soilphysics::maxcurv(x.range = range(box_widths), fun = function(x)a* x^b + c, method = "pd")
  
  return(max_curvature)
}
