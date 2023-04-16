#' lacunarity
#'
#' This function calculates the lacunarity value for both binary and continuous raster layers.This function 
#' was developed from the function gblemp of the R package lacunaritycovariance (Author: Kassel Liam Hingee), which
#' was extended to apply on continuous raster layers, and the box width will be increased with pixel width.
#'
#' @param init_boxwidth The box width for the beginning box put in the corner of a raster surface
#' @param end_boxwidth The box width for your identified last box
#' @param image The image formate of the raster layer
#' @param obserwin The optional observation window. If provided, it will be the intersection of the specified observation 
#' window (obserwin) and the non-NA pixels in the input image (image).
#' @param raster_type The selection of the raster type. It can be either "continuous" or "binary"
#' @return A dataframe includes box widths (box_width) and their corresponding lacunarity values.
#' @examples 
#' \dontrun{library(raster);library(maptools);library(RcppRoll)
#' study_area <- raster(study_area.tif)
#' im_area <- as.im.Raster(study_area)
#' LA <- gbl(30,10000,im_area,raster_type="continuous")
#' }
#' @import 
#' RcppRoll
#' spatstat
#' spatstat.geom
#' maptools
#' latticeExtra
#' @importFrom  spatstat.core fv fvnames
#' @export
##This function was developed from the function gblemp of the R package lacunaritycovariance (Author: Kassel Liam Hingee), which
## was extended to apply on continuous raster layers, and the box width will be increased with pixel width.
gbl <- function(init_boxwidth,end_boxwidth,image, obserwin = Frame(image), raster_type = c("continuous", "binary")){
  if (!is.im(image)){stop("input image must be of class image")}
  if (abs(image$xstep - image$ystep) > 1E-2 * image$xstep){stop("image pixels must be square")}
  #convert boxs to odd pixel amounts, taking into account that want a distance to edge box
  boxsizes <- seq(init_boxwidth,end_boxwidth,by=image$xstep)
  r <- round(boxsizes/ image$xstep) 
  r <- unique(r)
  distance<- r * image$xstep # the pixel amounts and multiple the width of each pixel 
  
  #compute observation mask
  obsobj <- image
  obsobj[is.finite(image$v)] <- TRUE
  if (class(obserwin) == "im"){obsobj <- eval.im(obserwin * obsobj)}
  obsobj <- as.owin(obsobj) #owin format may not be needed anymore
  if (class(obserwin) == "owin"){obsobj <- intersect.owin(obsobj, obserwin)}
  
  if (requireNamespace("RcppRoll") != TRUE){
    stop("RcppRoll package must be installed to calculate empirical gliding box lacunarity")
  }
  
  image[(complement.owin(intersect.owin(obserwin, Frame(image)), frame = Frame(image)))] <- NA  #make sure the pixels outside obserwin are set to NA so that reduce sampling happens naturally ##NOTE: this a time consuming operation that may never be needed
  if (raster_type == "continuous") {
    LAS <- mapply(gbl_grey, bowpix = r, MoreArgs = list(image = image, obserwin = obsobj), SIMPLIFY = FALSE)
  } else {
    LAS <- mapply(gbl_binary, bowpix = r, MoreArgs = list(image = image, obserwin = obsobj), SIMPLIFY = FALSE)
  }
  
  #converting results in fv objects
  muldf <- matrix(unlist(LAS), ncol = length(LAS[[1]]), byrow = TRUE)
  colnames(muldf) <- names(LAS[[1]])
  LASdf <- cbind(data.frame(box_width = distance), muldf)
  #recommended xlim:
  alim.min <- 1
  alim.max <- min(which(vapply(LASdf[, "lacunarity_values"], is.na, FUN.VALUE = TRUE)), nrow(LASdf))
  LASfv <- fv(LASdf,
           argu = "box_width",
             valu = "lacunarity_values",
             fmla = ".y ~ box_width",
             alim = c(LASdf[alim.min, "box_width"], LASdf[alim.max, "box_width"]),
             ylab = expression(lacunarity_values[gb]),
             unitname = unitname(image),
             labl = c("Box Width",
                       "Lacunarity value"),
              desc = c("the box width of gliding box", 
                    "Lacuanrity vaule corresponding to box widths"
              ),
             fname = "GBL"
  )
  fvnames(LASfv, a = ".") <- "lacunarity_values"
  return(LASfv)
}


##The following function calculates lacunarity for a box with side lengths (multiple pixels). 
##The RS version is automatically calculated by ignoring those boxes that have range or sums that include NA values
## for the continuous
gbl_grey <- function(image, bowpix, obserwin = Frame(image)){
  mat <- as.matrix(image)
  if ( (bowpix > nrow(mat)) | (bowpix > ncol(mat))){
    gbl.con <- NA
    sample_mean <- NA
    sampleran_2nd <- NA
  }
  else {
    glid.overrowsmin <- RcppRoll::roll_min(mat, bowpix)
    glid.overrowsmax <- RcppRoll::roll_max(mat, bowpix)
    glid.overrows<- glid.overrowsmax-glid.overrowsmin## range of the boxwidth switch by rows
    glid.overrowthencolsmin <- RcppRoll::roll_min(t(glid.overrows), bowpix) * image$xstep * image$ystep
    glid.overrowthencolsmax <- RcppRoll::roll_max(t(glid.overrows), bowpix) * image$xstep * image$ystep
    glid.overrowthencols <-glid.overrowthencolsmax- glid.overrowthencolsmin ## the range of the boxes switch by cols after rows
    sample_mean <- mean(glid.overrowthencols, na.rm = TRUE) #sample mean
    sampleran_2nd <- mean(glid.overrowthencols ^ 2, na.rm = TRUE) #biased sample second moment
    gbl.con <- sampleran_2nd / (sample_mean ^ 2) 
  }
  
  return(list(
    lacunarity_values = gbl.con
  ))
}
gbl_binary <- function(image, bowpix, obserwin = Frame(image)) {
  mat <- as.matrix(image)
  if ( (bowpix > nrow(mat)) | (bowpix > ncol(mat))){
    gbl.bin <- NA
    sample_mean <- NA
    sampleran_2nd <- NA
  }
  else {
    glid.overrows <- RcppRoll::roll_sum(mat, bowpix) ## sum of the boxwidth switch by rows
    glid.overthencols <- RcppRoll::roll_sum(t(glid.overrows), bowpix) * image$xstep * image$ystep
    sample_mean <- mean(glid.overthencols, na.rm = TRUE) #sample mean
    sampleran_2nd <- mean(glid.overthencols ^ 2, na.rm = TRUE) #biased sample second moment
    gbl.bin <- sampleran_2nd / (sample_mean ^ 2) 
  }
  
  return(list(
    lacunarity_values = gbl.bin
  ))
}




