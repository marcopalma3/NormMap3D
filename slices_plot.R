#' Plotting axial slices of brain images
#'
#' @param image.vec A vector containing the values of the images at all the voxel within the mask
#' @param dims A vector with 3 elements for the dimensions
#' @param voxels A vector containing the  of all the voxel within the mask
#' @param col.threshold In a diverging palette, the central value (corresponding to white)
#' @param mask.under The brain mask
#' @param legend.range A vector with the minimum and maximum image value
#' @param z.slice A vector containing the axial slices to be plot (they must be within the third dimension of dims)
#' @param alpha.val Transparency level (between 0 and 1)
#' @param ... Additional parameters from oro.nifti::overlay
#'
#' @return Selected axial slices of the 3D image
#' @export
#'
#' 
slices_plot <- function(image.vec, 
                        dims, 
                        voxels, 
                        col.threshold = 0, 
                        mask.under,
                        legend.range = range(image.vec, na.rm = TRUE), 
                        z.slice = seq(16, 136, by = 5),
                        alpha.val = 0.5,
                        ...){
  img <- rep(NA, prod(dims))
  img[voxels] <- image.vec
  
  ybr <- base::ifelse(diff(legend.range) > 100, 10,0.1) %>% 
    {seq(round(floor(legend.range[1]), log10(1/.)), 
         round(ceiling(legend.range[2]), log10(1/.)), 
         by = .)}
  
  rc1 <- colorRampPalette(colors = c("blue", "white"), space="Lab")(length(which(ybr <= col.threshold)))    
  rc2 <- colorRampPalette(colors = c("white", "red"), space="Lab")(length(which(ybr >= col.threshold)))
  rampcols <- c(rc1, rc2)
  
  
  if(max(ybr) < 1e-3){
    leg.text <- c(format(min(ybr), scientific = TRUE, digits = 3),
                  col.threshold,
                  format(max(ybr), scientific = TRUE, digits = 3))
  } else{
    leg.text <- round(c(min(ybr), col.threshold, max(ybr)), 3)
  }
  
  
  overlay(x = nifti(array(mask.under, dim = dims)), 
          y = nifti(array(img, dim = dims)), 
          plot.type="single", z = z.slice,
          zlim.y = range(ybr), 
          col.y = scales::alpha(rampcols, alpha.val), useRaster = TRUE, oma = c(4,0,0,0), ...)
  fields::image.plot(legend.only=TRUE, zlim = legend.range, 
                     col = scales::alpha(rampcols, alpha.val), horizontal = TRUE, 
                     legend.mar = 1.5, legend.cex=0.5, legend.width = 0.8, 
                     legend.args = list(text = leg.text,
                                        col="white", cex=0.7, side = 1, 
                                        at = c(min(ybr),col.threshold, max(ybr))))
} 
