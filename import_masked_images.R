# !/usr/bin/env Rscript
#$ -S /usr/local/packages/R-3.5.1/bin/Rscript
#$ -l h_rt=01:00:00
#$ -t 1:817
#$ -l h_vmem=4G
#$ -cwd
#$ -o $HOME/testlogs
#$ -e $HOME/testlogs

print("Start projection")

library(pacman)
p_load(stringr, tidyverse, magrittr)
p_load(fda, oro.nifti, spam)


args <- R.utils::commandArgs(asValue=TRUE, excludeReserved=TRUE)


mask_fname <- "<path-to-file>/smoothed_mask_125"
JobId <- as.numeric(Sys.getenv("SGE_TASK_ID"))
Files <- scan("imagelist.txt", what = "character")
img_PTID <- paste0(".*","<path-to-file>", "|", "_sc", ".*")   %>%        
  gsub(., "", Files[JobId])

print(img_PTID)


load(as.character(args$ws))

resize_image <- function(mask, img = mask){
  mask <- drop(mask)
  img <- drop(img)
  dims <- dim(mask)
  
  if(!all.equal(dim(img), dims)) stop("Mask and img must have the same dimensions!")
  
  nonzero_coord <- matrix(NA, nrow = length(dims), ncol = 2)
  rownames(nonzero_coord) <- paste0("Dim", 1:length(dims))
  colnames(nonzero_coord) <- c("first","last")
  
  
  for(i in 1:length(dims)) nonzero_coord[i,] <- range(which(apply(mask,i,sum) != 0))
  ###select first and last coordinate for which the sum is non zero
  
  resized_array <- (img[,,]*(mask[,,]>0))[nonzero_coord[1,"first"]:nonzero_coord[1,"last"],
                                          nonzero_coord[2,"first"]:nonzero_coord[2,"last"],
                                          nonzero_coord[3,"first"]:nonzero_coord[3,"last"]]
  
  return(list(array = resized_array, original_coord = nonzero_coord))
}




get_image_vector <- function(img_path, mask, voxel_grid_nonzero_mask, dimscan = dim(mask)){
  require(pacman)
  p_load(magrittr, fda, oro.nifti, spam, glmnet, Matrix)
  
  
  array(readBin(img_path, what = "int", n = prod(dimscan), size = 2, signed = FALSE), dim = dimscan) %>%
    ###oro.nifti::readNIfTI(img_path) %>%
    resize_image(mask = mask, img = .) %>%
    .$array %>%
    as.vector() %>%
    .[voxel_grid_nonzero_mask]
}


print("Functions loaded")

mask <- oro.nifti::readNIfTI(mask_fname)

print("Mask loaded")

result <- get_image_vector(Files[JobId], mask, voxel_grid_nonzero_mask) %>%
  data.frame(.) #%>%  round(., 8)

names(result) <- img_PTID

print("Image loaded")

cat("The name is ", names(result))

write.table(t(result), file = paste0("TBM_masked_fast/",img_PTID,".dat"), row.names = TRUE, col.names = FALSE)
