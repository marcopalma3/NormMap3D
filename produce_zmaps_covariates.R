# !/usr/bin/env Rscript
#$ -S /usr/local/packages/R-3.5.1/bin/Rscript
#$ -l h_rt=04:00:00
#$ -t 1:817
#$ -l h_vmem=4G
#$ -cwd
#$ -o $HOME/testlogs
#$ -e $HOME/testlogs

print("Start quantile map")

library(pacman)
p_load(ggplot2,tidyverse, oro.nifti, Matrix, sn, knitrProgressBar, pipeR, data.table,nortest, diptest, e1071)


args <- R.utils::commandArgs(asValue=TRUE, excludeReserved=TRUE)


JobId <- as.numeric(Sys.getenv("SGE_TASK_ID"))
Files <- scan("masked_imagelist.txt", what = "character")
img_PTID <- paste0("<path-to-file>",  "|", ".dat")   %>%       
  gsub(., "", Files[JobId])

print(img_PTID)
print(Files[JobId])

load(as.character(args$ws))

voxel_active <- which(allgrid_CN$mask > 0)
###quantmap_vec <- rep(NA, length(voxel_active))

subject_image <- Files[JobId] %>% 
  data.table::fread(., drop = 1, header = F)  %>% 
  as.numeric(.)


head(subject_image)
summary(subject_image)

cat("The subject vector is loaded. Number of active voxels:", length(subject_image), "\n")
#cat("The number of active voxels is", length(voxel_active), "\n")

#pb <- progress_estimated(length(voxel_active))


img <- rep(NA, nrow(allgrid_CN))
img[voxel_active] <- subject_image

indgrid <- with(allgrid_CN, which(knots*mask > 0))

quantmap_vec <- rep(NA, length(indgrid))


covariates <- filter(data_ADNI_bl_817, PTID == img_PTID) %>%
    mutate(PTGENDERMale = 1*(PTGENDER == "Male"), AGEminus70 = AGE - 70)  %>% 
    select(AGEminus70, PTGENDERMale) %>%
    as.numeric(.) %>% 
    c(1,.) 

print(covariates)

for(vox in 1:length(indgrid)){
#  update_progress(pb)
  if(vox %% 1000 == 0) cat(vox, "\t")
  quantmap_vec[vox] <- allgrid_CN[indgrid[vox], ] %>% 
    dplyr::mutate(mean = c(X.Intercept.CP., AGEminus70, PTGENDERMale) %*% covariates)  %>%
    dplyr::select(mean, s.d., gamma1)  %>%          ###select(ends_with("radial")) %>%
    as.numeric(.) %>%
    cp2dp(., "SN") %>%
    psn(img[indgrid[vox]], dp = .) %>%
    qnorm(.) %>% 
    round(., 8)
}

thr <- 7
quantmap_vec <- ifelse(abs(quantmap_vec) > thr, sign(quantmap_vec)*thr, quantmap_vec)


cat("The name is", img_PTID, "\n")


load(as.character(args$basis_ws))


mod1 <- allgrid_CN %>% 
  filter(., mask > 0) %>% 
  with(., which(knots == 1)) %>%
  radial_design_mat[., ] %>%
  glmnet::glmnet(.,
                 quantmap_vec,
                 lambda = 0,
                 intercept = TRUE)

cat("The dimensions of the coefficient vector are ", dim(coef(mod1)), ".\n")
cat("R^2 is", mod1$dev.ratio, ".\n")

coefvec <- mod1 %>%
  coefficients(.) %>%
  .[,1] %>%  ###intercept saved
  data.frame(.) %>%
  round(., 8)


names(coefvec) <- img_PTID

write.table(t(coefvec), file = paste0("qmaps_projections_fromgrid/",img_PTID,".dat"), row.names = TRUE, col.names = FALSE)






zmap <- predict(mod1, radial_design_mat) %>%     ###cbind(1, radial_design_mat) %*% coef(mod1) %>%
  as.vector(.)%>%
  data.table(.)


zmap <- ifelse(abs(zmap[[1]])>thr, sign(zmap[[1]])*thr, zmap[[1]]) %>%  data.table(.)

names(zmap) <- img_PTID

summary(zmap)

normstat <- rep(NA, 4)
normstat <- zmap[[1]] %>%
  {c(ad.test(.)$statistic, cvm.test(.)$statistic, dip.test(.)$statistic, e1071::skewness(.))} %>%
  as.numeric(.)

cat(img_PTID, normstat, "\n", sep = " ")


mean_zmap <- mean(zmap[[1]])
median_zmap <- median(zmap[[1]])
diffrange_zmap <- diff(range(zmap[[1]]))
var_zmap <- var(zmap[[1]])
IQR_zmap <- IQR(zmap[[1]])
expos_zmap <- sum(zmap[[1]] >= 3)
exneg_zmap <- sum(zmap[[1]] <= -3)
meanposblock_zmap <- sort(zmap[[1]], decreasing = TRUE) %>% {.[1:floor(length(.)/10000)]} %>% mean
meannegblock_zmap <- sort(zmap[[1]], decreasing = FALSE) %>% {.[1:floor(length(.)/10000)]} %>% mean
meanabsblock_zmap <- sort(abs(zmap[[1]]), decreasing = TRUE) %>% {.[1:floor(length(.)/10000)]} %>% mean
norm_zmap <- sqrt(sum(zmap[[1]]^2))

cat("PTID", "mean_zmap","median_zmap", "diffrange_zmap", "var_zmap", "IQR_zmap", 
    "expos_zmap", "exneg_zmap", "meanposblock_zmap", "meannegblock_zmap", "meanabsblock_zmap", "norm_zmap\n", sep = " ")
cat(img_PTID, c(mean_zmap, median_zmap, diffrange_zmap, var_zmap, IQR_zmap, 
    expos_zmap, exneg_zmap, meanposblock_zmap, meannegblock_zmap, meanabsblock_zmap, norm_zmap))



data.table::fwrite(zmap, file = paste0("zmaps_fromgrid/", img_PTID,"_fromgrid.dat"))
