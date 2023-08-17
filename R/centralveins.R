#' @title Central Vein Detection
#' @description This function obtains the probability that each lesion in a subject's deep white-matter has a central vein.
#' @param epi a T2*-EPI volume of class \code{nifti}.
#' @param t1 a T1-weighted volume of class \code{nifti}.
#' @param flair a T2-FLAIR volume of class \code{nifti}.
#' @param probmap an image of class \code{nifti}, containing the probability that each voxel
#' is a lesion voxel.
#' If a probability map is not included, the MIMoSA model will be applied (Valcarcel et al., 2017).
#' @param binmap a \code{nifti} mask in which voxels are classified as either lesion voxels
#' or not lesion voxels.
#' Note that mask should be in the same space as the probmap volume.
#' @param parallel is a logical value that indicates whether the user's computer
#' is Linux or Unix (i.e. macOS), and should run the code in parallel.
#' @param cores if parallel = TRUE, cores is an integer value that indicates how many cores
#' the function should be run on.
#' @param skullstripped a logical value reflecting whether or not the images have already been skull-stripped.
#' @param biascorrected a logical value reflecting whether or not the images have already been bias-corrected.
#' @param c3d a logical value reflecting whether or not the Convert3D imaging toolbox is installed.
#'
#' @importFrom ANTsRCore labelClusters
#' @importFrom neurobase niftiarr
#' @import mimosa
#' @importFrom extrantsr ants2oro oro2ants bias_correct registration fslbet_robust
#' @importFrom stats predict
#' @importFrom fslr fslsmooth fast fslerode
#' @return A list containing candidate.lesions (a nifti file with labeled lesions evaluated for CVS),
#' cvs.probmap (a nifti file in which candidate lesions are labeled with their CVS probability), and
#' cvs.biomarker (a numeric value representing the average CVS probability of a subject's lesions).
#' @examples \dontrun{
#' library(neurobase)
#' epi <- readnii("path/to/epi")
#' flair <- readnii("path/to/flair")
#' t1 <- readnii("path/to/t1")
#' cvs <- centralveins(
#'   epi = epi, t1 = t1, flair = flair,
#'   parallel = TRUE, cores = 4, c3d = T
#' )
#' }
#' @export
centralveins <- function(epi, t1, flair, mask, 
                         probmap = NULL, binmap = NULL, 
                         parallel = F, cores = 2, 
                         skullstripped = F, biascorrected = F, registered = F, 
                         c3d = F) {
  if (biascorrected == F) {
    epi <- bias_correct(epi, correction = "N4", reorient = F)
    t1 <- bias_correct(t1, correction = "N4", reorient = F)
    flair <- bias_correct(flair, correction = "N4", reorient = F)
  }
  if (registered == F) {
    flair <- registration(
      filename = flair, template.file = t1, typeofTransform = "Rigid",
      remove.warp = FALSE, outprefix = "fun"
    )
    flair <- flair$outfile
  }
  
  if (skullstripped == F) {
    t1_ss <- fslbet_robust(t1, correct = F)
    epi_ss <- fslbet_robust(epi, correct = F)
    flair_ss <- flair
    flair_ss[t1_ss == 0] <- 0
  } else {
    t1_ss <- t1
    epi_ss <- epi
    flair_ss <- flair
  }
  
  frangi <- frangi(image = epi_ss, mask = epi_ss != 0, 
                   parallel = parallel, cores = cores, c3d = c3d)
  frangi[frangi < 0] <- 0
  
  if (registered == F) {
    regs <- labelreg(epi, t1, frangi)
    epi_t1 <- regs$imagereg
    frangi_t1 <- regs$labelreg
  } else {
    epi_t1 <- epi
    frangi_t1 <- frangi
  }
  
  if (is.null(probmap)) {
    mimosa_data <- mimosa_data(
      brain_mask = mask,
      FLAIR = flair_ss,
      T1 = t1_ss,
      normalize = "Z",
      cores = cores,
      verbose = TRUE
    )
    
    mimosa_df <- mimosa_data$mimosa_dataframe
    mimosa_cm <- mimosa_data$top_voxels
    rm(mimosa_data)
    
    predictions <- predict(mimosa::mimosa_model_No_PD_T2,
                           newdata = mimosa_df, type = "response"
    )
    probmap <- niftiarr(mask, 0)
    probmap[mimosa_cm == 1] <- predictions
    probmap <- fslsmooth(probmap,
                         sigma = 1.25, mask = mask,
                         retimg = TRUE, smooth_mask = TRUE
    )
  }
  
  if (is.null(binmap)) {
    binmap <- probmap
    binmap[probmap >= 0.2] <- 1
    binmap[probmap < 0.2] <- 0
  }
  
  les <- lesioncenters(probmap, binmap, 
                       parallel = parallel, 
                       cores = cores, c3d = c3d)
  
  #csf <- fast(t1_orig, opts = "--nobias") # Doesn't work in my testing (FH)
  csf <- fuzzySpatialCMeansSegmentation(oro2ants(t1), mask = oro2ants(mask), 
                                        numberOfClusters = 3)$segmentationImage # Substitute for fast (1 = CSF)
  csf[csf != 1] <- 0
  csf <- ants2oro(labelClusters(csf, minClusterSize = 300))
  csf[csf > 0] <- 1
  csf <- (csf != 1)
  csf <- fslerode(csf, kopts = paste("-kernel boxv", 3), verbose = TRUE)
  csf <- (csf == 0) # End up with mask of large CSF clusters that are dilated
  
  labels <- les$lesioncenters # Already labeled...? FH
  #labels <- ants2oro(labelClusters(oro2ants(les$lesioncenters), minClusterSize = 27))
  
  for (j in 1:max(labels)) {
    if (sum(csf[labels == j]) > 0) {
      labels[labels == j] <- 0
    }
  }
  les <- labels > 0
  dtb <- dtboundary(les)
  
  labels <- ants2oro(labelClusters(oro2ants(les), minClusterSize = 27))
  probles <- labels
  avprob <- NULL
  maxles <- max(labels)
  for (j in 1:maxles) {
    frangsub <- frangi[labels == j]
    centsub <- dtb[labels == j]
    coords <- which(labels == j, arr.ind = T)
    prod <- frangsub * centsub
    score <- sum(prod)
    nullscores <- NULL
    for (k in 1:1000) {
      samp <- sample(1:length(centsub))
      centsamp <- centsub[samp]
      coordsamp <- coords[samp, ]
      sampprod <- frangsub * centsamp
      sampscore <- sum(sampprod)
      nullscores <- c(nullscores, sampscore)
    }
    lesprob <- sum(nullscores < score) / length(nullscores)
    avprob <- c(avprob, lesprob)
    probles[labels == j] <- lesprob
    
    print(paste0("Done with lesion ", j, " of ", maxles))
  }
  
  return(list(candidate.lesions = labels, cvs.probmap = probles, cvs.biomarker = mean(avprob)))
}
