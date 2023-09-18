#' @title Central Vein Detection
#' @description This function obtains the probability that each lesion in a subject's deep white-matter has a central vein.
#' @param epi T2*-EPI volume of class \code{antsImage}, class \code{nifti}, or path to ".nii.gz" file.
#' @param t1 T1-weighted volume of class \code{antsImage}, class \code{nifti}, or path to ".nii.gz" file.
#' @param flair T2-FLAIR volume of class \code{antsImage}, class \code{nifti}, or path to ".nii.gz" file.
#' @param prob_map volume of class \code{antsImage}, class \code{nifti},  or path to ".nii.gz" file containing the probability that each voxel is a lesion voxel.
#' If a probability map is not included, the MIMoSA model will be applied (Valcarcel et al., 2017).
#' @param bin_map binarized volume of class\code{antsImage}, class \code{nifti}, or path to ".nii.gz" file. Mask in which voxels are classified as either lesion voxels or not lesion voxels.
#' @param verbose logical value for printing pipeline updates (default = FALSE)
#' @param bias_correct a logical value reflecting whether bias correction still needs to be performed (default = TRUE)
#' @param register a logical value reflecting whether images still need to be registered (default = TRUE)
#' @param skull_strip a logical value reflecting whether skull stripping still needs to be performed (default = TRUE)
#'
#' @import ANTsR
#' @import ANTsRCore
#' @importFrom neurobase niftiarr
#' @import mimosa
#' @importFrom extrantsr ants2oro oro2ants fslbet_robust
#' @importFrom stats predict
#' @importFrom fslr fslsmooth fast fslerode
#' @return A list containing candidate.lesions (a nifti file with labeled lesions evaluated for CVS),
#' cvs.prob_map (a nifti file in which candidate lesions are labeled with their CVS probability), and
#' cvs.biomarker (a numeric value representing the average CVS probability of a subject's lesions).
#' @examples \dontrun{
#' library(ANTsRCore)
#' epi <- check_ants("path/to/epi")
#' flair <- check_ants("path/to/flair")
#' t1 <- check_ants("path/to/t1")
#' cvs <- centralveins(epi = epi, t1 = t1, flair = flair)
#' }
#' @export
#'
central_veins <- function(epi, t1, flair,
                          prob_map = NULL, bin_map = NULL,
                          verbose = FALSE,
                          bias_correct = TRUE, register = TRUE, skull_strip = TRUE) {
  # Reading file path or converting nifti to antsImage
  epi <- check_ants(epi)
  t1 <- check_ants(t1)
  flair <- check_ants(flair)

  antsSameMetadata <- function(ants1, ants2) {
    same_direction <- all(antsGetDirection(ants1) == antsGetDirection(ants2))
    same_origin <- all(antsGetOrigin(ants1) == antsGetOrigin(ants2))
    same_spacing <- all(antsGetSpacing(ants1) == antsGetSpacing(ants2))
    return(same_direction & same_origin & same_spacing)
  }

  # Processing prob_map and bin_map
  prob_map_space <- "none"
  if (!is.null(prob_map) & is.null(bin_map)) {
    warning("If prob_map is provided without bin_map. Thresholding prob_map at 0.2")
    bin_map <- prob_map >= 0.3
  }
  if (is.null(prob_map) & !is.null(bin_map)) {
    warning("bin_map cannot be provided without prob_map. prob_map and bin_map will be recalculated using MIMoSA")
    bin_map <- NULL
  }
  if (!is.null(prob_map) & !is.null(bin_map)) {
    if (class(prob_map) == "nifti") {
      prob_map <- oro2ants(prob_map)
    }
    if (class(bin_map) == "nifti") {
      bin_map <- oro2ants(bin_map)
    }
    if (!antsSameMetadata(prob_map, bin_map)) {
      warning("bin_map is not in same space as prob_map. Thresholding prob_map at 0.3")
      bin_map <- prob_map >= 0.3
    }
    if (antsSameMetadata(prob_map, epi)) {
      prob_map_space <- "epi"
    } else if (antsSameMetadata(prob_map, flair)) {
      prob_map_space <- "flair"
    } else if (antsSameMetadata(prob_map, t1)) {
      prob_map_space <- "t1"
    } else {
      warning("Prob_map is not in the same space as t1, flair, or epi. prob_map and bin_map will be recalculated using MIMoSA")
      prob_map <- NULL
      bin_map <- NULL
    }
  }

  # Bias correction
  if (bias_correct == T) {
    if (verbose) {
      print("Performing N4 bias correction")
    }
    epi <- n4BiasFieldCorrection(epi)
    t1 <- n4BiasFieldCorrection(t1)
    flair <- n4BiasFieldCorrection(flair)
  }

  # Registration
  if (register == F & !(antsSameMetadata(t1, flair) & antsSameMetadata(t1, epi))) {
    warning("T1, FLAIR, and EPI images are not in the same space. Re-registering images.")
    register <- T
  }
  if (register == T) {
    if (verbose) {
      print("Registering images")
    }

    t1_registration <- antsRegistration(
      fixed = epi, moving = t1,
      typeofTransform = "Rigid"
    )
    t1 <- antsApplyTransforms(
      fixed = epi, moving = t1,
      transformlist = t1_registration$fwdtransform,
      interpolator = "lanczosWindowedSinc"
    )
    flair_registration <- antsRegistration(
      fixed = epi, moving = flair,
      typeofTransform = "Rigid"
    )
    flair <- antsApplyTransforms(
      fixed = epi, moving = flair,
      transformlist = flair_registration$fwdtransform,
      interpolator = "lanczosWindowedSinc"
    )

    if (prob_map_space == "t1") {
      prob_map <- antsApplyTransforms(
        fixed = epi, moving = prob_map,
        transformlist = t1_registration$fwdtransform,
        interpolator = "lanczosWindowedSinc"
      )
      bin_map <- antsApplyTransforms(
        fixed = epi, moving = bin_map,
        transformlist = t1_registration$fwdtransform,
        interpolator = "nearestNeighbors"
      )
    }
    if (prob_map_space == "flair") {
      prob_map <- antsApplyTransforms(
        fixed = epi, moving = prob_map,
        transformlist = flair_registration$fwdtransform,
        interpolator = "lanczosWindowedSinc"
      )
      bin_map <- antsApplyTransforms(
        fixed = epi, moving = bin_map,
        transformlist = flair_registration$fwdtransform,
        interpolator = "nearestNeighbors"
      )
    }
  }

  # Skull stripping
  if (skull_strip == T) {
    if (verbose) {
      print("Performing FSL robust BET")
    }
    t1_ss <- fslbet_robust(ants2oro(t1), correct = F, verbose = F)
    mask <- t1_ss != 0

    t1 <- t1 * mask
    flair <- flair * mask
    epi <- epi * mask
  } else {
    mask <- ants2oro(t1 != 0)
  }

  csf <- fast(ants2oro(t1), opts = "--nobias", verbose = verbose)
  csf[csf != 1] <- 0
  csf <- ants2oro(labelClusters(csf, minClusterSize = 300))
  csf[csf > 0] <- 1
  csf <- (csf != 1)
  csf <- fslerode(csf, kopts = paste("-kernel boxv", 3), verbose = TRUE)
  csf <- (csf == 0)

  # MIMoSA
  if (is.null(prob_map)) {
    if (verbose) {
      print("Calculating MIMoSA prob_map")
    }
    mimosa_data <- mimosa_data(
      brain_mask = mask,
      FLAIR = ants2oro(flair),
      T1 = ants2oro(t1),
      normalize = "Z",
      verbose = verbose
    )

    predictions <- predict(mimosa::mimosa_model_No_PD_T2,
                           newdata = mimosa_data$mimosa_dataframe,
                           type = "response"
    )
    prob_map <- niftiarr(mask, 0)
    prob_map[mimosa_data$top_voxels == 1] <- predictions
    prob_map <- oro2ants(fslsmooth(prob_map,
                                   sigma = 1.25, mask = mask,
                                   retimg = TRUE, smooth_mask = TRUE
    ))
    bin_map <- prob_map >= 0.3
  }

  if (sum(bin_map) == 0) {
    warning("No lesions detected")
    return(NULL)
  }
  
  les <- label_lesion(prob_map, bin_map)
  labels <- antsImageClone(les)

  if (sum(csf * labels) > 0) {
    for (j in 1:max(labels)) {
      tmp_label_mask <- labels == j
      if (sum(csf * tmp_label_mask) > 0) {
        print(j)
        labels[tmp_label_mask] <- 0
      }
    }
  }

  lesion_mask <- antsImageClone(labels) > 0
  dtb <- dtboundary(ants2oro(lesion_mask))

  # Permutation testing
  labels <- labelClusters(lesion_mask, minClusterSize = 27)
  probles <- antsImageClone(labels)
  avprob <- NULL
  maxles <- max(labels)

  if (maxles == 0) {
    warning("No identified lesions were large enough to consider for CVS detection.")
    return(NULL)
  }

  # Calculate "lesionness"
  frangi_image <- frangi(image = ants2oro(epi), mask = mask)
  frangi_image[frangi_image < 0] <- 0

  for (j in 1:maxles) {
    frangsub <- frangi_image[labels == j]
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

  return(list(candidate.lesions = labels, cvs.prob_map = probles, cvs.biomarker = mean(avprob)))
}
