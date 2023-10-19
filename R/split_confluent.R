#' Split Confluent Lesions
#'
#' This function splits confluent lesions labeled as "i" in the labeled image
#' using the lesion center information and returns a labeled image with split lesions.
#'
#' @param i The label of the confluent lesions to be split.
#' @param labeled_image The labeled image containing lesions.
#' @param lesion_center_image The image containing lesion center information.
#'
#' @return A labeled image with split lesions.
#' @export
#'
#' @examples \dontrun{
#' labeled_image <- check_ants("labeled_image.nii.gz")
#' lesion_center_image <- check_ants("lesion_center_image.nii.gz")
#' split_lesion_image <- split_confluent(i = 1, labeled_image, lesion_center_image)
#' }
#'
#' @seealso
#' \code{\link{get_lesion_labels}} for obtaining lesion labels based on lesion centers.
#'
#' @description
#' This function takes a labeled image, lesion center image, and the label "i" of the confluent lesion
#' to be split. It uses the lesion center information to split the confluent lesions and returns a labeled
#' image with the split lesions labeled distinctly.
#'
#' @export
split_confluent <- function(i, labeled_image, lesion_center_image, mincluster) {
  lesion <- labeled_image == i
  centers_in_label <- lesion_center_image[lesion]
  n_centers <- unique(centers_in_label[centers_in_label != 0])

  if (length(centers_in_label) < mincluster) {
    return(NULL)
  }
  
  if (length(n_centers) <= 1) {
    return(lesion)
  }

  split_lesion <- get_lesion_labels(
    lesion,
    lesion_center_image * lesion
  )
  s <- table(split_lesion[split_lesion != 0])
  s_large <- s[s >= mincluster]
  
  for (j in 1:length(s_large)) { # relabel to start from 1 and count up continuously
    split_lesion[split_lesion == as.numeric(names(s_large)[j])] <- j
  }
  return(split_lesion)
}
