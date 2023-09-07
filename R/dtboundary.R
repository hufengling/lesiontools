#' @title Distance to Mask Boundary
#' @description This function finds the distance of each voxel to the nearest boundary in a given mask.
#' @param mask a 3D array or image of class \code{nifti}, containing a binary mask where 1 represents structure.
#'
#' @return A new image in which voxels have been assigned their distance to the nearest boundary.
#' @examples \dontrun{
#' library(neurobase)
#' lesion.mask <- readnii("path/to/mask")
#' dtb <- dtboundary(mask = lesion.mask)
#' }
#' @export
dtboundary <- function(mask) {
  get.d.boundary.exact.balloon <- function(v, mask, d.max = 30) {
    if (mask[v[1], v[2], v[3]] == 0)
      stop("ERROR! - voxel outside of mask...")
    r <- 1
    # expand balloon
    while (TRUE) {
      balloon <- 1 - mask[(v[1] - r):(v[1] + r), (v[2] - r):(v[2] + r), (v[3] - r):(v[3] + r)]
      if (sum(balloon > 0)) {
        which.outside <- which(balloon > 0, arr.ind = TRUE) - (r + 1)
        d.out <- sqrt(min(rowSums(which.outside^2)))
        break
      }
      if (r > d.max) {
        d.out <- 1000
        break
      }
      r <- r + 1
    }
    return(d.out)
  }
  which.mask.arrind <- which(mask > 0, arr.ind = TRUE)
  # For each voxel in the mask
  min.d <- rep(0, nrow(which.mask.arrind))
  for (i in 1:length(min.d)) {
    # Get minimum distance to boundary
    min.d[i] <- get.d.boundary.exact.balloon(which.mask.arrind[i, ], mask)
  }
  mask[mask > 0] <- min.d
  return(mask)
}
