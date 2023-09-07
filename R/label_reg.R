#' @title Register Label Image Based on Full Image
#' @description This function registers a full image to a fixed image, then applies the registration to a label or binary image in the same space as the full image.
#' @param full_image an image of class \code{antsImage}, in the same space as the label image, which will be registered to the fixed image.
#' @param label_image an image of class \code{antsImage} with limited structural information, in the same space as the full image, to which the full image registration will be applied.
#' @param fixed_image an image of class \code{antsImage}, to which the other images will be registered.
#' @param typeofTransform the type of registration desired; this value is passed onto the registration function from extrantsr.
#'
#' @importFrom ANTsRCore antsApplyTransforms
#' @importFrom ANTsRCore antsRegistration
#' @return A list containing image_reg (the registered version of fullimage) and label_reg (the registered version of labelimage).
#' @examples \dontrun{
#' flair <- readnii("path/to/flair")
#' t1 <- readnii("path/to/t1")
#' tissue.class <- fast(t1, opts = "--nobias")
#' registered <- labelreg(t1, tissue.class, flair)
#' t1_reg <- registered$image_reg
#' tissue.class_reg <- registered$label_reg
#' }
#' @export
label_reg <- function(full_image, label_image, fixed_image,
                      typeofTransform = "Rigid") {
  registration_obj <- antsRegistration(
    fixed = fixed_image, moving = full_image,
    typeofTransform = typeofTransform
  )
  img_to_fix <- antsApplyTransforms(
    fixed = fixed_image, moving = full_image,
    transformlist = registration_obj$fwdtransforms,
    interpolator = "lanczosWindowedSinc"
  )
  lab_to_fix <- antsApplyTransforms(
    fixed = fixed_image, moving = label_image,
    transformlist = registration_obj$fwdtransforms,
    interpolator = "nearestNeighbor"
  )
  return(list(image_reg = img_to_fix, label_reg = lab_to_fix))
}
