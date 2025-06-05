#' Decompose a tumor profile into COSMIC ID signatures
#' 
#' @param tumor.ref Matrix of mutation counts (samples × 83 ID types)
#' @param sample.id Sample name (must be in rownames)
#' @param signatures.ref Matrix of COSMIC ID signatures (signatures × 83 ID types)
#' @param signature.cutoff Minimum weight to report a signature
#' @param associated Vector of associated signatures. If given, will narrow the 
#'   signatures tested to only the ones listed.
#' @return A list with weights, reconstructed profile, residuals, and unknown
#' @export

whichIDSignatures <- function(tumor.ref,
                              sample.id,
                              signatures.ref,
                              associated = c(),
                              signature.cutoff = 0.06) {
  # Basic checks
  if (!sample.id %in% rownames(tumor.ref)) {
    stop(paste(sample.id, "not found in tumor.ref rownames"))
  }
  
  tumor <- tumor.ref[sample.id, , drop = FALSE]
  tumor <- as.matrix(tumor)
  
  if (sum(tumor) == 0) stop("Tumor profile is empty.")
  
  # Normalize tumor profile
  tumor <- tumor / sum(tumor)
  
  # Check that signature matrix matches tumor profile
  common <- intersect(colnames(tumor), colnames(signatures.ref))
  if (length(common) == 0) {
    stop("No overlapping mutation types between tumor and signatures.")
  }
  
  tumor <- tumor[, common, drop = FALSE]
  signatures <- signatures.ref[, common, drop = FALSE]
  
  #Take a subset of the signatures
  if(!is.null(associated)){
    signatures <- signatures[rownames(signatures) %in% associated,,drop=FALSE ]
  }
  
  # NNLS to solve weights
  nnls_fit <- nnls::nnls(t(as.matrix(signatures)), as.numeric(tumor))
  weights <- as.vector(nnls_fit$x)
  weights <- weights / sum(weights)
  
  # Filter with cutoff
  weights[weights < signature.cutoff] <- 0
  weights <- weights / sum(weights)
  names(weights) <- rownames(signatures)
  
  # Reconstruct profile and calculate difference
  product <- t(weights) %*% as.matrix(signatures)
  diff <- tumor - product
  
  unknown <- 1 - sum(weights, na.rm = TRUE)
  
  x <- matrix(0, 
              nrow = 1, 
              ncol = nrow(signatures.ref), 
              dimnames = list(sample.id, rownames(signatures.ref)))
  x <- as.data.frame(x)
  x[names(weights)] <- weights
  weights <- x
  
  return(list(
    weights = weights,
    tumor = tumor,
    product = product,
    diff = diff,
    unknown = unknown
  ))
}