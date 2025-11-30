################################################################################
# Goodman and Stern (1962) modal series solution for elastic-shelled spheres
################################################################################
#' Calculates the theoretical TS of an elastic-shelled sphere using the modal
#' series solution from Goodman and Stern (1962).
#' @param object GAS- or SBF-class object.
#' @details
#' Calculates the theoretical TS of an elastic-shelled sphere using an exact
#' modal series solution
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Goodman, R.R., and Stern, R. (1962). Reflection and transmission of sound by
#' elastic spherical shells. The Journal of the Acoustical Society of America,
#' 34, 338-344.
#' @export
MSS_goodman_stern <- function(object) {
  # Extract model parameters/inputs ============================================
  model <- extract(object, "model_parameters")$MSS_goodman_stern
  # Get requisite elastic properties ===========================================
  G <- model$shell$G
  lambda <- model$shell$lambda
  density_shell <- model$shell$density
  # Extract the required morphometrics =========================================
  radius_shell <- model$shell$radius
  radius_fluid <- model$fluid$radius
  # Get the internal fluid density =============================================
  density_fluid <- model$fluid$density
  # Calculate shell sound speeds in the longitudinal and transverse directions =
  sound_speed_longitudinal <- sqrt((lambda + 2 * G) / density_shell)
  sound_speed_transversal <- sqrt(G / density_shell)
  # Calculate the associated wavenumbers =======================================
  kL <- k(model$parameters$acoustics$frequency, sound_speed_longitudinal)
  kT <- k(model$parameters$acoustics$frequency, sound_speed_transversal)
  # Calculate the reindexed ka values ==========================================
  ka_matrix <- calculate_ka_matrix(
    model$parameters$acoustics$frequency,
    model$medium$sound_speed,
    model$fluid$sound_speed,
    sound_speed_longitudinal,
    sound_speed_transversal,
    radius_shell,
    radius_fluid
  )
  # Get the modal series iterators =============================================
  m_limit <- model$parameters$m_limit
  m <- model$parameters$m
  # Expand the wavenumber matrix over the modal series limits ==================
  ka_matrix_m <- lapply(rownames(ka_matrix), function(ka) {
    modal_matrix(ka_matrix[ka, ], m_limit)
  })
  # Add the names ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  names(ka_matrix_m) <- rownames(ka_matrix)
  # Pre-calculate and cache the Bessel function outputs ========================
  bessel_cache <- .calculate_bessel_cache(ka_matrix_m, m)
  # Compute the alpha coefficient matrix =======================================
  alpha <- .goodman_stern_alpha(
    bessel_cache, ka_matrix_m, m,
    lambda, G, model$medium$density,
    density_shell, density_fluid
  )
  # Create the boundary conditions matrices ====================================
  boundary_matrices <- .goodman_stern_boundaries(
    alpha, ka_matrix, m
  )
  # Compute the modal series coefficient (b_m) using Cramer's Rule =============
  b_m <- lapply(seq_along(boundary_matrices), function(freq_idx) {
    vapply(seq_along(m), function(m_idx) {
      boundary_m <- boundary_matrices[[freq_idx]][[m_idx]]
      # Calculate determinants using QR decomposition for numerical stability ++
      det_numerator <- prod(
        abs(Re(diag(qr(boundary_m$A_numerator)$qr)))
      )
      det_denominator <- prod(
        abs(Re(diag(qr(boundary_m$A_denominator)$qr)))
      )
      # Apply Cramer's rule: b_m = det(A_numerator) / det(A_denominator) +++++++
      if (det_denominator == 0) {
        warning(
          paste(
            "Zero denominator at frequency",
            freq_idx, "mode", m[m_idx]
          )
        )
        return(NA)
      }

      det_numerator / det_denominator
    }, FUN.VALUE = numeric(1))
  })
  # Calculate the linear scattering coefficient, f_bs ==========================
  f_bs <- -1i / model$parameters$acoustics$k_sw *
    vapply(seq_along(b_m), function(ka_idx) {
      sum((-1)^m * (2 * m + 1) * b_m[[ka_idx]])
    }, FUN.VALUE = numeric(1))
  # Convert to sigma_bs ========================================================
  sigma_bs <- abs(f_bs)^2
  # Define MSS slot for ESS-type scatterer =====================================
  methods::slot(object, "model")$MSS_goodman_stern <- data.frame(
    frequency = model$parameters$acoustics$frequency,
    ka_shell = model$parameters$acoustics$k_sw * radius_shell,
    ka_fluid = model$parameters$acoustics$k_sw * radius_fluid,
    f_bs = f_bs,
    sigma_bs = sigma_bs,
    TS = db(sigma_bs)
  )
  # Return object ==============================================================
  object
}
