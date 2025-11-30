################################################################################
# Kirchoff-Ray Mode approximation
################################################################################
#' Calculates the theoretical TS using Kirchoff-ray Mode approximation.
#'
#' @param object Desired object/animal shape. Must be class "SBF".
#' @usage
#' KRM(object)
#' @details
#' Calculates the theoretical TS using the Kirchoff-ray Mode model.
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Clay C.S. and Horne J.K. (1994). Acoustic models of fish: The Atlantic cod
#' (Gadus morhua). Journal of the Acoustical Society of AMerica, 96, 1661-1668.
#' @export
KRM <- function(object) {
  # Detect object class ========================================================
  scatterer_type <- class(object)
  # Extract model parameter inputs =============================================
  model <- extract(object, "model_parameters")$KRM
  # Extract body parameters ====================================================
  body <- extract(object, "body")
  # Calculate reflection coefficient for medium-body interface =================
  R12 <- reflection_coefficient(model$medium, model$body)
  # Calculate transmission coefficient and its reverse =========================
  T12T21 <- 1 - R12 * R12
  # Sum across body position vector ============================================
  rpos <- switch(scatterer_type,
                 FLS = rbind(
                   x = body$rpos[1, ],
                   w = c(
                     body$radius[2],
                     body$radius[2:(length(body$radius) - 1)],
                     body$radius[(length(body$radius)) - 1]
                   ) * 2,
                   zU = c(
                     body$radius[2],
                     body$radius[2:(length(body$radius) - 1)],
                     body$radius[(length(body$radius)) - 1]
                   ),
                   zL = -c(
                     body$radius[2],
                     body$radius[2:(length(body$radius) - 1)],
                     body$radius[(length(body$radius)) - 1]
                   )
                 ),
                 SBF = body$rpos
  )
  body_rpos_sum <- along_sum(rpos, model$parameters$ns_b)
  # Approximate radius of body cylinders =======================================
  a_body <- switch(scatterer_type,
                   FLS = body_rpos_sum[2, ] / 4,
                   SBF = body_rpos_sum[2, ] / 4
  )
  # Combine wavenumber (k) and radii to calculate "ka" =========================
  ka_body <- matrix(
    data = rep(a_body,
               each = length(model$parameters$acoustics$k_sw)
    ),
    ncol = length(a_body),
    nrow = length(model$parameters$acoustics$k_b)
  ) * model$parameters$acoustics$k_b
  # Convert c-z coordinates to required u-v rotated coordinates ================
  uv_body <- acousticTS::body_rotation(
    body_rpos_sum,
    body$rpos,
    body$theta,
    length(model$parameters$acoustics$k_sw)
  )
  # Calculate body empirical phase shift function ==============================
  body_dorsal_sum <- switch(scatterer_type,
                            FLS = matrix(
                              data = rep(body_rpos_sum[3, ],
                                         each = length(model$parameters$acoustics$k_sw)
                              ),
                              ncol = length(body_rpos_sum[3, ]),
                              nrow = length(model$parameters$acoustics$k_sw)
                            ) / 2,
                            SBF = matrix(
                              data = rep(body_rpos_sum[3, ],
                                         each = length(model$parameters$acoustics$k_sw)
                              ),
                              ncol = length(body_rpos_sum[3, ]),
                              nrow = length(model$parameters$acoustics$k_sw)
                            ) / 2
  )
  Psi_b <- -pi * model$parameters$acoustics$k_b * body_dorsal_sum /
    (2 * (model$parameters$acoustics$k_b * body_dorsal_sum + 0.4))
  # Estimate natural log function (phase, etc.) ================================
  exp_body <- exp(-2i * model$parameters$acoustics$k_sw * uv_body$vbU) -
    T12T21 * exp(-2i * model$parameters$acoustics$k_sw * uv_body$vbU +
                   2i * model$parameters$acoustics$k_b *
                   (uv_body$vbU - uv_body$vbL) + 1i * Psi_b)
  # Resolve summation term =====================================================
  body_summation <- sqrt(ka_body) * uv_body$delta_u
  # Calculate linear scattering length (m) =====================================
  f_body <- rowSums(-((1i * (R12 / (2 * sqrt(pi)))) *
                        body_summation * exp_body))
  if (scatterer_type == "FLS") {
    # Define KRM slot for FLS-type scatterer ===================================
    methods::slot(object, "model")$KRM <- data.frame(
      frequency = model$parameters$acoustics$frequency,
      ka = model$parameters$acoustics$k_sw *
        stats::median(a_body, na.rm = TRUE),
      f_bs = f_body,
      sigma_bs = abs(f_body) * abs(f_body),
      TS = 20 * log10(abs(f_body))
    )
  } else if (scatterer_type == "SBF") {
    #### Repeat process for bladder ============================================
    # Extract bladder parameters ===============================================
    bladder <- acousticTS::extract(object, "bladder")
    # Calculate reflection coefficient for bladder =============================
    R23 <- acousticTS::reflection_coefficient(
      body,
      bladder
    )
    # Sum across body/swimbladder position vectors =============================
    bladder_rpos_sum <- acousticTS::along_sum(
      bladder$rpos,
      model$parameters$ns_sb
    )
    # Approximate radii of bladder discrete cylinders ==========================
    a_bladder <- bladder_rpos_sum[2, ] / 4
    # Combine wavenumber (k) and radii to calculate "ka" for bladder ===========
    ka_bladder <- matrix(
      data = rep(a_bladder, each = length(model$parameters$acoustics$k_sw)),
      ncol = length(a_bladder),
      nrow = length(model$parameters$acoustics$k_sw)
    ) * model$parameters$acoustics$k_sw
    # Calculate Kirchoff approximation empirical factor, A_sb ==================
    A_sb <- ka_bladder / (ka_bladder + 0.083)
    # Calculate empirical phase shift for a fluid cylinder, Psi_p ==============
    Psi_p <- ka_bladder / (40 + ka_bladder) - 1.05
    # Convert x-z coordinates to requisite u-v rotated coordinates =============
    uv_bladder <- acousticTS::bladder_rotation(
      bladder_rpos_sum,
      bladder$rpos,
      bladder$theta,
      length(model$parameters$acoustics$k_sw)
    )
    # Estimate natural log functions  ==========================================
    exp_bladder <- exp(-1i * (2 * model$parameters$acoustics$k_b *
                                uv_bladder$v + Psi_p)) * uv_bladder$delta_u
    # Calculate the summation term =============================================
    bladder_summation <- A_sb * sqrt((ka_bladder + 1) * sin(bladder$theta))
    # Estimate backscattering length, f_fluid/f_soft ===========================
    f_bladder <- rowSums(-1i * (R23 * T12T21) / (2 * sqrt(pi)) *
                           bladder_summation * exp_bladder)
    # Estimate total backscattering length, f_bs ===============================
    f_bs <- f_body + f_bladder
    # Define KRM slot for FLS-type scatterer ===================================
    methods::slot(object, "model")$KRM <- data.frame(
      frequency = model$parameters$acoustics$frequency,
      ka = model$parameters$acoustics$k_sw *
        stats::median(a_body, na.rm = TRUE),
      f_body = f_body,
      f_bladder = f_bladder,
      f_bs = f_bs,
      sigma_bs = abs(f_bs) * abs(f_bs),
      TS = 20 * log10(abs(f_bs))
    )
  }
  # Return object ==============================================================
  object
}
