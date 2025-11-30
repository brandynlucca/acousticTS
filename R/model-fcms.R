################################################################################
# Fine cylinder modal series solution
################################################################################
FCMS <- function(object) {
  # Extract model parameters/inputs ============================================
  model_params <- extract(object, "model_parameters")$FCMS
  parameters <- model_params$parameters
  acoustics <- parameters$acoustics
  medium <- model_params$medium
  body <- model_params$body
  # Pull out maximum 'm' specifically ==========================================
  m_max <- acoustics$m_max
  # Pre-allocate Neumann factors ===============================================
  nu <- lapply(
    seq_len(nrow(acoustics)),
    function(f) neumann(0:acoustics$m_max[f])
  )
  # Compute kL and ka ==========================================================
  k1L <- body$length_body * acoustics$k_sw
  k1a <- acoustics$k_sw * sin(body$theta) * body$radius_body
  k2a <- acoustics$k_sw * sin(body$theta) / body$h * body$radius_body
  # Combine material properties ================================================
  gh <- body$g * body$h
  # Resolve modal series coefficient calculation method ========================
  if (parameters$Bm_method == "Bm_fluid") {
    .bm_wrap <- function(f) {
      # Get m vector ===========================================================
      m <- 0:m_max[f]
      # Convert to modal matrices ==============================================
      k1a_m <- modal_matrix(k1a[f], max(m))
      k2a_m <- modal_matrix(k2a[f], max(m))
      # Compute numerator for the boundary coefficient Cm ======================
      Cm_num <- (jcd(m, k2a_m)*yc(m, k1a_m)) /
        (jc(m, k2a_m)*jcd(m, k1a_m)) -
        gh*(ycd(m, k1a_m)/jcd(m, k1a_m))
      # Compute denominator for the boundary coefficient Cm ====================
      Cm_denom <- (jcd(m, k2a_m)*jc(m, k1a_m)) /
        (jc(m, k2a_m)*jcd(m, k1a_m)) - gh
      # Resolve Cm =============================================================
      Cm <- Cm_num / Cm_denom
      # Return Bm
      1i^(m+1) * (-nu[[f]] * 1i^(m) / (1 + 1i * Cm))
    }
  } else if (parameters$Bm_method == "Bm_rigid") {
    .bm_wrap <- function(f) {
      # Get m vector ===========================================================
      m <- 0:m_max[f]
      # Convert to modal matrix ================================================
      k1a_m <- modal_matrix(k1a[f], max(m))
      # Return Bm ==============================================================
      (-1)^m * nu[[f]] * (jcd(m, k1a_m) / hcd(m, k1a_m))
    }
  } else {
    .bm_wrap <- function(f) {
      # Get m vector ===========================================================
      m <- 0:m_max[f]
      # Convert to modal matrix ================================================
      k1a_m <- modal_matrix(k1a[f], max(m))
      # Return Bm ==============================================================
      (-1)^m * nu[[f]] * (jc(m, k1a_m) / hc(m, k1a_m))
    }
  }
  fc_bm <- lapply(seq_len(nrow(acoustics)), FUN = function(f) .bm_wrap(f))
  # Compute the linear scattering coefficient, f_bs ============================
  f_bs <- vapply(
    seq_len(nrow(acoustics)),
    FUN = function(f) {
      if (parameters$Bm_method == "Bm_fluid") {
        -body$length_body / pi *
          (sin(k1L[f] * cos(body$theta_body)) /
             (k1L[f] * cos(body$theta_body))) *
          sum(fc_bm[[f]])
      } else {
        1i * body$length_body / pi *
          (sin(k1L[f] * cos(body$theta_body)) /
             (k1L[f] * cos(body$theta_body))) *
          sum(fc_bm[[f]])
      }
    },
    FUN.VALUE = complex(1)
  )
  # Calculate backscatter and return ===========================================
  # Compute sigma_bs +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  o_bs <- abs(f_bs)
  methods::slot(object, "model")$FCMS <- data.frame(
    frequency = acoustics$frequency,
    f_bs = f_bs,
    sigma_bs = o_bs,
    TS = 20 * log10(o_bs)
  )
  object
}
