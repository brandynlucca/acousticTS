################################################################################
# Sphere modal series solution
################################################################################
SPHMS <- function(object) {
  # Extract model parameters/inputs ============================================
  model_params <- extract(object, "model_parameters")$SPHMS
  parameters <- model_params$parameters
  acoustics <- parameters$acoustics
  medium <- model_params$medium
  body <- model_params$body
  # Pull out maximum 'm' specifically ==========================================
  m_limit <- acoustics$m_limit
  # Compute Delta (shell thickness) =============================================
  b_delta <- body$radius_shell - body$radius_fluid
  # Compute ka and kb ==========================================================
  k1a <- acoustics$k_sw * body$radius_shell
  k2a <- acoustics$k_shell * body$radius_shell
  k2b <- acoustics$k_shell * body$radius_fluid
  k3b <- acoustics$k_fluid * body$radius_fluid
  # Define the boundary expansion coefficient method ===========================
  .bm_wrap <- switch(
    parameters$Bm_method,
    Bm_rigid = function(f) .sphms_bm_rigid(k1a[f], m_limit[f]),
    Bm_pressure_release = function(f) .sphms_bm_prelease(k1a[f], m_limit[f]) ,
    Bm_fluid = function(f) {
      .sphms_bm_fluid(k1a[f], k2a[f], body$g31, body$h31, m_limit[f])
    },
    Bm_shelled_pressure_release = function(f) {
      .sphms_bm_shelled_prelease(
        k1a[f], k2a[f], k2b[f], body$g31, body$h31, m_limit[f]
      )
    },
    Bm_shelled_fluid = function(f) {
      .sphms_bm_shelled_fluid(
        k1a[f], k2a[f], k2b[f], k3b[f], body$g21, body$g32, body$h21, body$h32,
        m_limit[f]
      )
    }
  )
  sph_bm <- vapply(seq_len(nrow(acoustics)),
                   FUN = function(f) .bm_wrap(f),
                   FUN.VALUE = complex(1))
  # Compute the linear scattering coefficient, f_bs ============================
  f_bs <- -1i / acoustics$k_sw * sph_bm
  # Calculate backscatter and return ===========================================
  # Compute sigma_bs +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  o_bs <- abs(f_bs)
  methods::slot(object, "model")$SPHMS <- data.frame(
    frequency = acoustics$frequency,
    f_bs = f_bs,
    sigma_bs = o_bs,
    TS = 20 * log10(o_bs)
  )
  object
}

.sphms_bm_rigid <- function(k1a, m_limit) {
  # Expand k1a into matrix =====================================================
  k1a_m <- modal_matrix(k1a, m_limit)
  # Expand modal series iterator ===============================================
  m <- 0:m_limit
  # Calculate the expansion coefficient, A_m ===================================
  Am <- -(jsd(m_limit, k1a_m) / hsd(m_limit, k1a_m))
  # Return the entire boundary modal term ======================================
  sum((2 * m + 1) * (-1)^m * Am)
}

.sphms_bm_prelease <- function(k1a, m_limit) {
  # Expand k1a into matrix =====================================================
  k1a_m <- modal_matrix(k1a, m_limit)
  # Expand modal series iterator ===============================================
  m <- 0:m_limit
  # Calculate the expansion coefficient, A_m ===================================
  Am <- -(js(m_limit, k1a_m) / hs(m_limit, k1a_m))
  # Return the entire boundary modal term ======================================
  sum((2 * m + 1) * (-1)^m * Am)
}

.sphms_bm_fluid <- function(k1a, k2a, g31, h31, m_limit) {
  # Expand ka into matrices ====================================================
  k1a_m <- modal_matrix(k1a, m_limit)
  k3a_m <- modal_matrix(k3a, m_limit)
  # Expand modal series iterator ===============================================
  m <- 0:m_limit
  # Get material properties product ============================================
  gh <- g31 * h31
  # Calculate expansion ceofficient, C_m =======================================
  # Numerator ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  cm_num <- (jsd(m, k3a_m) * ys(m, k1a_m)) / (js(m, k3a_m) * jsd(m, k1a_m)) -
    (gh * ysd(m, k1a_m) / jsd(m, k1a_m))
  # Denominator ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  cm_denom <- (jsd(m, k3a_m) * js(m, k1a_m)) / (js(m, k3a_m) * jsd(m, k1a_m)) -
    gh
  # Quotient +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  cm <- cm_num / cm_denom
  # Calculate the expansion coefficient, A_m ===================================
  Am <- (-1 / (1 + 1i * cm))
  # Return the entire boundary modal term ======================================
  sum((2 * m + 1) * (-1)^m * Am)
}

.sphms_bm_shelled_prelease <- function(k1a, k2a, k2b, g31, h31, m_limit) {
  # Expand ka into matrices ====================================================
  k1a_m <- modal_matrix(k1a, m_limit)
  k2a_m <- modal_matrix(k2a, m_limit)
  k2b_m <- modal_matrix(k2b, m_limit)
  # Expand modal series iterator ===============================================
  m <- 0:m_limit
  # Get material properties product ============================================
  gh <- g31 * h31
  # Calculate the expansion coefficient, A_m ===================================
  # Compute elements for boundary matrix +++++++++++++++++++++++++++++++++++++++
  a11 <- -hs(m, k1a_m)
  a12 <- js(m, k1a_m)
  a21 <- -g21 * h21 * hsd(m, k1a_m)
  b1 <- js(m, k1a_m)
  b2 <- jsd(m, k1a_m) * g21 * h21
  d1 <- js(m, k2a_m) * ys(m, k2b_m) - js(m, k2b_m) * ys(m, k2a_m)
  d2 <- jsd(m, k2a_m) * ys(m, k2b_m) - js(m, k2b_m) * ysd(m, k2a_m)
  # Calculate coefficient ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Am <- (b1 * d2 - d1 * b2) / (a11 * d2 - d1 * a21)
  # Return the entire boundary modal term ======================================
  sum((2 * m + 1) * (-1)^m * Am)
}

.sphms_bm_shelled_fluid <- function(
    k1a, k2a, k2b, k3b, g21, g32, h21, h32, m_limit
) {
  # Expand modal series iterator ===============================================
  m <- 0:m_limit
  # Expand ka into matrices ====================================================
  k1a_m <- modal_matrix(k1a, m_limit)
  k2a_m <- modal_matrix(k2a, m_limit)
  k2b_m <- modal_matrix(k2b, m_limit)
  k3b_m <- modal_matrix(k3b, m_limit)
  # Get material properties product ============================================
  gh <- g31 * h31
  # Calculate the expansion coefficient, A_m ===================================
  # Compute elements for boundary matrix +++++++++++++++++++++++++++++++++++++++
  a11 <- -hs(m, k1a_m)
  a12 <- js(m, k2a_m)
  a13 <- ys(m, k2a_m)
  a21 <- -g21 * h21 * hsd(m, k1a_m)
  a22 <- jsd(m, k2a_m)
  a23 <- ysd(m, k2a_m)
  a32 <- js(m, k2b_m) * jsd(m, k3b_m) - g32 * h32 * jsd(m, k2b_m) * js(m, k3b_m)
  a33 <- ys(m, k2b_m) * jsd(m, k3b_m) - g32 * h32 * ysd(m, k2b_m) * js(m, k3b_m)
  b1 <- js(m, k1a_m)
  b2 <- jsd(m, k1a_m) * g21 * h21
  # Calculate coefficient ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Am <- (b1*a22*a33 + a13*b2*a32 - a12*b2*a33 - b1*a23*a32) /
    (a11*a22*a33 + a13*a21*a32 - a12*a21*a33 - a11*a23*a32)
  # Return the entire boundary modal term ======================================
  sum((2 * m + 1) * (-1)^m * Am)
}
