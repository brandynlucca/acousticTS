################################################################################
# Prolate spheroid modal series solution
################################################################################
#' Calculates the theoretical target strength of a prolate spheroid using a
#' modal series solution
#'
#' @param object Prolate spheroid object.
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Furusawa, M. (1988). Prolate spheroidal models for predicting general trends
#' of fish target strength. Journal of the Acoustical Society of Japan, 9:
#' 13-24.
#'
#' @export
PSMS <- function(object) {
  # Extract model parameters/inputs ============================================
  model_params <- extract(object, "model_parameters")$PSMS
  parameters <- model_params$parameters
  acoustics <- parameters$acoustics
  medium <- model_params$medium
  body <- model_params$body
  # Pre-allocate Neumann factors ===============================================
  nu <- lapply(
    seq_len(nrow(acoustics)),
    function(f) neumann(0:acoustics$m_max[f])
  )
  # Pre-allocate azimuthal phase relations =====================================
  azimuth <- lapply(
    seq_len(nrow(acoustics)),
    function(f) {
      cos(0:acoustics$m_max[f] * (body$phi_body - body$phi_scatter))
    }
  )
  # Compute the expansion matrix Amn ===========================================
  if (parameters$Amn_method == "Amn_liquid_simplify") {
    .amn_wrap <- function(f) {
      # Compute using the simplified expression from Eq. 5 (Furusawa, 1988) ++++
      liquid_spheroidal_simplified_expansion(
        acoustics$m_max[f], acoustics$n_max[f], acoustics$chi_sw[f],
        acoustics$chi_body[f], body$xi, medium$density, body$density
      )$amn
    }
  } else if (parameters$Amn_method == "Amn_liquid") {
    .amn_wrap <- function(f) {
      # Compute the kernel matrices for a liquid-filled spheroidal scatterer +++
      kernels <- liquid_spheroidal_kernels(
        acoustics$m_max[f], acoustics$n_max[f], acoustics$chi_sw[f],
        acoustics$chi_body[f], body$theta_body, body$xi,
        body$density, medium$density
      )
      # Solve for the expansion coefficient Amn ++++++++++++++++++++++++++++++++
      solve_liquid_spheroidal_Amn(
        kernels$K1_kernel, kernels$K3_kernel
      )
    }
  } else if (parameters$Amn_method == "Amn_pressure_release") {
    .amn_wrap <- function(f) {
      # Compute the external radial incident wave function +++++++++++++++++++++
      R_incident <- radial_external_incoming_matrix(
        acoustics$m_max[f],
        acoustics$n_max[f], acoustics$chi_sw[f], body$xi
      )
      # Compute the external radial scattering wave function +++++++++++++++++++
      R_scattering <- radial_external_scattering_matrix(
        acoustics$m_max[f],
        acoustics$n_max[f], acoustics$chi_sw[f], body$xi
      )
      # Compute the expansion matrix +++++++++++++++++++++++++++++++++++++++++++
      - R_incident$value / R_scattering$value
    }
  } else {
    .amn_wrap <- function(f) {
      # Compute the first deriv of the external radial incident wave function ++
      dR_incident <- radial_external_incoming_matrix(
        acoustics$m_max[f],
        acoustics$n_max[f], acoustics$chi_sw[f], body$xi
      )
      # Compute the first deriv of the external radial scattering wave function
      dR_scattering <- radial_external_scattering_matrix(
        acoustics$m_max[f],
        acoustics$n_max[f], acoustics$chi_sw[f], body$xi
      )
      # Compute the expansion matrix +++++++++++++++++++++++++++++++++++++++++++
      - dR_incident$derivative / dR_scattering$derivative
    }
  }
  ps_amn <- lapply(
    seq_len(nrow(acoustics)),
    FUN = function(f) .amn_wrap(f)
  )
  # Compute the linear scattering coefficient, f_bs ============================
  f_bs <- vapply(
    seq_len(nrow(acoustics)),
    FUN = function(f){
      .psms_fbs(
        nu[[f]], azimuth[[f]], acoustics$m_max[f], acoustics$n_max[f],
        acoustics$chi_sw[f],
        body$theta_body, body$theta_scatter, ps_amn[[f]],
        ifelse(
          parameters$Amn_method %in% c("Amn_liquid", "Amn_liquid_simplified"),
          TRUE,
          FALSE
        )
      )
    },
    FUN.VALUE = complex(1)
  )
  # Calculate backscatter and return ===========================================
  # Compute sigma_bs +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  o_bs <- abs(-2i / acoustics$k_sw * f_bs)
  methods::slot(object, "model")$PSMS <- data.frame(
    frequency = acoustics$frequency,
    f_bs = f_bs,
    sigma_bs = o_bs,
    TS = 20 * log10(o_bs)
  )
  object
}

.psms_fbs <- function(nu, azimuth, m_max, n_max, chi_sw, theta_body,
                      theta_scatter, ps_amn, list2mat = FALSE) {
  # Pre-compute the product of the regular and scatter Smn matrices ============
  smn_matrix <- outer(
    0:m_max, 0:n_max,
    Vectorize(function(m, n) {
      if (n < m) return(NA)
      # Regular
      Smn(m, n, chi_sw, cos(theta_body), normalize = TRUE)$value *
        # Scattering
        Smn(m, n, chi_sw, cos(theta_scatter), normalize = TRUE)$value
    })
  )

  # Convert list to matrix if defined ==========================================
  if (list2mat) {
    # Get the matrix dimensions and pre-allocated 'Amn_matrix' +++++++++++++++++
    smn_matrix_dims <- dim(smn_matrix)
    nrows <- smn_matrix_dims[1]
    ncols <- smn_matrix_dims[2]
    Amn_matrix <- matrix(0 + 0i,
                         nrow = nrows,
                         ncol = ncols)
    # Reshape Amn list into matrix +++++++++++++++++++++++++++++++++++++++++++++
    seq_dim <- seq_len(min(length(ps_amn), nrows))
    Amn_rows <- vapply(
      seq_dim, FUN = function(r) {
        vec <- complex(length = ncols)
        am_vec <- as.vector(ps_amn[[r]])
        if (length(am_vec) > 0L) {
          end_col <- min(smn_matrix_dims[2], r + length(am_vec) - 1L)
          cols <- r:end_col
          vec[cols] <- am_vec[seq_len(length(cols))]
        }
        vec
      }, FUN.VALUE = complex(smn_matrix_dims[2]), USE.NAMES = FALSE
    )
    Amn_matrix[seq_dim, ] <- t(Amn_rows)
  } else {
    Amn_matrix <- ps_amn
  }
  # Calculate the linear scattering coefficient ================================
  sum(nu * smn_matrix * Amn_matrix * azimuth, na.rm = TRUE)
}
