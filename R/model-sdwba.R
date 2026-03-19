################################################################################
# Stochastic distorted wave Born approximation (SDWBA)
################################################################################
#' Stochastic distorted wave Born approximation (SDWBA) for weak scatterers
#'
#' @description
#' Calculates the far-field scattering amplitude and related quantities for
#' fluid-like, weak scatterers using the stochastic distorted wave Born
#' approximation (SDWBA), as described by Demer and Conti (2003). The SDWBA
#' extends the deterministic DWBA by incorporating stochastic phase variability
#' to account for unresolved structural complexity and dynamic variability in
#' biological scatterers.
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model = "sdwba",
#'   n_iterations,
#'   n_segments_init,
#'   phase_sd_init,
#'   length_init,
#'   frequency_init,
#'   sound_speed_sw,
#'   density_sw
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{n_iterations}}{Number of stochastic realizations for averaging
#'   target strength predictions.}
#'   \item{\code{n_segments_init}}{Reference number of body segments.}
#'   \item{\code{phase_sd_init}}{Reference phase deviation (radians).}
#'   \item{\code{length_init}}{Reference body length (m).}
#'   \item{\code{frequency_init}}{Reference frequency (Hz).}
#'   \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#' }
#'
#' @section Theory:
#' The SDWBA is derived under the weak scattering assumption, where the
#' differences in compressibility (\eqn{\kappa}) and density (\eqn{\rho})
#' between the scatterer and the surrounding fluid are small enough to
#' linearize the acoustic scattering problem (see \code{\link{DWBA}}). in this
#' regime, multiple scattering within the body is neglected, and the total
#' scattered field is approximated as the coherent sum of first-order
#' contributions from individual body segments.
#'
#' The key extension introduced by the SDWBA is the inclusion of stochastic
#' phase variability to represent unresolved morphological complexity, internal
#' inhomogeneity, and dynamic effects such as body flexure and orientation
#' variability. The linear scattering coefficient is written as:
#'
#'  \deqn{
#'    f_{bs}(\theta) = \sum\limits_{j=1}^N f_{bs}^{(j)}(\theta)
#'    \exp(i \varphi_j),
#'  }
#'
#' where \eqn{N} is the number of body segments, \eqn{f_{bs}^{(j)}} is the
#' contribution from segment \eqn{j}, and \eqn{\varphi_j} is a random phase
#' perturbation drawn independently for each segment.
#'
#' The phase perturbations are assumed to follow a zero-mean Gaussian
#' distribution with variance related to the effective signal-to-noise ratio
#' (SNR) of the scattering process. The minimum expected phase variance due to
#' noise is given by:
#'
#'  \deqn{
#'    \mathbb{V}(\varphi_j) = \frac{1}{2 \mathrm{SNR}},
#'  }
#'
#' though in practice larger variances are used to account for additional
#' physical sources of phase decorrelation not explicitly modeled.
#'
#' The expected backscattering cross-section is obtained by ensemble averaging
#' over multiple stochastic realizations:
#'  \deqn{
#'    \langle \sigma_{bs}(\theta) \rangle
#'    = \mathbb{E}\!\left[ \left| f_{bs}(\theta) \right|^2 \right]
#'    \approx \frac{1}{M} \sum_{m=1}^{M}
#'    \left| f_{bs}^{(m)}(\theta) \right|^2,
#'  }
#'
#' and the expected target strength, \eqn{\mathbb{E}[TS(\theta)]}, is computed
#' from this mean.
#'
#' To ensure consistency across frequencies and body sizes, the SDWBA enforces
#' scale invariance by preserving the product of the phase standard deviation,
#' \eqn{\mathrm{sd}_\varphi}, and frequency, \eqn{f}:
#'
#'  \deqn{
#'    \mathrm{sd}_{\varphi}(f)\, f =
#'    \mathrm{sd}_{\varphi_0}\, f_0,
#'  }
#'
#' and by scaling the number of segments to maintain constant spatial
#' resolution relative to acoustic wavelength:
#'
#'  \deqn{
#'    N(f, L) = N_0 \frac{f L}{f_0 L_0}.
#'  }
#'
#' The phase standard deviation at arbitrary frequency and length is then:
#'
#'  \deqn{
#'    \mathrm{sd}_{\varphi}(f, L) =
#'    \mathrm{sd}_{\varphi_0}
#'    \frac{N_0 L}{N(f, L) L_0}.
#'  }
#'
#' These scaling relationships ensure that stochastic decorrelation effects
#' remain physically consistent across different acoustic and geometric regimes.
#'
#' @section Implementation:
#' The implementation extracts geometric and acoustic parameters from the input
#' object, constructs the required rotation and wavenumber matrices, and
#' evaluates the DWBA contribution for each segment. For each stochastic
#' realization, random phase perturbations are applied, and the resulting
#' backscattering amplitudes are averaged over all realizations to estimate the
#' expected target strength.
#'
#' @seealso
#' See the \link[=boundary_conditions]{boundary conditions documentation} for
#' more details on weak scattering assumptions,
#' \code{\link{target_strength}}, \code{\link{FLS}}, \code{\link{DWBA}}
#'
#' @references
#' Conti, D.A., and Conti, S.G. (2006). Improved parameterization of the SDWBA
#' for estimating krill target strength. ICES Journal of Marine Science, 63:
#' 928-935.
#'
#' Demer, D.A., and Conti, S.G. (2003). Reconciling theoretical versus
#' empirical target strengths of krill: effects of phase variability on the
#' distorted-wave Born approximation. ICES Journal of Marine Science, 60:
#' 429-434.
#'
#' Stanton, T.K., Chu, D., and Wiebe, P.H. (1998). Sound scattering by several
#' zooplankton groups. II. Scattering models. The Journal of the Acoustical
#' Society of America, 103, 236-253.
#'
#' @name SDWBA
#' @aliases sdwba SDWBA sdwba_curved SDWBA_CURVED
#' @docType data
#' @keywords models acoustics
NULL

#' Initialize FLS-class object for SDWBA modeling
#' @param object FLS-class object.
#' @param frequency Frequency vector (Hz).
#' @param sound_speed_sw Seawater sound speed (m/s).
#' @param density_sw Seawater density (kg/m³).
#' @param n_iterations Number of iterations to repeat SDWBA
#' @param n_segments_init Reference number of body segments
#' @param phase_sd_init Reference phase deviation
#' @param length_init Reference body length
#' @param frequency_init Reference frequency
#' @noRd
sdwba_initialize <- function(object,
                             frequency,
                             sound_speed_sw = 1500,
                             density_sw = 1026,
                             n_iterations = 100,
                             n_segments_init = 14,
                             phase_sd_init = sqrt(2) / 2,
                             length_init = 38.35e-3,
                             frequency_init = 120e3) {
  # Parse shape ================================================================
  shape <- acousticTS::extract(object, "shape_parameters")
  # Parse body =================================================================
  body <- acousticTS::extract(object, "body")
  contrasts <- .derive_contrasts(body, sound_speed_sw, density_sw)
  body$h <- body_h <- contrasts$h
  body$g <- body_g <- contrasts$g
  # Define medium parameters ===================================================
  medium_params <- data.frame(
    sound_speed = sound_speed_sw,
    density = density_sw
  )
  # Define model parameters recipe =============================================
  model_params <- list(
    acoustics = data.frame(
      frequency = frequency,
      # Wavenumber (medium) ====================================================
      k_sw = acousticTS::wavenumber(
        frequency,
        sound_speed_sw
      ),
      # Wavenumber (fluid) =====================================================
      k_f = acousticTS::wavenumber(
        frequency,
        body_h * sound_speed_sw
      )
    ),
    n_segments = shape$n_segments
  )
  # Define stochastic recipe ===================================================
  # First calculate new resampled body shape resolution ++++++++++++++++++++++++
  N_f <- ceiling(n_segments_init * (frequency / frequency_init) *
                   (shape$length / length_init))
  N_f_vec <- ifelse(N_f > n_segments_init,
                    N_f,
                    n_segments_init
  )
  N_f_idx <- unique(
    N_f_vec
  )
  # Now we calculate the new phase standard deviation value ++++++++++++++++++++
  phase_sd <- phase_sd_init * (n_segments_init / N_f_vec) *
    (shape$length / length_init)
  # Create stochastic recipe +++++++++++++++++++++++++++++++++++++++++++++++++++
  stochastic_params <- lapply(seq_along(N_f_idx),
                              FUN = function(i) {
                                idx <- which(N_f_vec == N_f_idx[i])
                                object_new <- sdwba_resample(object,
                                                             n_segments = N_f_idx[i]
                                )
                                body <- acousticTS::extract(object_new, "body")
                                n_segments <- N_f_idx[i]
                                phase_sd <- phase_sd[i]
                                acoustics <- model_params$acoustics[idx, ]
                                list(
                                  meta_params = data.frame(
                                    n_iterations = n_iterations,
                                    phase_sd = phase_sd,
                                    N0 = n_segments_init,
                                    f0 = frequency_init,
                                    L0 = length_init,
                                    p0 = phase_sd_init
                                  ),
                                  body_params = body,
                                  n_segments = n_segments,
                                  acoustics = acoustics
                                )
                              }
  )
  # Add model parameters slot to scattering object =============================
  methods::slot(
    object,
    "model_parameters"
  )$SDWBA <- list(
    parameters = stochastic_params,
    medium = medium_params,
    body = body
  )
  # Add model results slot to scattering object ================================
  methods::slot(
    object,
    "model"
  )$SDWBA <- data.frame(
    frequency = frequency,
    f_bs = rep(
      NA,
      length(frequency)
    ),
    sigma_bs = rep(
      NA,
      length(frequency)
    ),
    TS = rep(
      NA,
      length(frequency)
    )
  )
  # Output =====================================================================
  return(object)
}


#' Summarize stochastic SDWBA realizations from deterministic segment integrals
#' @param segment_integrals Complex matrix with one row per frequency and one
#' column per body segment.
#' @param phase_sd Phase standard deviation in radians.
#' @param n_iterations Number of stochastic realizations.
#' @keywords internal
#' @noRd
sdwba_stochastic_summary <- function(segment_integrals, phase_sd, n_iterations) {
  n_k <- nrow(segment_integrals)
  n_segments <- ncol(segment_integrals)
  phase_cyl <- vapply(
    seq_len(n_k),
    FUN = function(x) {
      cyl_phase <- array(
        vapply(
          seq_len(n_segments),
          FUN = function(y) {
            phase_rng <- stats::rnorm(n_iterations, mean = 0, sd = 1)
            segment_integrals[x, y] * exp(1i * phase_rng * phase_sd)
          },
          FUN.VALUE = complex(n_iterations)
        ),
        dim = c(n_iterations, n_segments)
      )
      rowSums(cyl_phase, na.rm = TRUE)
    },
    FUN.VALUE = complex(n_iterations)
  )
  phase_cyl <- array(phase_cyl, dim = c(n_iterations, n_k))
  data.frame(
    f_bs = colMeans(phase_cyl),
    sigma_bs = colMeans(.sigma_bs(phase_cyl)),
    TS_mean = db(colMeans(.sigma_bs(phase_cyl))),
    TS_sd = db(apply(.sigma_bs(phase_cyl), 2, stats::sd))
  )
}

#' Evaluate one resampled SDWBA group
#' @param sub_params One element from the SDWBA resampling recipe.
#' @param theta Incident angle in radians.
#' @param R Weak-scattering contrast factor.
#' @param h Sound-speed contrast.
#' @keywords internal
#' @noRd
sdwba_resampled_backscatter <- function(sub_params, theta, R, h) {
  rpos <- rbind(
    sub_params$body_params$rpos[1:3, ],
    sub_params$body_params$radius
  )
  segment_integrals <- dwba_segment_integrals_cpp(
    rpos = rpos,
    k_sw = sub_params$acoustics$k_sw,
    theta = theta,
    h = h,
    R = R
  )
  sdwba_stochastic_summary(
    segment_integrals = segment_integrals,
    phase_sd = sub_params$meta_params$phase_sd,
    n_iterations = sub_params$meta_params$n_iterations
  )
}


#' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' using the stochastic distorted Born wave approximation (DWBA) model.
#' @param object FLS-class scatterer.
#' @noRd
SDWBA <- function(object) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract(object, "model_parameters")$SDWBA
  body <- acousticTS::extract(object, "body")
  theta <- body$theta
  # Material properties calculation ============================================
  g <- body$g
  h <- body$h
  R <- 1 / (g * h * h) + 1 / g - 2
  SDWBA_resampled <- function(i) {
    sub_params <- model$parameters[[i]]
    sdwba_resampled_backscatter(
      sub_params = sub_params,
      theta = theta,
      R = R,
      h = h
    )
  }
  # Generate results dataframe by collating resampled results ==================
  results <- do.call(
    "rbind",
    lapply(seq_along(model$parameters),
           FUN = function(i) SDWBA_resampled(i)
    )
  )
  # Update scatterer object ====================================================
  methods::slot(object, "model")$SDWBA$f_bs <- results$f_bs
  methods::slot(object, "model")$SDWBA$sigma_bs <- results$sigma_bs
  methods::slot(object, "model")$SDWBA$TS <- results$TS_mean
  methods::slot(object, "model")$SDWBA$TS_sd <- results$TS_sd
  # Return object ==============================================================
  object
}

#' Initialize FLS-class object for SDWBA modeling
#' @param object FLS-class object.
#' @param object FLS-class object.
#' @param frequency Frequency vector (Hz).
#' @param sound_speed_sw Seawater sound speed (m/s).
#' @param density_sw Seawater density (kg/m³).
#' @param n_iterations Number of iterations to repeat SDWBA
#' @param n_segments_init Reference number of body segments
#' @param phase_sd_init Reference phase deviation
#' @param length_init Reference body length
#' @param frequency_init Reference frequency
#' @noRd
sdwba_curved_initialize <- function(object,
                                    frequency,
                                    sound_speed_sw = 1500,
                                    density_sw = 1026,
                                    n_iterations = 100,
                                    n_segments_init = 14,
                                    phase_sd_init = sqrt(2) / 2,
                                    length_init = 38.35e-3,
                                    frequency_init = 120e3) {
  # Parse shape ================================================================
  shape <- acousticTS::extract(object, "shape_parameters")
  # Parse body =================================================================
  body <- acousticTS::extract(object, "body")
  # Bend body ==================================================================
  object_copy <- object
  object_copy <- brake(object_copy,
                       radius_curvature = body$radius_curvature_ratio,
                       mode = "ratio"
  )
  # Define medium parameters ===================================================
  medium_params <- data.frame(
    sound_speed = sound_speed_sw,
    density = density_sw
  )
  # Define model parameters recipe =============================================
  model_params <- list(
    acoustics = data.frame(
      frequency = frequency,
      # Wavenumber (medium) ====================================================
      k_sw = acousticTS::wavenumber(
        frequency,
        sound_speed_sw
      ),
      # Wavenumber (fluid) =====================================================
      k_f = acousticTS::wavenumber(
        frequency,
        body$h * sound_speed_sw
      )
    ),
    n_segments = shape$n_segments
  )
  # Define stochastic recipe ===================================================
  # First calculate new resampled body shape resolution ++++++++++++++++++++++++
  N_f <- ceiling(n_segments_init * (frequency / frequency_init) *
                   (shape$length / length_init))
  N_f_vec <- ifelse(N_f > n_segments_init,
                    N_f,
                    n_segments_init
  )
  N_f_idx <- unique(
    N_f_vec
  )
  # Now we calculate the new phase standard deviation value ++++++++++++++++++++
  phase_sd <- phase_sd_init * (n_segments_init / N_f_vec) *
    (shape$length / length_init)
  # Create stochastic recipe +++++++++++++++++++++++++++++++++++++++++++++++++++
  stochastic_params <- lapply(seq_along(N_f_idx),
                              FUN = function(i) {
                                idx <- which(N_f_vec == N_f_idx[i])
                                object_new <- sdwba_resample(object_copy,
                                                             n_segments = N_f_idx[i]
                                )
                                body <- extract(object_new, "body")
                                n_segments <- N_f_idx[i]
                                phase_sd <- phase_sd[i]
                                acoustics <- model_params$acoustics[idx, ]
                                list(
                                  meta_params = data.frame(
                                    n_iterations = n_iterations,
                                    phase_sd = phase_sd,
                                    N0 = n_segments_init,
                                    f0 = frequency_init,
                                    L0 = length_init,
                                    p0 = phase_sd_init
                                  ),
                                  body_params = body,
                                  n_segments = n_segments,
                                  acoustics = acoustics
                                )
                              }
  )
  # Add model parameters slot to scattering object =============================
  methods::slot(
    object,
    "model_parameters"
  )$SDWBA_curved <- list(
    parameters = stochastic_params,
    medium = medium_params,
    body = body
  )
  # Add model results slot to scattering object ================================
  methods::slot(
    object,
    "model"
  )$SDWBA_curved <- data.frame(
    frequency = frequency,
    f_bs = rep(
      NA,
      length(frequency)
    ),
    sigma_bs = rep(
      NA,
      length(frequency)
    ),
    TS = rep(
      NA,
      length(frequency)
    )
  )
  # Output =====================================================================
  return(object)
}

#' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' using the stochastic distorted Born wave approximation (SDWBA) model
#' @param object FLS-class scatterer.
#' @noRd
SDWBA_curved <- function(object) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract(object, "model_parameters")$SDWBA_curved
  body <- acousticTS::extract(object, "body")
  # Parse position matrices ====================================================
  theta <- model$body$theta
  # Material properties calculation ============================================
  g <- body$g
  h <- body$h
  R <- 1 / (g * h * h) + 1 / g - 2
  SDWBA_resampled_c <- function(i) {
    sub_params <- model$parameters[[i]]
    sdwba_resampled_backscatter(
      sub_params = sub_params,
      theta = theta,
      R = R,
      h = h
    )
  }
  # Generate results dataframe by collating resampled results ==================
  results <- do.call(
    "rbind",
    lapply(seq_along(model$parameters),
           FUN = function(i) SDWBA_resampled_c(i)
    )
  )
  # Update scatterer object ====================================================
  methods::slot(object, "model")$SDWBA_curved$f_bs <- results$f_bs
  methods::slot(object, "model")$SDWBA_curved$sigma_bs <- results$sigma_bs
  methods::slot(object, "model")$SDWBA_curved$TS <- results$TS_mean
  methods::slot(object, "model")$SDWBA_curved$TS_sd <- results$TS_sd
  # Return object ==============================================================
  object
}
