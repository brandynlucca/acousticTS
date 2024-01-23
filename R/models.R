################################################################################
# Primary scattering models for fluid-like scatterers (FLS)
################################################################################
#' Calculates the theoretical TS of a target using the deformed cylinder
#' model (DCM).
#'
#' @param object FLS-class scatterer.
#' @usage
#' DCM( object )
#' @return
#' Target strength (TS, dB re. 1 \eqn{m^2}) of a FLS-class object.
#' #' @references
#' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several
#' zooplankton groups. II. Scattering models. Journal of the Acoustical Society
#' of America, 103(1), 236-253.
#' @export
DCM <- function( object ) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract( object , "model_parameters" )$DCM
  body <- acousticTS::extract( object , "body" )
  body$density <- model$medium$density * body$g
  body$sound_speed <- model$medium$sound_speed * body$h
  # Multiply acoustic wavenumber by body radius ================================
  k1a <- model$parameters$acoustics$k_sw * model$body$radius
  k2a <- model$parameters$acoustics$k_b *model$ body$radius
  # Calculate Reflection Coefficient ===========================================
  R12 <- acousticTS::reflection_coefficient( model$medium , body )
  # Calculate Transmission Coefficients ========================================
  T12 <- 1 - R12 ; T21 <- 1 + R12
  # Calculate phase shift term for the ray model ===============================
  mu <- ( -pi / 2 * k1a ) / ( k1a + 0.4 )
  # Calculate Interference Term between echoes from front/back interfaces ======
  I0 <- 1 - T12 * T21 * base::exp( 4i * k2a ) * base::exp( 1i * mu * k1a )
  # Calculate half-width of directivity pattern ================================
  w12 <- base::sqrt( model$body$radius_curvature * model$body$radius)
  # Apply shift to orientation required for this model =========================
  theta_shift <- body$theta - pi / 2
  # Calculate Taylor expansion of directivity function =========================
  D_theta <- base::exp( - model$body$alpha_B *
                          ( 2 * theta_shift * model$body$radius_curvature /
                              model$body$length ) ^ 2 )
  # Calculate linear backscatter function, f_bs ================================
  f_bs <- 0.5 * w12 * R12 * base::exp( - 2i * k1a ) * I0 * D_theta
  # Update scatterer object ====================================================
  methods::slot( object, "model" )$DCM <- base::data.frame(
    frequency = model$parameters$acoustics$frequency ,
    ka = k1a ,
    f_bs = f_bs ,
    sigma_bs = base::abs( f_bs ) ^ 2 ,
    TS = 10 * base::log10( base::abs( f_bs ) ^ 2 )
  )
  # Return object ==============================================================
  return( object )
}
#' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' using the distorted Born wave approximation (DWBA) model.
#'
#' @param object FLS-class scatterer.
#' @references
#' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several
#' zooplankton groups. II. Scattering models. J. Acoust. Soc. Am., 103, 236-253.
#'
#' Demer, D.A., and Conti, S.G. 2003. Reconciling theoretical versus empirical
#' target strengths of krill: effects of phase variability on the distorted-wave
#' Born approximation. ICES J. Mar. Sci., 60: 429-434.
#' @import stats
#' @export
DWBA <- function( object ) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract( object , "model_parameters" )$DWBA
  body <- acousticTS::extract( object , "body" )
  theta <- body$theta
  r0 <- body$rpos[ 1 : 3 , ]
  # Material properties calculation ============================================
  g <- body$g ; h <- body$h
  R <- 1 / ( g * h * h ) + 1 / g - 2
  # Calculate rotation matrix and update wavenumber matrix =====================
  rotation_matrix <- matrix(  c( cos( theta ) ,
                                 0.0 ,
                                 sin( theta ) ) ,
                              1 )
  k_sw_rot <- model$parameters$acoustics$k_sw %*% rotation_matrix
  # Calculate Euclidean norms ==================================================
  k_sw_norm <- acousticTS::vecnorm( k_sw_rot )
  # Update position matrices  ==================================================
  rpos <- rbind( r0 , a = model$body$radius )
  # Calculate position matrix lags  ============================================
  rpos_diff <- t.default( diff( t.default( rpos ) ) )
  # Multiply wavenumber and body matrices ======================================
  rpos_diff_k <- t.default( sapply( 1 : length( k_sw_norm ) ,
                                    FUN = function( x ) {
                                      colSums( rpos_diff[ 1 : 3 , ] *
                                                 k_sw_rot[ x , ] )
                                    } ) )
  # Calculate Euclidean norms ==================================================
  rpos_diff_norm <- sqrt( colSums( rpos_diff[ 1 : 3 , ] * rpos_diff[ 1 : 3 , ] ) )
  # Estimate angles between body cylinders =====================================
  alpha <- acos( rpos_diff_k / ( k_sw_norm %*% t.default( rpos_diff_norm ) ) )
  beta <- abs( alpha - pi / 2 )
  # Define integrand ===========================================================
  integrand <- function( s , x ) {
    rint_mat <- s * rpos_diff + rpos[ , 1 : ncol( rpos_diff ) ]
    rint_k1_h_mat <- k_sw_rot[ x , ] %*% rint_mat[ 1 : 3 , ] / h
    bessel <- jc( 1 , 2 * ( k_sw_norm[ x ] * rint_mat[ 4 , ] / h *
                              cos( beta[ x , ] ) ) ) / cos( beta[ x , ] )
    fb_a <- k_sw_norm[ x ] / 4 * R * rint_mat[ 4 , ] * 
      exp( 2i * rint_k1_h_mat ) * bessel * rpos_diff_norm
    return( sum( fb_a , na.rm = T ) )
  }
  # Vectorize integrand function ===============================================
  integrand_vec <- Vectorize( integrand )
  # Calculate linear scatter response ==========================================
  f_bs <- rep( NA , length( k_sw_norm ) )
  f_bs <- sapply( 1 : length( k_sw_norm ) ,
                  FUN = function( x ) {
                    # Real ===============================================
                    Ri <- integrate( function( s ) {
                      Re( integrand_vec( s , x ) ) } , 0 , 1 )$value
                    # Real ===============================================
                    Ii <- integrate( function( s ) {
                      Im( integrand_vec( s , x ) ) } , 0 , 1 )$value
                    # Return =============================================
                    return( sqrt( Ri ^ 2 +  Ii ^ 2 ) ) } )
  # Update scatterer object ====================================================
  slot( object, "model" )$DWBA <- data.frame(
    frequency = model$parameters$acoustics$frequency ,
    ka = model$parameters$acoustics$k_sw * median( model$body$radius ) ,
    f_bs = f_bs ,
    sigma_bs = abs( f_bs ) * abs( f_bs ) ,
    TS = 10 * log10( abs( f_bs ) * abs( f_bs ) )
  )
  # Return object ==============================================================
  return( object )
}
#' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' using the distorted Born wave approximation (DWBA) model.
#'
#' @param object FLS-class scatterer.
#' @references
#' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several
#' zooplankton groups. II. Scattering models. J. Acoust. Soc. Am., 103, 236-253.
#'
#' Demer, D.A., and Conti, S.G. 2003. Reconciling theoretical versus empirical
#' target strengths of krill: effects of phase variability on the distorted-wave
#' Born approximation. ICES J. Mar. Sci., 60: 429-434.
#' @import stats
#' @export
DWBA_curved <- function( object ) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract( object , "model_parameters" )$DWBA_curved
  body <- acousticTS::extract( object , "body" )
  # Parse position matrices ====================================================
  r0 <- model$body$rpos[ 1 : 3 , ]
  theta <- model$body$theta
  # Material properties calculation ============================================
  g <- body$g ; h <- body$h
  R <- 1 / ( g * h * h ) + 1 / g - 2
  # Calculate rotation matrix and update wavenumber matrix =====================
  rotation_matrix <- matrix( c( cos( theta ) ,
                                0.0 ,
                                sin( theta ) ) ,
                             1 )
  k_sw_rot <- model$parameters$acoustics$k_sw %*% rotation_matrix
  # Calculate Euclidean norms ==================================================
  k_sw_norm <- acousticTS::vecnorm( k_sw_rot )
  # Update position matrices  ==================================================
  rpos <- rbind( r0 , a = model$body$radius )
  # Calculate position matrix lags  ============================================
  rpos_diff <- t.default( diff( t.default( rpos ) ) )
  # Multiply wavenumber and body matrices ======================================
  rpos_diff_k <- t.default( sapply( 1 : length( k_sw_norm ) ,
                                    FUN = function( x ) {
                                      colSums( rpos_diff[ 1 : 3 , ] *
                                                 k_sw_rot[ x , ] )
                                    } ) )
  # Calculate Euclidean norms ==================================================
  rpos_diff_norm <- sqrt( colSums( rpos_diff[ 1 : 3 , ] * rpos_diff[ 1 : 3 , ] ) )
  # Estimate angles between body cylinders =====================================
  alpha <- acos( rpos_diff_k / ( k_sw_norm %*% t.default( rpos_diff_norm ) ) )
  beta <- abs( alpha - pi / 2 )
  # Define integrand ===========================================================
  integrand <- function( s , x ) {
    rint_mat <- s * rpos_diff + rpos[ , 1 : ncol( rpos_diff ) ]
    rint_k1_h_mat <- k_sw_rot[ x , ] %*% rint_mat[ 1 : 3 , ] / h
    bessel <- jc( 1 , 2 * ( k_sw_norm[ x ] * rint_mat[ 4 , ] / h *
                              cos( beta[ x , ] ) ) ) / cos( beta[ x , ] )
    fb_a <- k_sw_norm[ x ] / 4 * R * rint_mat[ 4 , ] *
      exp( 2i * rint_k1_h_mat ) * bessel * rpos_diff_norm
    return( sum( fb_a , na.rm = T ) )
  }
  # Vectorize integrand function ===============================================
  integrand_vec <- Vectorize( integrand )
  # Calculate linear scatter response ==========================================
  f_bs <- rep( NA , length( k_sw_norm ) )
  f_bs <- sapply( 1 : length( k_sw_norm ) ,
                  FUN = function( x ) {
                    # Real ===============================================
                    Ri <- integrate( function( s ) {
                      Re( integrand_vec( s , x ) ) } , 0 , 1 )$value
                    # Real ===============================================
                    Ii <- integrate( function( s ) {
                      Im( integrand_vec( s , x ) ) } , 0 , 1 )$value
                    # Return =============================================
                    return( sqrt( Ri ^ 2 +  Ii ^ 2 ) ) } )
  # Update scatterer object ====================================================
  slot( object, "model" )$DWBA_curved <- data.frame(
    frequency = model$parameters$acoustics$frequency ,
    ka = model$parameters$acoustics$k_sw * median( model$body$radius ) ,
    f_bs = f_bs ,
    sigma_bs = abs( f_bs ) * abs( f_bs ) ,
    TS = 10 * log10( abs( f_bs ) * abs( f_bs ) )
  )
  # Return object ==============================================================
  return( object )
}
#' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' using the stochastic distorted Born wave approximation (DWBA) model.
#'
#' @param object FLS-class scatterer.
#' @references
#' Demer, D.A., and Conti, S.G. 2003. Reconciling theoretical versus empirical
#' target strengths of krill: effects of phase variability on the distorted-wave
#' Born approximation. ICES J. Mar. Sci., 60: 429-434.
#' @import stats
#' @export
SDWBA <- function( object ) {
  # Extract model parameters/inputs ============================================
  model <- extract( object , "model_parameters" )$SDWBA
  body <- extract( object , "body" )
  theta <- body$theta
  # Material properties calculation ============================================
  g <- body$g ; h <- body$h
  R <- 1 / ( g * h * h ) + 1 / g - 2
  # Calculate rotation matrix and update wavenumber matrix =====================
  rotation_matrix <- matrix( c( cos( theta ) ,
                                0.0 ,
                                sin( theta ) ) ,
                             1 )
  SDWBA_resampled <- function( i ) {
    # Parse phase-specific parameters ==========================================
    sub_params <- model$parameters[[ i ]]
    # Calculate rotation matrix and update wavenumber matrix ===================
    k_sw_rot <- sub_params$acoustics$k_sw %*% rotation_matrix
    # Calculate Euclidean norms ================================================
    k_sw_norm <- vecnorm( k_sw_rot )
    # Update position matrices  ================================================
    r0 <- rbind( sub_params$body_params$rpos[ 1 : 3 , ] , 
                 sub_params$body_params$radius )
    # Calculate position matrix lags  ==========================================
    r0_diff <- t.default( diff.default( t.default( r0 ) ) )
    # Multiply wavenumber and body matrices ====================================
    r0_diff_k <- t.default( sapply( 1 : length( k_sw_norm ) ,
                                    FUN = function( x ) {
                                      colSums( r0_diff[ 1 : 3 , ] *
                                                 k_sw_rot[ x , ] )
                                    } ) )
    # Calculate Euclidean norms ================================================
    r0_diff_norm <- sqrt( colSums( r0_diff[ 1 : 3 , ] * r0_diff[ 1 : 3 , ] ) )
    # Estimate angles between body cylinders ===================================
    alpha <- acos( r0_diff_k / ( k_sw_norm %*% t.default( r0_diff_norm ) ) )
    beta <- abs( alpha - pi / 2 )
    # Call in metrics ==========================================================
    phase_sd <- sub_params$meta_params$p0 
    r0_diff_h <- r0_diff / h 
    r0_h <- r0 / h
    # Define integrand =========================================================
    integrand <- function( s , x , y ) {
      rint_mat <- s * r0_diff_h[ , y ] + r0_h[ , y ]
      rint_k1_h_mat <- k_sw_rot[ x , ] %*% rint_mat[ 1 : 3 ]
      bessel <- jc( 1 , 2 * ( k_sw_norm[ x ] * rint_mat[ 4 ] *
                                cos( beta[ x , y ] ) ) ) / cos( beta[ x , y ] )
      fb_a <- k_sw_norm[ x ] / 4 * R * rint_mat[ 4 ] * h *
        exp( 2i * rint_k1_h_mat ) * bessel * r0_diff_norm[ y ]
      return( sum( fb_a , na.rm = T ) )
    } 
    # Vectorize integrand function =============================================
    integrand_vectorized <- Vectorize( integrand )
    stochastic_TS <- function( n_k , n_segments , n_iterations ) {
      sapply( 1 : n_k , 
              FUN = function( x ) {
                cyl_phase <- sapply( 1 :  n_segments ,
                                     FUN = function( y ) {
                                       phase_integrate( x , y , 
                                                        n_iterations , 
                                                        integrand_vectorized ,
                                                        phase_sd )
                                     }
                )
                cyl_sum_phase <- rowSums( cyl_phase , na.rm = T )
                cyl_sum_phase
              }
      ) -> phase_cyl
      data.frame( f_bs = colMeans( phase_cyl ) ,
                  sigma_bs = colMeans( sigma_bs( phase_cyl ) ) ,
                  TS_mean = db( colMeans( sigma_bs( phase_cyl ) ) ) ,
                  TS_sd = db( apply( sigma_bs( phase_cyl ) , 2 , sd ) ) )
    }
    # Calculate linear scatter response ========================================
    backscatter_df <- stochastic_TS( n_k = length( k_sw_norm ) , 
                                     n_segments = sub_params$n_segments - 1 , 
                                     n_iterations = sub_params$meta_params$n_iterations )
    return( backscatter_df )
  }
  # Generate results dataframe by collating resampled results ==================
  results <- do.call( "rbind" ,
                      lapply( 1 : length( model$parameters ) ,
                              FUN = function( i ) SDWBA_resampled( i )
                      ) )
  # Update scatterer object ====================================================
  slot( object , "model" )$SDWBA$f_bs <- results$f_bs
  slot( object , "model" )$SDWBA$sigma_bs <- results$sigma_bs
  slot( object , "model" )$SDWBA$TS <- results$TS_mean
  slot( object , "model" )$SDWBA$TS_sd <- results$TS_sd
  # Return object ==============================================================
  return( object )
}
#' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' using the stochastic distorted Born wave approximation (DWBA) model using 
#' Eq. (6) froom Stanton et al. (1998).
#' @param object FLS-class scatterer.
#' @references
#' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several
#' zooplankton groups. II. Scattering models. J. Acoust. Soc. Am., 103, 236-253.
#' @import stats
#' @rdname SDWBA_curved
#' @export
SDWBA_curved <- function( object ) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract( object , "model_parameters" )$SDWBA_curved
  body <- acousticTS::extract( object , "body" )
  # Parse position matrices ====================================================
  theta <- model$body$theta
  # Material properties calculation ============================================
  g <- body$g ; h <- body$h
  R <- 1 / ( g * h * h ) + 1 / g - 2
  # Calculate rotation matrix and update wavenumber matrix =====================
  rotation_matrix <- matrix( c( cos( theta ) ,
                                0.0 ,
                                sin( theta ) ) ,
                             1 )
  SDWBA_resampled_c <- function( i ) {
    # Parse phase-specific parameters ==========================================
    sub_params <- model$parameters[[ i ]]
    # Calculate rotation matrix and update wavenumber matrix ===================
    k_sw_rot <- sub_params$acoustics$k_sw %*% rotation_matrix
    # Calculate Euclidean norms ================================================
    k_sw_norm <- vecnorm( k_sw_rot )
    # Update position matrices  ================================================
    r0 <- rbind( sub_params$body_params$rpos[ 1 : 3 , ] , 
                 sub_params$body_params$radius )
    # Calculate position matrix lags  ==========================================
    r0_diff <- t.default( diff.default( t.default( r0 ) ) )
    # Multiply wavenumber and body matrices ====================================
    r0_diff_k <- t.default( sapply( 1 : length( k_sw_norm ) ,
                                    FUN = function( x ) {
                                      colSums( r0_diff[ 1 : 3 , ] *
                                                 k_sw_rot[ x , ] )
                                    } ) )
    # Calculate Euclidean norms ================================================
    r0_diff_norm <- sqrt( colSums( r0_diff[ 1 : 3 , ] * r0_diff[ 1 : 3 , ] ) )
    # Estimate angles between body cylinders ===================================
    alpha <- acos( r0_diff_k / ( k_sw_norm %*% t.default( r0_diff_norm ) ) )
    beta <- abs( alpha - pi / 2 )
    # Call in metrics ==========================================================
    phase_sd <- sub_params$meta_params$p0 
    r0_diff_h <- r0_diff / h 
    r0_h <- r0 / h
    # Define integrand =========================================================
    integrand_c <- function( s , x , y ) {
      rint_mat <- s * r0_diff_h[ , y ] + r0_h[ , y ]
      rint_k1_h_mat <- k_sw_rot[ x , ] %*% rint_mat[ 1 : 3 ]
      bessel <- jc( 1 , 2 * ( k_sw_norm[ x ] * rint_mat[ 4 ] *
                                cos( beta[ x , y ] ) ) ) / cos( beta[ x , y ] )
      fb_a <- k_sw_norm[ x ] / 4 * R * rint_mat[ 4 ] * h *
        exp( 2i * rint_k1_h_mat ) * bessel * r0_diff_norm[ y ]
      return( sum( fb_a , na.rm = T ) )
    } 
    # Vectorize integrand function =============================================
    integrand_vectorized_c <- Vectorize( integrand_c )
    stochastic_TS_c <- function( n_k , n_segments , n_iterations ) {
      sapply( 1 : n_k , 
              FUN = function( x ) {
                cyl_phase <- sapply( 1 :  n_segments ,
                                     FUN = function( y ) {
                                       phase_integrate( x , y , 
                                                        n_iterations , 
                                                        integrand_vectorized_c ,
                                                        phase_sd )
                                     }
                )
                cyl_sum_phase <- rowSums( cyl_phase , na.rm = T )
                cyl_sum_phase
              }
      ) -> phase_cyl
      data.frame( f_bs = colMeans( phase_cyl ) ,
                  sigma_bs = colMeans( sigma_bs( phase_cyl ) ) ,
                  TS_mean = db( colMeans( sigma_bs( phase_cyl ) ) ) ,
                  TS_sd = db( apply( sigma_bs( phase_cyl ) , 2 , sd ) ) )
    }
    # Calculate linear scatter response ========================================
    backscatter_df <- stochastic_TS_c( n_k = length( k_sw_norm ) , 
                                       n_segments = sub_params$n_segments - 1 , 
                                       n_iterations = sub_params$meta_params$n_iterations )
    return( backscatter_df )
  }
  # Generate results dataframe by collating resampled results ==================
  results <- do.call( "rbind" ,
                      lapply( 1 : length( model$parameters ) ,
                              FUN = function( i ) SDWBA_resampled_c( i )
                      ) )
  # Update scatterer object ====================================================
  slot( object , "model" )$SDWBA_curved$f_bs <- results$f_bs
  slot( object , "model" )$SDWBA_curved$sigma_bs <- results$sigma_bs
  slot( object , "model" )$SDWBA_curved$TS <- results$TS_mean
  slot( object , "model" )$SDWBA_curved$TS_sd <- results$TS_sd
  # Return object ==============================================================
  return( object )
}
################################################################################
# Primary scattering model for an elastic solid sphere (CAL)
################################################################################
#' Calculates theoretical TS of a solid sphere of a certain material at a given
#' frequency.
#'
#' @description
#' This function is a wrapper around TS_calculate(...) that parametrizes the
#' remainder of the model, while also doing simple calculations that do not
#' need to be looped. This function provides a TS estimate at a given frequency.
#' @param object CAL-class object.
#' @usage
#' calibration(object)
#' @return The theoretical acoustic target strength (TS, dB re. 1 \eqn{m^2})  of a
#' solid sphere at a given frequency.
#' @references
#' MacLennan D. N. (1981). The theory of solid spheres as sonar calibration
#' targets. Scottish Fisheries Research No. 22, Department of Agriculture and
#' Fisheries for Scotland.
#' @export
calibration <- function( object ) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract( object ,
                                "model_parameters" )$calibration
  ### Now we solve / calculate equations =======================================
  # Equations 6a -- weight wavenumber by radius ================================
  ka_sw <- model$parameters$acoustics$k_sw * model$body$radius
  ka_l <- model$parameters$acoustics$k_l * model$body$radius
  ka_t <- model$parameters$acoustics$k_t * model$body$radius
  # Set limit for iterations ===================================================
  m_limit <- round( max( ka_sw ) ) + 10
  # Create modal series number vector ==========================================
  m <- 0 : m_limit
  # Convert these vectors into matrices ========================================
  ka_sw_m <- modal_matrix( ka_sw , m_limit )
  ka_l_m <- modal_matrix( ka_l , m_limit )
  ka_t_m <- modal_matrix( ka_t , m_limit )
  # Calculate Legendre polynomial ==============================================
  Pl <- Pn( m , cos( model$body$theta ) )
  # Calculate spherical Bessel functions of first kind =========================
  js_mat <- js( m , ka_sw_m )
  js_mat_l <- js( m , ka_l_m )
  js_mat_t <- js( m , ka_t_m )
  # Calculate spherical Bessel functions of second kind ========================
  ys_mat <- ys( m , ka_sw_m )
  # Calculate first derivative of spheric Bessel functions of first kind =======
  jsd_mat <- jsd( m , ka_sw_m )
  jsd_mat_l <- jsd( m , ka_l_m )
  jsd_mat_t <- jsd( m , ka_t_m )
  # Calculate first derivative of spheric Bessel functions of second kind ======
  ysd_mat <- ysd( m , ka_sw_m )
  # Calculate density contrast =================================================
  g <- model$body$density / model$medium$density
  # Tangent functions ==========================================================
  tan_sw <- - ka_sw_m * jsd_mat / js_mat
  tan_l <- - ka_l_m * jsd_mat_l / js_mat_l
  tan_t <- - ka_t_m * jsd_mat_t / js_mat_t
  tan_beta <- - ka_sw_m * ysd_mat / ys_mat
  tan_diff <- - js_mat / ys_mat
  # Difference terms ===========================================================
  along_m <- ( m * m + m )
  tan_l_add <- tan_l + 1
  tan_t_div <- along_m - 1 - ka_t_m * ka_t_m / 2 + tan_t
  numerator <- ( tan_l / tan_l_add ) - ( along_m / tan_t_div )
  denominator1 <- ( along_m - ka_t_m * ka_t_m / 2 + 2 * tan_l ) / tan_l_add
  denominator2 <- along_m * ( tan_t + 1 ) / tan_t_div
  denominator <- denominator1 - denominator2
  ratio <- -0.5 * ( ka_t_m * ka_t_m ) * numerator / denominator
  # Additional trig functions ==================================================
  phi <- - ratio / g
  eta_tan <- tan_diff * ( phi + tan_sw ) / ( phi + tan_beta )
  cos_eta <- 1 / sqrt( 1 + eta_tan * eta_tan )
  sin_eta <- eta_tan * cos_eta
  # Fill in rest of Hickling (1962) equation ===================================
  f_j <- colSums(
    ( 2 * m + 1 ) * Pl[ m + 1 ] * ( sin_eta * ( 1i * cos_eta - sin_eta ) )
  )
  # Calculate linear backscatter coefficient ===================================
  f_bs <- abs( - 2i * f_j / ka_sw ) * model$body$radius / 2
  sigma_bs <- f_bs * f_bs
  TS <- 10 * log10( sigma_bs )
  # Add results to scatterer object ============================================
  slot( object ,
        "model" )$calibration <- data.frame(
          frequency = model$parameters$acoustics$frequency ,
          ka = ka_sw ,
          f_bs = f_bs ,
          sigma_bs = sigma_bs ,
          TS = TS
        )
  # Return object ==============================================================
  return( object )
}
################################################################################
################################################################################
# GAS-FILLED SCATTERERS
################################################################################
################################################################################
# Anderson (1950) gas-filled fluid-sphere model
################################################################################
#' Calculates the theoretical TS of a fluid sphere using an exact modal series
#' solution proposed by Andersen (1950).
#' @param object GAS- or SBF-class object.
#' @details
#' Calculates the theoretical TS of a fluid sphere at a given frequency using
#' an exact modal series solution.
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Anderson, V.C. (1950). Sound scattering from a fluid sphere. Journal of the
#' Acoustical Society of America, 22, 426-431.
#' @export
MSS_anderson <- function( object ) {
  # Extract model parameters/inputs ============================================
  model <- extract( object , "model_parameters" )$MSS_anderson
  # Multiple acoustic wavenumber by radius =====================================
  ## Medium ====================================================================
  k1a <- model$parameters$acoustics$k_sw * model$body$radius
  ## Fluid within sphere =======================================================
  k2a <- model$parameters$acoustics$k_f * model$body$radius
  # Set limit for iterations ===================================================
  m_limit <- model$parameters$ka_limit
  m <- 0 : m_limit
  # Convert ka vectors to matrices =============================================
  ka1_m <- modal_matrix( k1a , m_limit )
  ka2_m <- modal_matrix( k2a , m_limit )
  # Calculate modal series coefficient, b_m (or C_m) ===========================
  ## Material properties term ==================================================
  gh <- model$body$g * model$body$h
  # Numerator term =============================================================
  N1 <- ( jsd( m , ka2_m ) * ys( m , ka1_m ) ) /
    ( js( m , ka2_m ) * jsd( m , ka1_m ) )
  N2 <- ( gh * ysd( m , ka1_m ) / jsd( m , ka1_m ) )
  CN <- N1 - N2
  # Denominator term ===========================================================
  D1 <- ( jsd( m , ka2_m ) * js( m , ka1_m ) ) /
    ( js( m , ka2_m ) * jsd( m , ka1_m ) )
  D2 <- gh
  CD <- D1 - D2
  # Finalize modal series coefficient ==========================================
  C_m <- CN / CD
  b_m <- -1 / ( 1 + 1i * C_m )
  # Calculate linear scatter response ==========================================
  ## Sum across columns to complete modal series summation over frequency range=
  f_j <- colSums( ( 2 * m + 1 ) * ( -1 ) ^ m * b_m )
  f_sphere <- -1i / model$parameters$acoustics$k_sw * f_j
  # Calculate linear scattering coefficient, sigma_bs ==========================
  sigma_bs <- abs( f_sphere ) * abs( f_sphere )
  # Calculate TS ===============================================================
  TS <- 10 * log10( sigma_bs )
  # Add results to scatterer object ============================================
  slot( object ,
        "model" )$MSS_anderson <- data.frame(
          frequency = model$parameters$acoustics$frequency ,
          ka = k1a ,
          f_bs = f_sphere ,
          sigma_bs = sigma_bs ,
          TS = TS )
  # Return object ==============================================================
  return( object )
}
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
KRM <- function( object ) {
  # Detect object class ========================================================
  scatterer_type <- class( object )
  # Extract model parameter inputs =============================================
  model <- extract( object , "model_parameters" )$KRM
  # Extract body parameters ====================================================
  body <- extract( object , "body" )
  # Calculate reflection coefficient for medium-body interface =================
  R12 <- reflection_coefficient( model$medium , model$body )
  # Calculate transmission coefficient and its reverse =========================
  T12T21 <- 1 - R12  * R12
  # Sum across body position vector ============================================
  rpos <- switch( scatterer_type ,
                  FLS = rbind( x = body$rpos[ 1 , ] ,
                               w = c( body$radius[ 2 ] ,
                                      body$radius[ 2 : ( length( body$radius ) - 1 ) ] ,
                                      body$radius[ ( length( body$radius ) ) - 1 ] ) * 2 ,
                               zU =c( body$radius[ 2 ] ,
                                      body$radius[ 2 : ( length( body$radius ) - 1 ) ] ,
                                      body$radius[ ( length( body$radius ) ) - 1 ] ) ,
                               zL = - c( body$radius[ 2 ] ,
                                         body$radius[ 2 : ( length( body$radius ) - 1 ) ] ,
                                         body$radius[ ( length( body$radius ) ) - 1 ] ) ) ,
                  SBF = body$rpos )
  body_rpos_sum <- along_sum( rpos , model$parameters$ns_b )
  # Approximate radius of body cylinders =======================================
  a_body <- switch( scatterer_type ,
                    FLS = body_rpos_sum[ 2 , ] / 4 ,
                    SBF = body_rpos_sum[ 2 , ] / 4 )
  # Combine wavenumber (k) and radii to calculate "ka" =========================
  ka_body <- matrix( data = rep( a_body ,
                                 each = length( model$parameters$acoustics$k_sw ) ) ,
                     ncol = length( a_body ) ,
                     nrow = length( model$parameters$acoustics$k_b ) ) * model$parameters$acoustics$k_b
  # Convert c-z coordinates to required u-v rotated coordinates ================
  uv_body <- acousticTS::body_rotation( body_rpos_sum ,
                                        body$rpos ,
                                        body$theta ,
                                        length( model$parameters$acoustics$k_sw ) )
  # Calculate body empirical phase shift function ==============================
  body_dorsal_sum <- switch( scatterer_type ,
                             FLS = matrix( data = rep( body_rpos_sum[ 3 , ] ,
                                                       each = length( model$parameters$acoustics$k_sw ) ) ,
                                           ncol = length( body_rpos_sum[ 3 , ] ) ,
                                           nrow = length( model$parameters$acoustics$k_sw ) ) / 2 ,
                             SBF = matrix( data = rep( body_rpos_sum[ 3 , ] ,
                                                       each = length( model$parameters$acoustics$k_sw ) ) ,
                                           ncol = length( body_rpos_sum[ 3 , ] ) ,
                                           nrow = length( model$parameters$acoustics$k_sw ) ) / 2 )
  Psi_b <- - pi * model$parameters$acoustics$k_b * body_dorsal_sum /
    (2 * ( model$parameters$acoustics$k_b * body_dorsal_sum + 0.4 ) )
  # Estimate natural log function (phase, etc.) ================================
  exp_body <- exp( - 2i * model$parameters$acoustics$k_sw * uv_body$vbU ) -
    T12T21 * exp( - 2i * model$parameters$acoustics$k_sw * uv_body$vbU +
                    2i * model$parameters$acoustics$k_b *
                    ( uv_body$vbU - uv_body$vbL ) + 1i * Psi_b )
  # Resolve summation term =====================================================
  body_summation <- sqrt( ka_body ) * uv_body$delta_u
  # Calculate linear scattering length (m) =====================================
  f_body <- rowSums( - ( ( 1i * ( R12 / ( 2 * sqrt( pi ) ) ) ) *
                           body_summation * exp_body ) )
  if ( scatterer_type == "FLS" ) {
    # Define KRM slot for FLS-type scatterer ===================================
    slot( object , "model" )$KRM <- data.frame( frequency = model$parameters$acoustics$frequency ,
                                                ka = model$parameters$acoustics$k_sw * median( a_body , na.rm = T ) ,
                                                f_bs = f_body ,
                                                sigma_bs = abs( f_body ) * abs( f_body ) ,
                                                TS = 20 * log10( abs( f_body ) ) )
  } else if ( scatterer_type == "SBF" ) {
    #### Repeat process for bladder ============================================
    # Extract bladder parameters ===============================================
    bladder <- acousticTS::extract( object , "bladder" )
    # Calculate reflection coefficient for bladder =============================
    R23 <- acousticTS::reflection_coefficient( body ,
                                               bladder )
    # Sum across body/swimbladder position vectors =============================
    bladder_rpos_sum <- acousticTS::along_sum( bladder$rpos ,
                                               model$parameters$ns_sb )
    # Approximate radii of bladder discrete cylinders ==========================
    a_bladder <- bladder_rpos_sum[ 2 , ] / 4
    # Combine wavenumber (k) and radii to calculate "ka" for bladder ===========
    ka_bladder <- matrix( data = rep( a_bladder , each = length( model$parameters$acoustics$k_sw ) ) ,
                          ncol = length( a_bladder ) ,
                          nrow = length( model$parameters$acoustics$k_sw ) ) * model$parameters$acoustics$k_sw
    # Calculate Kirchoff approximation empirical factor, A_sb ==================
    A_sb <- ka_bladder / ( ka_bladder + 0.083 )
    # Calculate empirical phase shift for a fluid cylinder, Psi_p ==============
    Psi_p <- ka_bladder / ( 40 + ka_bladder ) - 1.05
    # Convert x-z coordinates to requisite u-v rotated coordinates =============
    uv_bladder <- acousticTS::bladder_rotation( bladder_rpos_sum ,
                                                bladder$rpos ,
                                                bladder$theta ,
                                               length( model$parameters$acoustics$k_sw ) )
    # Estimate natural log functions  ==========================================
    exp_bladder <- exp( - 1i * ( 2 * model$parameters$acoustics$k_b *
                                   uv_bladder$v + Psi_p ) ) * uv_bladder$delta_u
    # Calculate the summation term =============================================
    bladder_summation <- A_sb * sqrt( ( ka_bladder + 1 ) * sin( bladder$theta ) )
    # Estimate backscattering length, f_fluid/f_soft ===========================
    f_bladder <- rowSums( - 1i * ( R23 * T12T21 ) / ( 2 * sqrt( pi ) ) *
                            bladder_summation * exp_bladder )
    # Estimate total backscattering length, f_bs ===============================
    f_bs <- f_body + f_bladder
    # Define KRM slot for FLS-type scatterer ===================================
    slot( object , "model" )$KRM <- data.frame( frequency = model$parameters$acoustics$frequency ,
                                                               ka = model$parameters$acoustics$k_sw * median( a_body , na.rm = T ) ,
                                                               f_body = f_body ,
                                                               f_bladder = f_bladder ,
                                                               f_bs = f_bs ,
                                                               sigma_bs = abs( f_bs ) * abs( f_bs ) ,
                                                               TS = 20 * log10( abs( f_bs ) ) )
  }
  # Return object ==============================================================
  return( object )
}
################################################################################
# Primary scattering model for an elastic shelled scatterers (ESS)
################################################################################
################################################################################
# Ray-based high pass approximation
################################################################################

#' Calculates the theoretical TS of a shelled organism using the non-modal
#' High Pass (HP) model
#'
#' @param object Desired animal object (Elastic Shelled).
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Lavery, A.C., Wiebe, P.H., Stanton, T.K., Lawson, G.L., Benfield, M.C.,
#' Copley, N. 2007. Determining dominant scatterers of sound in mixed
#' zooplankton popuilations. The Journal of the Acoustical Society of America,
#' 122(6): 3304-3326.
#'
#' @export
high_pass_stanton <- function( object ) {
    # Extract model parameters/inputs ==========================================
    model_params <- extract( object , "model_parameters" )$high_pass_stanton
    medium <- model_params$medium
    acoustics <- model_params$acoustics
    shell <- model_params$shell
    # Multiply acoustic wavenumber by body radius ==============================
    k1a <- acoustics$k_sw * shell$radius
    # Calculate Reflection Coefficient =========================================
    R <- ( shell$h * shell$g - 1 ) / ( shell$h * shell$g + 1 )
    # Calculate backscatter constant, alpha_pi
    alpha_pi <- (1 - shell$g * ( shell$h * shell$h ) ) / ( 3 * shell$g * ( shell$h * shell$h ) ) +
      ( 1 - shell$g ) / ( 1 + 2 * shell$g )
    # Define approximation constants, G_c and F_c ================================
    F_c <- 1; G_c <- 1
    # Caclulate numerator term ===================================================
    num <- ( ( shell$radius * shell$radius ) * ( k1a * k1a * k1a * k1a ) * 
               ( alpha_pi * alpha_pi ) * G_c)
    # Caclulate denominator term =================================================
    dem <- ( 1 + ( 4 * ( k1a * k1a * k1a * k1a ) * ( alpha_pi * alpha_pi ) ) / 
               ( ( R * R ) * F_c ) )
    # Calculate backscatter and return
    f_bs <- num / dem
    slot( object , "model" )$high_pass_stanton <- data.frame( f_bs = f_bs ,
                                                              sigma_bs = abs( f_bs ),
                                                              TS = 10 * log10( abs( f_bs ) ) )
      return( object )
}