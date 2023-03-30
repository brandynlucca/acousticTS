################################################################################
# Generic for "plot(...)" for each scattering class object
################################################################################
#' Method for what is printed for objects.
#'
#' @param x Scatterer-class object.
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param ... Additional plot inputs
#' @export
setMethod( f = "plot" ,
          methods::signature( x = "scatterer" ) ,
          definition = function( x ,
                                 type = "shape" ,
                                 nudge_y = 1.1 ,
                                 nudge_x = 1.05 ,
                                 aspect_ratio = "manual" ,
                                 x_units = "frequency",
                                 y_units = "TS" , ... ) {
            # Detect scatterer type ============================
            sc_type <- base::class( x )

            switch( sc_type,
                    CAL = cal_plot( x , type , nudge_y , nudge_x , x_units , ... ) ,
                    SBF = sbf_plot( x , type , nudge_y , nudge_x , x_units , ... ) ,
                    FLS = fls_plot( x , type , nudge_y , nudge_x , aspect_ratio , x_units , y_units ) ,
                    GAS = gas_plot( x , type , nudge_y , nudge_x , x_units , ... ) )
          } )
################################################################################
# Methods for "plot(...)" for each scattering class object
################################################################################
#' Base plotting color palette
#' @export
plot_colors <- base::c( "black" ,
                        "azure4" ,
                        "firebrick3" ,
                        "orangered3" ,
                        "sienna3" ,
                        "royalblue4" ,
                        "steelblue3" ,
                        "darkolivegreen" ,
                        "cadetblue" )
################################################################################
#' Plotting for CAL-class objects
#' @param object CAL-class object.
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param ... Additional plot inputs
#' @export
#' @import graphics
#' @export
cal_plot <- function( object ,
                      type = "shape" ,
                      nudge_y = 1.01 ,
                      nudge_x = 1.01 ,
                      x_units = "frequency" , ...) {
  # Retrieve default plot window parameters ====================================
  opar <- par( no.readonly = TRUE )
  on.exit( par( opar ) )
  if( type == "shape" ) {
    # Extract body shape information ===========================================
    body <- acousticTS::extract( object ,
                                 "body" )
    # Define plot margins ======================================================
    par( ask = FALSE ,
         oma = base::c( 1 , 1 , 1 , 0 ) ,
         mar = base::c( 3 , 4.5 , 1 , 2 ) )
    # Begin plotting ===========================================================
    plot( x = body$rpos[ , 1 ] ,
          y = body$rpos[ , 2 ] ,
          type = 'l' ,
          ylab = "Semi-minor diameter (m)" ,
          xlab = "Semi-major dimaeter (m)" ,
          lwd = 4 ,
          cex.lab = 1.2 ,
          cex.axis = 1.2 ,
          ylim = base::c( min( body$rpos[ , 3 ] ) * ( 1 - (1 - nudge_y ) ) ,
                          max( -body$rpos[ , 3 ] ) * nudge_y ) )
    # Add lower perimeter of shape =============================================
    lines( x = body$rpos[ , 1 ] ,
           y = body$rpos[ , 3 ] ,
           lty = 1 ,
           lwd = 4 )
    # Add body segments ========================================================
    segments( x0 = body$rpos[ , 1 ] ,
              x1 = body$rpos[ , 1 ] ,
              y0 = body$rpos[ , 2 ] ,
              y1 = body$rpos[ , 3 ] ,
              lty = 3 ,
              lwd = 1.25 )
  } else if (type == "model") {
    if(length(extract(object, "model")) == 0) {
      stop("ERROR: no model results detected in object.")
    } else {
      # Extract body shape information ============================
      shape <- extract( object , "body" )
      # Extract model results ====================================
      TS <- extract( object , "model" )$calibration$TS
      x_axis_domain <- extract(object, "model_parameters")$calibration$parameters
      if( x_units == "frequency" ) {
        x_axis <- x_axis_domain$acoustics$frequency * 1e-3
        x_lab <- "Frequency (kHz)"
      } else if (x_units == "k_sw") {
        x_axis <- x_axis_domain$acoustics$k_sw * shape$radius
        x_lab <- expression(italic(k[sw]*a))
      } else if (x_units == "k_l") {
        x_axis <- x_axis_domain$acoustics$k_l * shape$radius
        x_lab <- expression(italic(k[l]*a))
      } else{
        x_axis <- x_axis_domain$acoustics$k_t * shape$radius
        x_lab <- expression(italic(k[t]*a))
      }
      # Plot results ===============================================
      par(ask = F,
          mar = c(4, 4.5, 1, 1))
      plot(x = x_axis,
           y = TS,
           type = 'l',
           xlab = x_lab,
           ylab = expression(Target~strength~(dB~re.~1~m^2)),
           lwd = 2.5,
           cex.lab = 1.2,
           cex.axis = 1.2,
           xlim = c( base::abs( min( x_axis ) * ( 1 - ( 1 + nudge_x ) ) ) ,
                     max(x_axis) * (nudge_x)),
           ylim = c(min(TS) * (1 - (1 - nudge_y)),
                    max(TS) * (1 + (1 - nudge_y))),
           xaxs = "i",
           yaxs = "i")
    }
  }
  invisible()
}
################################################################################
# Methods for "plot(...)" for each scattering class object
################################################################################
#' Plotting for FLS-class objects
#' @param object FLS-class object.
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param ... Additional plot inputs
#' @import graphics
#' @import stats
#' @import grDevices
#' @export
fls_plot <- function( object,
                      type = "shape" ,
                      nudge_y = 1.05 ,
                      nudge_x = 1.01 ,
                      aspect_ratio = "manual" ,
                      x_units = "frequency" ,
                      y_units = "TS" , ... ) {
  # Retrieve default plot window parameters ====================================
  opar <- graphics::par( no.readonly = TRUE )
  base::on.exit( graphics::par( opar ) )
  if( type == "shape" ) {
    # Extract body shape information ===========================================
    body <- acousticTS::extract( object ,
                                 "body" )
    # Define plot margins ======================================================
    par( ask = FALSE ,
         oma = base::c( 1 , 1 , 1 , 0 ) ,
         mar = base::c( 5.0 , 4.5 , 1.5 , 2 ) )
    # Center shape =============================================================
    body$rpos[ 3 , ] <- body$rpos[ 3 , ] - stats::median( body$rpos[ 3 , ] )
    # Adjust axes ==============================================================
    if ( aspect_ratio == "manual" ) {
      vert_lims <- base::c( min( body$rpos[ 3 , ] - body$radius ) * ( 1 - ( 1 - nudge_y ) ) ,
                            max( body$rpos[ 3 , ] + body$radius ) * nudge_y )
    } else {
      vert_lims <- base::c( - base::max( body$rpos[ 1 , ]  )  * 0.10 ,
                            base::max( body$rpos[ 1 ,  ]  )  * 0.10 )
    }
    # Begin plotting ===========================================================
    graphics::plot( x = body$rpos[ 1 , ] ,
                    y = body$rpos[ 3 , ] ,
                    type = 'l' ,
                    lwd = 4 ,
                    cex.lab = 1.2 ,
                    cex.axis = 1.2 ,
                    xlab = "Length (mm)" ,
                    ylab = "Thickness (mm)" ,
                    ylim = vert_lims )
    # Add lower perimeter of shape =============================================
    graphics::lines( x = body$rpos[ 1 , ] ,
                     y = body$rpos[ 3 , ] - body$radius ,
                     lty = 1 ,
                     lwd = 4 )
    # Add upper perimeter of shape =============================================
    graphics::lines( x = body$rpos[ 1 , ] ,
                     y = body$rpos[ 3 , ] + body$radius ,
                     lty = 1 ,
                     lwd = 4 )
    # Add body segments ========================================================
    graphics::segments( x0 = body$rpos[ 1 , ] ,
                        x1 = body$rpos[ 1 , ] ,
                        y0 = body$rpos[ 3 , ] - body$radius ,
                        y1 = body$rpos[ 3 , ] + body$radius ,
                        lty = 3 ,
                        lwd = 1.25 )
  } else if ( type == "model" ) {
    # Detect model selection ===================================================
    models <- acousticTS::extract( object , "model" )
    model_names <- base::names( models )
    if ( base::length( model_names ) == 0 )
      stop( "ERROR: no model results detected in object." )
    # Extract body shape information ===========================================
    shape <- acousticTS::extract( object , "body" )
    # Append model name ========================================================
    models <- base::lapply( 1 : base::length( model_names ) ,
                            FUN = function( x ) {
                              base::transform( models[[x]] ,
                                               model = model_names[x] ) } )
    # Convert into a data.frame ================================================
    models_df <- base::do.call( "rbind" , models )
    # Define x-axis domain =====================================================
    # x_axis <- base::switch( x_units ,
    #                         frequency = base::unique( models_df$frequency ) * 1e-3 ,
    #                         k_sw = base::unique( models_df$ka ) )
    x_axis <- models_df[ , base::which( base::colnames( models_df )  == x_units ) ]
    x_mat <- base::split( x_axis , models_df$model )
    y_axis <- models_df[ , base::which( base::colnames( models_df ) == y_units ) ]
    y_mat <- base::split( y_axis , models_df$model )
    col_axis <- plot_colors[ base::as.numeric( base::as.factor( models_df$model ) ) ]
    x_lab <- base::switch( x_units ,
                           frequency = "Frequency (Hz)" ,
                           k_sw = base::expression(italic(k[sw]*a) ) )
    # Define plot margins ======================================================
    graphics::par( ask = FALSE ,
                   mar = base::c( 4.0 , 4.5 , 1.0 , 1.0 ) )
    # Initiate plotting ========================================================
    graphics::plot( x = base::seq( from = base::min( x_axis ) ,
                                   to = base::max( x_axis ) ,
                                   length.out = 2 ) ,
                    y = base::seq( from = base::min( y_axis ) ,
                                   to = base::max( y_axis ) ,
                                   length.out = 2 ) ,
                    xlim = base::c( base::min( x_axis ) * ( 1 - nudge_x ),
                                    base::max( x_axis ) * ( nudge_x ) ) ,
                    ylim = base::c( base::min( y_axis ) * ( 1 - ( 1 - nudge_y ) ) ,
                                    base::max( y_axis ) * ( 1 + ( 1 - nudge_y ) ) ) ,
                    xlab = x_lab ,
                    ylab = expression( "Target"~"strength"~("dB"~"re."~1~"m"^2) ) ,
                    yaxs = "i" ,
                    xaxs = "i" ,
                    type = 'n' )
   base::invisible( base::mapply( lines , x_mat , y_mat ,
                                  col = plot_colors[ 1 : base::length( model_names ) ] ,
                                  lwd = 2 ) )
   graphics::legend( "topright" ,
                     title = base::expression( bold("TS"~"model") ) ,
                     title.adj = 0.05 ,
                     legend = model_names ,
                     lty = base::rep( 1 , base::length( model_names ) ) ,
                     col = plot_colors[ 1 : base::length( model_names ) ] ,
                     cex = 1.05 )
   }
  base::invisible( )
}
#' Plotting for GAS-class objects
#' @param object GAS-class object.
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param ... Additional plot inputs
#' @import graphics
#' @import stats
#' @import grDevices
#' @export
gas_plot <- function( object ,
                      type = "shape" ,
                      nudge_y = 1.00 ,
                      nudge_x = 1.00 ,
                      x_units = "frequency" , ...) {
  # Retrieve default plot window parameters ====================================
  opar <- graphics::par( no.readonly = TRUE )
  base::on.exit( graphics::par( opar ) )
  if( type == "shape" ) {
    # Extract body shape information ===========================================
    body <- acousticTS::extract( object ,
                                 "body" )
    # Define plot margins ======================================================
    par( ask = FALSE ,
         oma = base::c( 1 , 1 , 1 , 0 ) ,
         mar = base::c( 3 , 4.5 , 1 , 2 ) )
    # Begin plotting ===========================================================
    plot( x = body$rpos[ , 1 ] ,
          y = body$rpos[ , 2 ] ,
          type = 'l' ,
          ylab = "Semi-minor diameter (m)" ,
          xlab = "Semi-major dimaeter (m)" ,
          lwd = 4 ,
          cex.lab = 1.2 ,
          cex.axis = 1.2 ,
          ylim = base::c( min( body$rpos[ , 3 ] ) * ( 1 - (1 - nudge_y ) ) ,
                          max( -body$rpos[ , 3 ] ) * nudge_y ) )
    # Add lower perimeter of shape =============================================
    lines( x = body$rpos[ , 1 ] ,
           y = body$rpos[ , 3 ] ,
           lty = 1 ,
           lwd = 4 )
    # Add body segments ========================================================
    segments( x0 = body$rpos[ , 1 ] ,
              x1 = body$rpos[ , 1 ] ,
              y0 = body$rpos[ , 2 ] ,
              y1 = body$rpos[ , 3 ] ,
              lty = 3 ,
              lwd = 1.25 )
  } else if ( type == "model" ) {
    # Detect model selection ===================================================
    models <- acousticTS::extract( object , "model" )
    model_names <- base::names( models )
    if ( base::length( model_names ) == 0 )
      stop( "ERROR: no model results detected in object." )
    # Extract body shape information ===========================================
    shape <- acousticTS::extract( object , "body" )
    # Append model name ========================================================
    models <- base::lapply( 1 : base::length( model_names ) ,
                            FUN = function( x ) {
                              base::transform( models[[x]] ,
                                               model = model_names[x] ) } )
    # Convert into a data.frame ================================================
    models_df <- base::do.call( "rbind" , models )
    # Define x-axis domain =====================================================
    x_axis <- base::switch( x_units ,
                            frequency = base::unique( models_df$frequency ) * 1e-3 ,
                            k_sw = base::unique( models_df$ka ) )
    x_lab <- base::switch( x_units ,
                           frequency = "Frequency (kHz)" ,
                           k_sw = base::expression(italic("k"["sw"]*"a") ) )
    # Define plot margins ======================================================
    graphics::par( ask = FALSE ,
                   mar = base::c( 4.0 , 4.5 , 1.0 , 1.0 ) )
    # Initiate plotting ========================================================
    graphics::plot( x = x_axis ,
                    y = TS ,
                    type = 'l' ,
                    xlab = x_lab ,
                    ylab = expression( Target~strength~(dB~re.~1~m^2) ) ,
                    lwd = 2.5 ,
                    col = plot_colors[ 1 : base::length( model_names ) ] ,
                    cex.lab = 1.3 ,
                    cex.axis = 1.15 ,
                    xlim = base::c( base::min( x_axis ) * ( 1 - nudge_x ),
                                    base::max( x_axis ) * ( nudge_x ) ) ,
                    ylim = base::c( base::min( TS ) * ( 1 - ( 1 - nudge_y ) ) ,
                                    base::max( TS ) * ( 1 + ( 1 - nudge_y ) ) ) ,
                    yaxs = "i" ,
                    xaxs = "i" )
    graphics::legend( "topright" ,
                      title = base::expression( bold(TS~model) ) ,
                      title.adj = 0.05 ,
                      legend = model_names ,
                      lty = base::rep( 1 , base::length( model_names ) ) ,
                      col = plot_colors[ 1 : base::length( model_names ) ] ,
                      cex = 1.05 )
  }
  base::invisible( )
}
#' Plotting for SBF-class objects
#' @param object SBF-class object.
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @export
sbf_plot <- function(object,
                     type = "shape",
                     nudge_y = 1.05,
                     nudge_x = 1.01,
                     x_units = "frequency") {
  # Retrieve default plot window parameters ===================
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  if(type == "shape") {
    # Extract body shape information ============================
    body <- extract(object, "body")$rpos
    # Extract bladder information ===============================
    bladder <- extract(object, "bladder")$rpos
    # Set up plot limits for each shape component ===============
    ylow_lim <- min(body[4, ]) * (1 - (1 - nudge_y))
    yhi_lim <- max(body[3, ]) * nudge_y
    xlow_lim <- min(body[1, ]) * (1 - (1 - nudge_x))
    xhi_lim <- max(body[1, ]) * nudge_x
    # Plot dorsoventral view =====================================
    par(ask = F,
        mfrow = c(2, 1),
        mar = c(2, 4.5, 1, 2),
        oma = c(1.5, 0.5, 0.5, 0))
    # Dorsal =====================================================
    plot(x = body[1, ],
         y = body[3, ],
         type = 'l',
         ylim = c(ylow_lim, yhi_lim),
         xlim = c(xlow_lim, xhi_lim),
         xlab = "",
         ylab = "Height (m)",
         lwd = 4,
         xaxt = "n",
         cex.lab = 1.2,
         cex.axis = 1.2)
    # Ventral ====================================================
    lines(x = body[1, ],
          y = body[4, ],
          lty = 1,
          lwd = 4)
    # Connect anterior/posterior gaps =============================
    segments(x0 = c(min(body[1, ]), max(body[1, ])),
             x1 = c(min(body[1, ]), max(body[1, ])),
             y0 = c(body[4, ][body[1, ] == min(body[1, ])],
                    body[4, ][body[1, ] == max(body[1, ])]),
             y1 = c(body[3, ][body[1, ] == min(body[1, ])],
                    body[3, ][body[1, ] == max(body[1, ])]),
             lty = 1,
             lwd = 4,
             col = 'black')
    # Plot body segments  =========================================
    segments(x0 = body[1, ], x1 = body[1, ], y0 = body[4, ], y1 = body[3, ],
             lty = 3,
             lwd = 1,
             col = 'gray30')
    # Ventral bladder ==============================================
    lines(x = bladder[1, ],
          y = bladder[4, ],
          lty = 1,
          lwd = 3.5,
          col = 'red')
    # Dorsal bladder ================================================
    lines(x = bladder[1, ],
          y = bladder[3, ],
          lty = 1,
          lwd = 3.5,
          col = 'red')
    # Connect anterior/posterior gaps =============================
    segments(x0 = c(min(bladder[1, ]), max(bladder[1, ])),
             x1 = c(min(bladder[1, ]), max(bladder[1, ])),
             y0 = c(bladder[4, ][bladder[1, ] == min(bladder[1, ])],
                    bladder[4, ][bladder[1, ] == max(bladder[1, ])]),
             y1 = c(bladder[3, ][bladder[1, ] == min(bladder[1, ])],
                    bladder[3, ][bladder[1, ] == max(bladder[1, ])]),
             lty = 1,
             lwd = 3.5,
             col = 'red')
    # Plot bladder segments ==========================================
    segments(x0 = bladder[1, ], x1 = bladder[1, ],
             y0 = bladder[4, ], y1 = bladder[3, ],
             lty = 3,
             lwd = 1,
             col = 'red')
    # Create legend ==================================================
    legend(x = "bottomright",
           legend = c("Body","Resonant feature"),
           lty = c(1, 1),
           lwd = c(4, 3.5),
           col = c('black', 'red'),
           cex = 0.95)
    # End dorsoventral view =====================================
    left_limit <- -max(body[2, ] / 2) * (1 - (1 - nudge_y))
    right_limit <- max(body[2, ] / 2) * (1 - (1 - nudge_y))
    # Plot ventral view =========================================
    plot(x = body[1, ],
         y = body[2, ] / 2,
         type = 'l',
         ylim = c(left_limit, right_limit),
         ylab = "Width (m)",
         lwd = 4,
         cex.lab = 1.2,
         cex.axis = 1.2)
    # Opposite side ====================================================
    lines(x = body[1, ],
          y = -body[2, ] / 2,
          lty = 1,
          lwd = 4)
    # Close ends ====================================================
    segments(x0 = c(min(body[1, ]), max(body[1, ])),
             x1 = c(min(body[1, ]), max(body[1, ])),
             y0 = c(-body[2, ][body[1, ] == min(body[1, ])] / 2,
                    -body[2, ][body[1, ] == max(body[1, ])] / 2),
             y1 = c(body[2, ][body[1, ] == min(body[1, ])] / 2,
                    body[2, ][body[1, ] == max(body[1, ])] / 2),
             lty = 1,
             lwd = 4,
             col ='black')
    # Plot body segments ==========================================
    segments(x0 = body[1, ], x1 = body[1, ],
             y0 = -body[2, ] / 2, y1 = body[2, ] / 2,
             lty = 3,
             lwd = 1,
             col = 'gray30')
    # Left bladder ==================================================
    lines(x = bladder[1, ],
          y = bladder[2, ] / 2,
          lty = 1,
          lwd = 3.5,
          col = 'red')
    # Right bladder =================================================
    lines(x = bladder[1, ],
          y = -bladder[2, ]/2,
          lty = 1,
          lwd = 3.5,
          col = 'red')
    # Close end of bladders =========================================
    segments(x0 = c(min(bladder[1, ]), max(bladder[1, ])),
             x1 = c(min(bladder[1, ]), max(bladder[1, ])),
             y0 = c(bladder[2, ][bladder[1, ] == min(bladder[1, ])]/ 2,
                    bladder[2, ][bladder[1, ] == max(bladder[1, ])] / 2),
             y1 = c(-bladder[2, ][bladder[1, ] == min(bladder[1, ])] / 2,
                    -bladder[2, ][bladder[1, ] == max(bladder[1, ])] / 2),
             lty = 1,
             lwd = 3.5,
             col = 'red')
    # Plot bladder segments ==========================================
    segments(x0 = bladder[1 ,], x1 = bladder[1, ],
             y0 = -bladder[2, ] / 2, y1 = bladder[2, ] / 2,
             lty = 3,
             lwd = 1,
             col ='red')
    # Re-add x-axis text ==========================================
    mtext('Along-body axis (m)',
          side = 1,
          outer = TRUE,
          line = 0,
          cex = 1.2)
  } else if(type == "model") {
    if(length(extract(object, "model")) == 0) {
      stop("ERROR: no model results detected in object.")
    } else {
      # Extract body shape information ============================
      shape <- extract(object, "body")
      # Extract model results ====================================
      TS <- extract(object, "model")$KRM$TS
      x_axis_domain <- extract(object, "model_parameters")$KRM$parameters
      if(x_units == "frequency") {
        x_axis <- x_axis_domain$acoustics$frequency * 1e-3
        x_lab <- "Frequency (kHz)"
      } else if(x_units == "k_b") {
        x_axis <- x_axis_domain$acoustics$k_sw * max(shape$rpos[1, ])
        x_lab <- expression(italic(k[b]*L))
      }
      # Plot results ===============================================
      par(ask = F,
          mar = c(4, 4.5, 1, 1))
      plot(x = x_axis,
           y = TS,
           type = 'l',
           xlab = x_lab,
           ylab = expression(Target~strength~(dB~re.~1~m^2)),
           lwd = 2.5,
           cex.lab = 1.2,
           cex.axis = 1.2,
           xlim = c(min(x_axis) * (1 - nudge_x),
                    max(x_axis) * (nudge_x)),
           ylim = c(min(TS) * (1 - (1 - nudge_y)),
                    max(TS) * (1 + (1 - nudge_y))),
           xaxs = "i",
           yaxs = "i")
    }
  }
  invisible()
}
