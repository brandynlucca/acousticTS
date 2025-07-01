################################################################################
################################################################################
# GENERIC PLOT FUNCTION FOR BODY SHAPES AND MODELS
################################################################################
################################################################################
#' Method for what is printed for objects.
#' @param x Scatterer-class object.
#' @param y Ignored (required for plot method signature).
#' @param ... Additional arguments passed to plotting functions
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param aspect_ratio Aspect ratio setting ( defaults to "manual" for nudge_y
#' and nudge_x to apply; otherwise, input "auto").
#' @param y_units y-axis data selection (e.g. TS, sigma_bs -- defaults to TS)
#' @param ... Additional plot inputs
#' @rdname plot.scatterer
#' @export
setGeneric( "plot" , function( x , y , ... ) standardGeneric( "plot" ) )
#' @rdname plot.scatterer
#' @export
setMethod( f = "plot" ,
          signature = c( x = "scatterer" , y = "missing" ) ,
          definition = function( x ,
                                 y ,
                                 type = "shape" ,
                                 nudge_y = 1.1 ,
                                 nudge_x = 1.05 ,
                                 aspect_ratio = "manual" ,
                                 x_units = "frequency",
                                 y_units = "TS" , ... ) {
            # Detect scatterer type ============================================
            sc_type <- class( x )

            switch( sc_type,
                    CAL = cal_plot( x , type , nudge_y , nudge_x , x_units , ... ) ,
                    SBF = sbf_plot( x , type , nudge_y , nudge_x , x_units , ... ) ,
                    FLS = fls_plot( x , type , nudge_y , nudge_x , aspect_ratio , x_units , y_units ) ,
                    GAS = gas_plot( x , type , nudge_y , nudge_x , x_units , y_units , ... ) )
          } )
################################################################################
#' Base plotting color palette
#' @description 
#' Color palette vector referenced for plotting multiple models simultaneously.
#' @rdname model_palette
#' @export
model_palette <- c( "black" , "azure4" , "royalblue4" , "firebrick3" , "steelblue3" ,
                    "orangered3" , "darkolivegreen" , "sienna3" , "cadetblue" )
################################################################################
################################################################################
# ELASTIC SHELLED SCATTERERS
################################################################################
################################################################################
#' Plotting for CAL-class objects
#' @param object CAL-class object.
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param ... Additional plot inputs
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
################################################################################
# FLUID-LIKE SCATTERERS
################################################################################
################################################################################
#' Plotting for FLS-class objects
#' @param object FLS-class object.
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param aspect_ratio Aspect ratio setting ( defaults to "manual" for nudge_y
#' and nudge_x to apply; otherwise, input "auto").
#' @param y_units y-axis data selection (e.g. TS, sigma_bs -- defaults to TS)
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
    body$rpos[ 3 , ] <- body$rpos[ 3 , ] - median( body$rpos[ 3 , ] )
    # Sort index ===============================================================
    if ( body$rpos[ 1 , 1 ] > body$rpos[ 1, ncol( body$rpos ) ] ) {
      body$rpos <- body$rpos[ , rev( seq_len( ncol( body$rpos ) ) ) ,
                             drop = FALSE ]
      if ( !is.null( body$radius ) ) body$radius <- rev( body$radius )
      # sort_idx <- order(body$rpos[1, ])
      # body$rpos <- body$rpos[, sort_idx, drop = FALSE]
      # if (!is.null(body$radius)) body$radius <- body$radius[sort_idx]
    }
    # Adjust axes ==============================================================
    if ( aspect_ratio == "manual" ) {
      vert_lims <- base::c( min( body$rpos[ 3 , ] - body$radius ) * ( 1 - ( 1 - nudge_y ) ) ,
                            max( body$rpos[ 3 , ] + body$radius ) * nudge_y )
    } else {
      vert_lims <- base::c( - base::max( body$rpos[ 1 , ]  )  * 0.10 ,
                            base::max( body$rpos[ 1 ,  ]  )  * 0.10 )
    }
    # Begin plotting ===========================================================
    # graphics::plot( x = body$rpos[ 1 , ] ,
    #                 y = body$rpos[ 3 , ] ,
    #                 type = 'l' ,
    #                 lwd = 4 ,
    #                 cex.lab = 1.2 ,
    #                 cex.axis = 1.2 ,
    #                 xlab = "Length (m)" ,
    #                 ylab = "Thickness (m)" ,
    #                 ylim = vert_lims )
    plot( body$rpos[ 1 , ] , body$rpos[ 3 , ] , type = 'n',
          xlab = "Length (m)", ylab = "Thickness (m)",
          ylim = c( min(body$rpos[ 3 , ] - body$radius ) * 1.1 ,
                    max(body$rpos[ 3 , ] + body$radius ) * 1.1 ) )
    # Add lower perimeter of shape =============================================
    # graphics::lines( x = body$rpos[ 1 , ] ,
    #                  y = body$rpos[ 3 , ] - body$radius ,
    #                  lty = 1 ,
    #                  lwd = 4 )
    # Add upper perimeter of shape =============================================
    # graphics::lines( x = body$rpos[ 1 , ] ,
    #                  y = body$rpos[ 3 , ] + body$radius ,
    #                  lty = 1 ,
    #                  lwd = 4 )
    # Add body segments ========================================================
    # graphics::segments( x0 = body$rpos[ 1 , ] ,
    #                     x1 = body$rpos[ 1 , ] ,
    #                     y0 = body$rpos[ 3 , ] - body$radius ,
    #                     y1 = body$rpos[ 3 , ] + body$radius ,
    #                     lty = 3 ,
    #                     lwd = 1.25 )
    # Draw angled "cylinders" for each segment =================================
    n_segments <- acousticTS::extract( object , 
                                       "shape_parameters" )$n_segments - 1
    # count <- 0
    # rc <- c()
    # ---- Iterate =============================================================
    for ( i in 1 : ( n_segments ) ) {
      # ---- Leading coordinates
      x0 <- body$rpos[ 1 , i ]
      y0 <- body$rpos[ 3 , i ]
      # ---- Trailing coordinates
      x1 <- body$rpos[1 , i + 1 ]
      y1 <- body$rpos[3 , i + 1 ]
      # ---- Thickness
      r <- body$radius[ i ]
      # rc <- c(rc, r)
      if ( r == 0 ) next
      # cat(sprintf("Segment %d [%d, %d]: (%.6f, %.6f) -> (%.6f, %.6f), r=%.6f\n", 
      #             i, order(body$rpos[1,])[i], order(body$rpos[1,])[i+1], x0, y0, x1, y1, r))
      # ---- Direction vector
      dx <- x1 - x0
      dy <- y1 - y0
      len <- sqrt( dx ^2 + dy ^2 )
      # ---- Perpendicular vector (unit)
      px <- - dy / len
      py <- dx / len
      # ---- Polygon nodes
      xA <- x0 + r * px
      yA <- y0 + r * py
      xB <- x1 + r * px
      yB <- y1 + r * py
      xC <- x1 - r * px
      yC <- y1 - r * py
      xD <- x0 - r * px
      yD <- y0 - r * py
      # ---- Draw polygon
      polygon(
        x = c( xA , xB , xC , xD ),
        y = c( yA , yB , yC , yD ),
        col = adjustcolor( "gray50" , alpha.f = 0.6 ) ,
        border = "black", lwd = 1
      )
      # count <- count + 1
    }
    # cat("Polygons drawn:", count, "\n")
    # Draw centerline and points
    lines(body$rpos[1, ], body$rpos[3, ], lwd = 3, col = "gray90")
    points(body$rpos[1, ], body$rpos[3, ], pch = 1, col = "black", cex = 0.8)
  } else if ( type == "model" ) {
    # Detect model selection ===================================================
    models <- extract( object , "model" )
    model_names <- names( models )
    if ( base::length( model_names ) == 0 )
      stop( "ERROR: no model results detected in object." )
    # Extract body shape information ===========================================
    shape <- extract( object , "body" )
    # Append model name ========================================================
    models <- lapply( 1 : length( model_names ) ,
                      FUN = function( x ) {
                        models[[ x ]] <- models[[ x ]][ c( "frequency" , "TS" ) ]
                        transform( models[[ x ]] ,
                                   model = model_names[ x ] ) } )
    # Convert into a data.frame ================================================
    models_df <- do.call( "rbind" , models )
    # Define x-axis domain =====================================================
    # x_axis <- base::switch( x_units ,
    #                         frequency = base::unique( models_df$frequency ) * 1e-3 ,
    #                         k_sw = base::unique( models_df$ka ) )
    x_axis <- models_df[ , base::which( base::colnames( models_df )  == x_units ) ]
    x_mat <- base::split( x_axis , models_df$model )
    y_axis <- models_df[ , base::which( base::colnames( models_df ) == y_units ) ]
    y_mat <- base::split( y_axis , models_df$model )
    col_axis <- model_palette[ base::as.numeric( base::as.factor( models_df$model ) ) ]
    x_lab <- switch( x_units ,
                     frequency = "Frequency (Hz)" ,
                     k_sw = expression(italic(k[sw]*a) ) )
    # Define plot margins ======================================================
    graphics::par( ask = FALSE ,
                   mar = base::c( 5.0 , 5.5 , 2.0 , 3.5 ) )
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
                    xaxt = 'n' ,
                    type = 'n' ,
                    cex.axis = 1.3 ,
                    cex.lab = 1.5 )
    atx <- seq(par("xaxp")[1], par("xaxp")[2], (par("xaxp")[2] - par("xaxp")[1])/par("xaxp")[3])
    axis( 1 , at = atx , 
          labels = format( atx , scientific = F ) ,
          cex.axis = 1.3 )
    nm_names <- names( x_mat )
    invisible( mapply( lines , x_mat , y_mat ,
                       col = model_palette[ 1 : base::length( model_names ) ] ,
                       lwd = 4 ) )
    legend( "bottomright" ,
            title = expression( bold("TS"~"model") ) ,
            title.adj = 0.05 ,
            legend = nm_names ,
            lty = rep( 1 , length( nm_names ) ) ,
            col = model_palette[ 1 : length( nm_names  ) ] ,
            cex = 1.05 ,
            lwd = 4 )
   }
  base::invisible( )
}
################################################################################
################################################################################
# GAS-BEARING SCATTERERS
################################################################################
################################################################################
#' Plotting for GAS-class objects
#' @param object GAS-class object.
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param y_units y-axis data selection (e.g. TS, sigma_bs -- defaults to TS).
#' @param ... Additional plot inputs
#' @import graphics
#' @import stats
#' @import grDevices
#' @export
gas_plot <- function( object ,
                      type = "shape" ,
                      nudge_y = 1.01 ,
                      nudge_x = 1.01 ,
                      x_units = "frequency" ,
                      y_units = "TS" , ... ) {
  # Retrieve default plot window parameters ====================================
  opar <- par( no.readonly = TRUE )
  on.exit( par( opar ) )
  if( type == "shape" ) {
    # Extract body shape information ===========================================
    body <- extract( object , "body" )
    # Define plot margins ======================================================
    par( ask = FALSE ,
         oma = c( 1 , 1 , 1 , 0 ) ,
         mar = c( 3 , 4.5 , 1 , 2 ) )
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
    models <- extract( object , "model" )
    model_names <- names( models )
    if ( length( model_names ) == 0 )
      stop( "ERROR: no model results detected in object." )
    # Extract body shape information ===========================================
    shape <- extract( object , "body" )
    # Append model name ========================================================
    models <- lapply( 1 : length( model_names ) ,
                      FUN = function( x ) {
                        models[[ x ]] <- models[[ x ]][ c( "frequency" , "TS" ) ]
                        transform( models[[ x ]] ,
                                   model = model_names[ x ] ) } )
    # Convert into a data.frame ================================================
    models_df <- do.call( "rbind" , models )
    # Define x-axis domain =====================================================
    # x_axis <- base::switch( x_units ,
    #                         frequency = base::unique( models_df$frequency ) * 1e-3 ,
    #                         k_sw = base::unique( models_df$ka ) )
    x_axis <- models_df[ , base::which( base::colnames( models_df )  == x_units ) ]
    x_mat <- base::split( x_axis , models_df$model )
    y_axis <- models_df[ , base::which( base::colnames( models_df ) == y_units ) ]
    y_mat <- base::split( y_axis , models_df$model )
    col_axis <- model_palette[ base::as.numeric( base::as.factor( models_df$model ) ) ]
    x_lab <- switch( x_units ,
                     frequency = "Frequency (Hz)" ,
                     k_sw = expression(italic(k[sw]*a) ) )
    # Define plot margins ======================================================
    graphics::par( ask = FALSE ,
                   mar = base::c( 5.0 , 5.5 , 2.0 , 3.5 ) )
    # Initiate plotting ========================================================
    plot( x = seq( from = min( x_axis ) ,
                   to = max( x_axis ) ,
                   length.out = 2 ) ,
          y = seq( from = min( y_axis ) ,
                   to = max( y_axis ) ,
                   length.out = 2 ) ,
          xlim = c( min( x_axis ) * ( 1 - nudge_x ),
                    max( x_axis ) * ( nudge_x ) ) ,
          ylim = c( min( y_axis ) * ( 1 - ( 1 - nudge_y ) ) ,
                    max( y_axis ) * ( 1 + ( 1 - nudge_y ) ) ) ,
          xlab = x_lab ,
          ylab = expression( "Target"~"strength"~("dB"~"re."~1~"m"^2) ) ,
          yaxs = "i" ,
          xaxs = "i" ,
          xaxt = 'n' ,
          type = 'n' ,
          cex.axis = 1.3 ,
          cex.lab = 1.5 )
    atx <- seq(par("xaxp")[1], par("xaxp")[2], (par("xaxp")[2] - par("xaxp")[1])/par("xaxp")[3])
    axis( 1 , at = atx , 
          labels = format( atx , scientific = F ) ,
          cex.axis = 1.3 )
    nm_names <- names( x_mat )
    invisible( mapply( lines , x_mat , y_mat ,
                       col = model_palette[ 1 : base::length( model_names ) ] ,
                       lwd = 4 ) )
    legend( "bottomright" ,
            title = expression( bold("TS"~"model") ) ,
            title.adj = 0.05 ,
            legend = nm_names ,
            lty = rep( 1 , length( nm_names ) ) ,
            col = model_palette[ 1 : length( nm_names  ) ] ,
            cex = 1.05 ,
            lwd = 4 )
  }
  invisible( )
}
################################################################################
#' Plotting for SBF-class objects
#' @param object SBF-class object.
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @export
sbf_plot <- function(object,
                     type = "shape" ,
                     nudge_y = 1.05 ,
                     nudge_x = 1.01 ,
                     x_units = "frequency" ) {
  # Retrieve default plot window parameters ===================
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  if(type == "shape") {
    # Extract body shape information ============================
    body <- extract( object , "body" )$rpos
    # Extract bladder information ===============================
    bladder <- extract( object , "bladder" )$rpos
    # Set up plot limits for each shape component ===============
    ylow_lim <- min( body[ 4 , ] ) * ( 1 - ( 1 - nudge_y ) )
    yhi_lim <- max( body[ 3 , ] ) * nudge_y
    xlow_lim <- min( body[ 1 , ] ) * ( 1 - ( 1 - nudge_x ) )
    xhi_lim <- max( body[ 1 , ] ) * nudge_x
    # Plot dorsoventral view =====================================
    par(ask = F,
        mfrow = c(2, 1),
        mar = c(2, 4.5, 1, 2),
        oma = c(1.5, 0.5, 0.5, 0))
    # Dorsal =====================================================
    plot(x = body[ 1 , ] ,
         y = body[ 3 , ] ,
         type = 'l',
         ylim = c( ylow_lim , yhi_lim ) ,
         xlim = c( xlow_lim , xhi_lim ) ,
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
