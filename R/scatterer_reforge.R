################################################################################
# FORGE FUNCTIONS FOR MANIPULATING SCATTERER SHAPES
################################################################################
################################################################################
# PRIMARY FORGE GENERATION FUNCTION
################################################################################
#' Resize or reparameterize a scatterer object
#'
#' Generic function to resize or reparameterize a scatterer object.
#'
#' @param object A scatterer object.
#' @param ... Additional arguments passed to specific methods.
#' @export
setGeneric(
  "reforge",
  function(object, ...) {
    standardGeneric("reforge")
  }
)

#' Resizing function for swimbladdered targets
#' @param object SBF-class object.
#' @param body_scale Proportional scaling to the body length, width, and height
#' dimensions. When a single value is supplied, all dimensions are scaled using
#' the same scaling factor. Otherwise, this input must be a named numeric
#' vector.
#' @param body_target Target dimensions (m) for the body length, width, and
#' height dimensions. This input must be a named numeric vector.
#' @param swimbladder_scale Proportional scaling to the swimbladder length,
#' width, and height dimensions. When a single value is supplied, all
#' dimensions are scaled using the same scaling factor. Otherwise, this input
#' must be a named numeric vector.
#' @param swimbladder_target Target dimensions (m) for the swimbladder length,
#' width, and height dimensions. This input must be a named numeric vector.
#' @param swimbladder_inflation_factor Proportional swimbladder volume where
#' the swimbladder x-axis origin and terminus are both held constant.
#' @param maintain_ratio Maintain size ratio between body and
#' swimbladder.
#' @param isometric_body Logical; maintain isometric scaling for body.
#' @param isometric_swimbladder Logical; maintain isometric scaling for bladder.
#' @param n_segments_body Number of segments along the body.
#' @param n_segments_swimbladder Number of segments along the bladder.
#' @export
setMethod(
  "reforge",
  signature(object = "SBF"),
  function(object,
           body_scale = NULL,
           body_target = NULL,
           swimbladder_scale = NULL,
           swimbladder_target = NULL,
           isometric_body = TRUE,
           isometric_swimbladder = TRUE,
           maintain_ratio = TRUE,
           swimbladder_inflation_factor = 1.0,
           n_segments_body = NULL,
           n_segments_swimbladder = NULL) {
    ############################################################################
    # Validation ===============================================================
    if (is.null(body_scale) && is.null(swimbladder_scale) && 
        is.null(body_target) && is.null(swimbladder_target) && 
        is.null(n_segments_body) && is.null(n_segments_swimbladder) &&
      swimbladder_inflation_factor == 1.0) {
      stop(
        "Must specify at least one scaling, target, inflation factor, ",
        "or segment count parameter."
      )
    }
    if (!is.null(body_scale) && !is.null(body_target)) {
      stop("Specify only one of body_scale or body_target, not both.")
    }
    if (!is.null(swimbladder_scale) && !is.null(swimbladder_target)) {
      stop(
        "Specify only one of swimbladder_scale or swimbladder_target, not both."
      )
    }
    if (((!is.null(body_scale) && length(body_scale) > 1) || 
         (!is.null(body_target) && length(body_target) > 1)) &&
      isometric_body) {
      message(
        "Hidden State Change Warning:\n ",
        "Multiple axes specified in body_scale/body_target: ",
        "'isometric_body' will be ignored for those axes."
      )
    }
    if (((!is.null(swimbladder_scale) && length(swimbladder_scale) > 1) || 
         (!is.null(swimbladder_target) && length(swimbladder_target) > 1)) &&
      isometric_swimbladder) {
      message(
        "Hidden State Change Warning:\n ",
        "Multiple axes specified in swimbladder_scale/swimbladder_target: ",
        "'isometric_swimbladder' will be ignored for those axes."
      )
    }
    if ((!is.null(body_scale) || !is.null(body_target)) && 
        (!is.null(swimbladder_scale) || !is.null(swimbladder_target)) &&
      maintain_ratio) {
      maintain_ratio <- FALSE
      message(
        "Hidden State Change Warning:\n ",
        "Multiple axes specified for the body and swimbladder: ",
        "'maintain_ratio' will be ignored for those axes."
      )
    }
    ############################################################################
    # Helpers ==================================================================
    get_scale_vector <- function(scale = NULL, target = NULL, dims) {
      # Return nothing if target is not used
      if (is.null(target)) {
        return(scale)
      }

      # Override with target values if given
      out <- c()
      for (nm in names(target)) {
        out[nm] <- target[nm] / dims[nm]
      }
      out
    }

    validate_target <- function(target, name) {
      if (is.null(target)) {
        return(NULL)
      }
      if (!is.numeric(target) || 
          is.null(names(target)) || 
          any(names(target) == "")) {
        stop(paste(name, "must be a named numeric vector."))
      }
      valid_names <- c("length", "width", "height")
      if (!all(names(target) %in% valid_names)) {
        stop(paste("Invalid", name, "names. Use: length, width, height."))
      }
      target
    }

    validate_scale <- function(scale, name, isometry, iso_name) {
      valid_names <- c("length", "width", "height")
      if (is.null(scale)) {
        return(NULL)
      }
      if (is.numeric(scale) && length(scale) == 1 && isometry) {
        if (!all(names(scale) %in% valid_names)) {
          stop(paste("Invalid", name, "names. Use: length, width, height."))
        }
        return(c(
          length = unname(scale),
          width = unname(scale),
          height = unname(scale)
        ))
      }
      if (is.numeric(scale) &&
          length(scale) == 1 && !isometry &&
          is.null(names(scale))) {
        stop(
          paste(
            name,
            "must be a named vector with length/width/height when",
            iso_name,
            "is FALSE."
          )
        )
      }
      if (is.numeric(scale) && !is.null(names(scale))) {
        if (!all(names(scale) %in% valid_names)) {
          stop(paste("Invalid", name, "names. Use: length, width, height."))
        }
        result <- c(length = 1, width = 1, height = 1)
        result[names(scale)] <- scale
        return(result)
      }
      stop(
        paste(
          name,
          "must be a single number or named vector with length/width/height"
        )
      )
    }

    apply_scaling <- function(rpos, scales, dims) {
      if (is.null(scales)) {
        return(rpos)
      }

      # Scale length from anchor point (preserve relative position)
      if (scales["length"] != 1) {
        anchor <- rpos[1, 1]
        rpos[1, ] <- anchor + (rpos[1, ] - anchor) * scales["length"]
      }

      # Scale width and height directly
      if (nrow(rpos) >= 2 && scales["width"] != 1) {
        rpos[2, ] <- rpos[2, ] * scales["width"]
      }
      if (nrow(rpos) >= 3 && scales["height"] != 1) {
        rpos[3, ] <- rpos[3, ] * scales["height"]
        if (nrow(rpos) >= 4) rpos[4, ] <- rpos[4, ] * scales["height"]
      }

      rpos
    }

    interpolate_segments <- function(rpos, n_segments) {
      if (is.null(n_segments)) {
        return(rpos)
      }

      x_new <- seq(rpos[1, 1], rpos[1, ncol(rpos)], length.out = n_segments)
      rpos_new <- rbind(
        x_new,
        t(vapply(2:nrow(rpos), function(i) {
          stats::approx(x = rpos[1, ], y = rpos[i, ], xout = x_new)$y
        }, FUN.VALUE = numeric(length(x_new))))
      )
      rpos_new
    }

    check_swimbladder_containment <- function(rpos_b, rpos_sb) {
      # rpos_b: body rpos matrix
      # rpos_sb: swimbladder rpos matrix

      # Define a common x grid over the overlap of body and bladder
      x_grid <- seq(
        max(min(rpos_b[1, ]), min(rpos_sb[1, ])),
        min(max(rpos_b[1, ]), max(rpos_sb[1, ])),
        length.out = 200
      )

      interp <- function(x, y) stats::approx(x, y, xout = x_grid, rule = 2)$y

      body_y <- interp(rpos_b[1, ], rpos_b[2, ])
      body_zU <- interp(rpos_b[1, ], rpos_b[3, ])
      body_zL <- interp(rpos_b[1, ], rpos_b[4, ])

      bladder_y <- interp(rpos_sb[1, ], rpos_sb[2, ])
      bladder_zU <- interp(rpos_sb[1, ], rpos_sb[3, ])
      bladder_zL <- interp(rpos_sb[1, ], rpos_sb[4, ])

      # Check: width (y) and height (z) containment at each x
      contained <- (bladder_y <= body_y) & (bladder_y >= -body_y) &
        (bladder_zU <= body_zU) & (bladder_zL >= body_zL)

      if (!all(contained)) {
        warning("Swimbladder exceeds body bounds at some positions.")
      }
    }
    ############################################################################
    # Extract components =======================================================
    body <- acousticTS::extract(object, "body")
    bladder <- acousticTS::extract(object, "bladder")
    shape <- acousticTS::extract(object, "shape_parameters")
    rpos_b <- body$rpos
    rpos_sb <- bladder$rpos
    ############################################################################
    # Calculate current dimensions =============================================
    body_height <- max(rpos_b[3, ] - rpos_b[4, ])
    body_dims <- c(
      length = shape$body$length,
      width = max(rpos_b[2, ]),
      height = body_height
    )
    bladder_height <- max(rpos_sb[3, ] - rpos_sb[4, ])
    bladder_dims <- c(
      length = ifelse(is.null(shape$bladder$length),
        max(rpos_sb[1, ]) - min(rpos_sb[1, ]),
        shape$bladder$length
      ),
      width = max(rpos_sb[2, ]),
      height = bladder_height
    )
    ############################################################################
    # Calculate swimbladder origin relative to body position ===================
    original_body_length <- shape$body$length
    bladder_relative_start <- bladder$rpos[1, 1] / original_body_length
    ############################################################################
    # Process target parameters ================================================
    body_target <- validate_target(body_target, "body_target")
    body_scale <- get_scale_vector(body_scale, body_target, body_dims)
    swimbladder_target <- validate_target(
      swimbladder_target,
      "swimbladder_target"
    )
    swimbladder_scale <- get_scale_vector(
      swimbladder_scale,
      swimbladder_target,
      bladder_dims
    )
    ############################################################################
    # Process scaling parameters ===============================================
    body_scales <- validate_scale(
      body_scale,
      "body_scale",
      isometric_body,
      "isometric_body"
    )
    bladder_scales <- validate_scale(
      swimbladder_scale,
      "swimbladder_scale",
      isometric_swimbladder,
      "isometric_swimbladder"
    )
    ############################################################################
    # Apply ratio maintenance logic ============================================
    if (maintain_ratio) {
      if (!is.null(body_scales) && is.null(bladder_scales)) {
        bladder_scales <- body_scales
      } else if (is.null(body_scales) && !is.null(bladder_scales)) {
        body_scales <- bladder_scales
      }
    }
    ############################################################################
    # Interpolate segments first (before scaling)  =============================
    rpos_b <- interpolate_segments(rpos_b, n_segments_body)
    rpos_sb <- interpolate_segments(rpos_sb, n_segments_swimbladder)
    ############################################################################
    # Apply scaling ============================================================
    rpos_b <- apply_scaling(rpos_b, body_scales, body_dims)
    rpos_sb <- apply_scaling(rpos_sb, bladder_scales, bladder_dims)
    ############################################################################
    # Adjust swimbladder position within scaled body if needed =================
    if (!is.null(body_scales) && body_scales["length"] != 1) {
      new_body_length <- max(rpos_b[1, ])
      new_bladder_start <- bladder_relative_start * new_body_length

      # Shift entire bladder to new relative position
      current_bladder_start <- rpos_sb[1, 1]
      shift_amount <- new_bladder_start - current_bladder_start
      rpos_sb[1, ] <- rpos_sb[1, ] + shift_amount
    }
    ############################################################################
    # Apply bladder inflation factor ===========================================
    if (swimbladder_inflation_factor != 1.0) {
      # Preserve relative position to body
      x_bladder_origin <- bladder$rpos[1, 1] / max(body$rpos[1, ])
      xsb_start <- x_bladder_origin * max(rpos_b[1, ])
      xsb_offset <- rpos_sb[1, 1] - xsb_start

      rpos_sb[1, ] <- rpos_sb[1, ] - xsb_offset
      rpos_sb[2, ] <- rpos_sb[2, ] * swimbladder_inflation_factor
      rpos_sb[3, ] <- rpos_sb[3, ] * swimbladder_inflation_factor
      if (nrow(rpos_sb) >= 4) {
        rpos_sb[4, ] <- rpos_sb[4, ] *
          swimbladder_inflation_factor
      }
    }
    ############################################################################
    # Validate swimbladder containment =========================================
    check_swimbladder_containment(rpos_b, rpos_sb)
    ############################################################################
    # Update object ============================================================
    methods::slot(object, "body")$rpos <- rpos_b
    methods::slot(object, "bladder")$rpos <- rpos_sb
    methods::slot(object, "shape_parameters")$body$length <- max(rpos_b[1, ]) -
      min(rpos_b[1, ])
    methods::slot(object, "shape_parameters")$bladder$length <-
      max(rpos_sb[1, ]) - min(rpos_sb[1, ])
    methods::slot(object, "shape_parameters")$body$n_segments <- ncol(rpos_b)
    methods::slot(object, "shape_parameters")$bladder$n_segments <-
      ncol(rpos_sb)
    return(object)
  }
)
################################################################################
#' Reforge FLS-class object.
#' @param object FLS-class object.
#' @param length  New body length resize.
#' @param radius New radius size
#' @param n_segments New number of segments
#' @param length_radius_ratio_constant Keep length-to-radius ratio based on new
#' length
#' @export
setMethod(
  "reforge",
  signature(object = "FLS"),
  function(object,
           length,
           radius,
           length_radius_ratio_constant = TRUE,
           n_segments) {
    ###################################################################
    # Parse shape =====================================================
    shape <- acousticTS::extract(object, "shape_parameters")
    # Parse body ======================================================
    body <- acousticTS::extract(object, "body")
    # Determine rescaling factors =====================================
    # Determine new number of cylinders +++++++++++++++++++++++++++++++
    if (!missing(n_segments)) {
      x_new <- seq(
        from = body$rpos[1, 1],
        to = body$rpos[1, shape$n_segments + 1],
        length.out = n_segments + 1
      )
      rpos_new <- rbind(
        x_new,
        t(
          vapply(2:nrow(body$rpos),
            FUN = function(i) {
              stats::approx(
                x = body$rpos[1, ],
                y = body$rpos[i, ],
                xout = x_new
              )$y
            },
            FUN.VALUE = numeric(base::length(x_new))
          )
        )
      )
      radius_new <- stats::approx(
        x = body$rpos[1, ],
        y = body$radius,
        xout = x_new
      )$y
      # Update metadata +++++++++++++++++++++++++++++++++++++++++++++++
      methods::slot(object, "body")$rpos <- rpos_new
      methods::slot(object, "body")$radius <- radius_new
      methods::slot(object, "shape_parameters")$n_segments <- n_segments
    }
    # # Update material properties if needed ============================
    # if ( length( body$g ) > 1 ) {
    #  test = approx( x = body$rpos[ 1, ] ,
    #          y = c( body$g, body$g[ length( body$g ) ] ) ,
    #                 xout = x_new )$y
    # }
    # Determine new length ++++++++++++++++++++++++++++++++++++++++++++
    if (!missing(length)) {
      new_scale <- length / shape$length
      matrix_rescale <- diag(
        x = 1,
        nrow = nrow(body$rpos),
        ncol = nrow(body$rpos)
      ) * new_scale
      rpos_new <- t(t(body$rpos) %*% matrix_rescale)
      # New radius based on constant ratio or adjust ++++++++++++++++++
      if (length_radius_ratio_constant) {
        if (missing(radius)) {
          radius_new <- body$radius * new_scale
        } else {
          radius_rescale <- radius / shape$radius
          radius_new <- body$radius * radius_rescale
        }
        methods::slot(object, "body")$radius <- radius_new
        methods::slot(object, "shape_parameters")$radius <- max(radius_new)
      }
      # Update metadata +++++++++++++++++++++++++++++++++++++++++++++++
      methods::slot(object, "body")$rpos <- rpos_new
      methods::slot(object, "shape_parameters")$length <- max(rpos_new[1, ])
    }
    # Determine new radius ++++++++++++++++++++++++++++++++++++++++++++
    if (!missing(radius)) {
      new_scale <- radius / shape$radius
      vector_rescale <- body$radius * new_scale
      methods::slot(object, "body")$radius <- vector_rescale
      methods::slot(object, "shape_parameters")$radius <- max(vector_rescale)
    }
    # Return object ===================================================
    return(object)
  }
)

#' Get reforge parameters from known method signatures
#' @param object_class Character string of object class
#' @return Character vector of parameter names
#' @keywords internal
#' @noRd
.discover_reforge_params <- function(object_class) {
  switch(object_class,
    "FLS" = c(
      "length", "radius", "length_radius_ratio_constant",
      "n_segments"
    ),
    "SBF" = c(
      "body_scale", "swimbladder_scale", "body_target",
      "swimbladder_target", "swimbladder_inflation_factor",
      "isometric_body", "isometric_swimbladder", "maintain_ratio",
      "n_segments_body", "n_segments_swimbladder"
    ),
    character(0) # Default for unknown classes
  )
}
