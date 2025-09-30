#' Sample sardine shape with fully inflated swimbladder
#'
#' A pre-generated SBF scatterer containing all information required for
#' target strength modeling.
#'
#' @format A named list with the following components:
#' \describe{
#'   \item{metadata}{Relevant and identifying metadata (\code{list}).}
#'   \item{model_parameters}{Specified model parameters (\code{list}).}
#'   \item{model}{Model outputs and results (\code{list}).}
#'   \item{body}{A list with:
#'     \itemize{
#'       \item \code{rpos}: Position matrix (x, yw, zU, zL; m).
#'       \item \code{sound_speed}: Flesh sound speed (\eqn{c_{body}}, m/s).
#'       \item \code{density}: Flesh density (\eqn{\rho_{body}}, kg/m\eqn{^3}).
#'       \item \code{theta}: Orientation relative to transducer
#'       (\eqn{\theta_{body}}, radians).
#'     }
#'   }
#'   \item{bladder}{A list with:
#'     \itemize{
#'       \item \code{rpos}: Position matrix (x, yw, zU, zL; m).
#'       \item \code{sound_speed}: Bladder sound speed (\eqn{c_{bladder}}, m/s).
#'       \item \code{density}: Bladder density
#'       (\eqn{\rho_{bladder}}, kg/m\eqn{^3}).
#'       \item \code{theta}: Orientation relative to transducer
#'       (\eqn{\theta_{bladder}}, radians).
#'     }
#'   }
#'   \item{shape_parameters}{A named list with:
#'     \itemize{
#'       \item \code{body}: A list describing the body:
#'         \itemize{
#'           \item \code{length}: Body length (m).
#'           \item \code{ncyl}: Number of discrete cylinders.
#'           \item \code{theta_units}: Units for orientation angle.
#'           \item \code{length_units}: Units for length.
#'         }
#'       \item \code{bladder}: A list describing the swimbladder:
#'         \itemize{
#'           \item \code{length}: Bladder length (m).
#'           \item \code{ncyl}: Number of discrete cylinders.
#'           \item \code{theta_units}: Units for orientation angle.
#'           \item \code{length_units}: Units for length.
#'         }
#'     }
#'   }
#' }
#'
#' @usage data(sardine)
"sardine"

#' Sample code shape with fully inflated swimbladder.
#'
#' @format A pre-generated SBF scatterer containing all information required for
#' target strength modeling.
#' \describe{
#'   \item{metadata}{Relevant and identifying metadata (\code{list}).}
#'   \item{model_parameters}{Container for specified model parameters
#'   (\code{list}).}
#'   \item{model}{Model outputs and results (\code{list}).}
#'   \item{body}{List with:
#'     \itemize{
#'       \item \code{rpos}: Position matrix (x, yw, zU, zL; m).
#'       \item \code{sound_speed}: Flesh sound speed (\eqn{c_{body}}, m/s).
#'       \item \code{density}: Flesh density (\eqn{\rho_{body}}, kg/m\eqn{^3}).
#'       \item \code{theta}: Orientation relative to transducer
#'       (\eqn{\theta_{body}}, radians).
#'     }
#'   }
#'   \item{bladder}{List with:
#'     \itemize{
#'       \item \code{rpos}: Position matrix (x, yw, zU, zL; m).
#'       \item \code{sound_speed}: Flesh sound speed (\eqn{c_{bladder}}, m/s).
#'       \item \code{density}: Flesh density
#'       (\eqn{\rho_{bladder}}, kg/m\eqn{^3}).
#'       \item \code{theta}: Orientation relative to transducer
#'       (\eqn{\theta_{bladder}}, radians).
#'     }
#'   }
#'   \item{shape_parameters}{Named list with:
#'     \itemize{
#'       \item \code{body}: List with length (m), ncyl (int),
#'       theta_units (str), length_units (str).
#'       \item \code{bladder}: List with length (m), ncyl (int),
#'       theta_units (str), length_units (str).
#'     }
#'   }
#' }
#' @usage data(cod)
"cod"

#' Sample krill (Euphausia superba) shape taken from McGehee et al. (1998)
#'
#' A dataset containing a sample krill (Euphausia superba) body shape proposed
#'  by McGehee et al. (1998).
#'
#' @format A pre-generated FLS scatterer containing all information required
#' for target strength modeling.
#' \describe{
#'   \item{metadata}{Relevant and identifying metadata (\code{list}).}
#'   \item{model_parameters}{Container for specified model parameters
#'   (\code{list}).}
#'   \item{model}{Model outputs and results (\code{list}).}
#'   \item{body}{A list with elements:
#'     \itemize{
#'       \item \code{rpos}: Position matrix (x, y, z; m).
#'       \item \code{radius}: Radius of each discrete cylinder (m).
#'       \item \code{g}: Body density contrast relative to medium.
#'       \item \code{h}: Sound speed contrast relative to medium.
#'       \item \code{theta}: Orientation angle (\eqn{\theta_{body}}, radians).
#'     }
#'   }
#'   \item{shape_parameters}{A list with:
#'     \itemize{
#'       \item \code{length}: Body length (m).
#'       \item \code{ncyl}: Number of cylinders.
#'       \item \code{theta_units}: Unit of orientation.
#'       \item \code{length_units}: Unit of length.
#'     }
#'   }
#' }
#' @usage data(krill)
"krill"
