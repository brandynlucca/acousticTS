#' Sample sardine shape with fully inflated swimbladder
#'
#' A pre-generated SBF scatterer containing all information required for
#' target strength modeling. The packaged geometry follows the sardine entry
#' distributed through the NOAA Fisheries KRM reference collection and
#' archived in the \code{echoSMs} resources. The object metadata identifies the
#' target as \emph{Sardinops sagax caerulea} following Conti and Demer (2003).
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
#' @source NOAA Fisheries KRM model reference collection
#'   (\url{https://www.fisheries.noaa.gov/data-tools/krm-model}) and the
#'   archived \code{echoSMs} resource set
#'   (\url{https://github.com/ices-tools-dev/echoSMs}).
#' @references
#' Conti, S.G., and Demer, D.A. (2003). Wide-bandwidth acoustical
#' characterization of anchovy and sardine from reverberation measurements in
#' an echoic tank. \emph{ICES Journal of Marine Science}, 60, 617-624.
#' \doi{10.1016/S1054-3139(03)00056-0}
#'
#' @usage data(sardine)
"sardine"

#' Sample cod shape with fully inflated swimbladder
#'
#' @format A pre-generated SBF scatterer containing all information required for
#' target strength modeling. The packaged object corresponds to the historical
#' Atlantic cod example used in the KRM literature and matches the Cod D case
#' archived in the \code{echoSMs} KRM shape collection.
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
#' @source Historical KRM cod shape collection archived in \code{echoSMs}
#'   (\url{https://github.com/ices-tools-dev/echoSMs}) and associated with the
#'   Clay and Horne Atlantic cod example.
#' @references
#' Clay, C.S., and Horne, J.K. (1994). Acoustic models of fish: The Atlantic
#' cod (\emph{Gadus morhua}). \emph{The Journal of the Acoustical Society of
#' America}, 96, 1661-1668. \doi{10.1121/1.410245}
#'
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

#' Benchmark model outputs from Jech et al. (2015)
#'
#' A packaged list of benchmark target-strength spectra for the canonical
#' sphere, shell-sphere, finite-cylinder, and prolate-spheroid comparison cases
#' used throughout the package validation workflow. The stored values follow
#' the benchmark definitions assembled by Jech et al. (2015) and distributed in
#' the \code{echoSMs} reference resources.
#'
#' @source \code{echoSMs} benchmark resources, especially the target-definition
#'   and Jech benchmark tables archived at
#'   \url{https://github.com/ices-tools-dev/echoSMs}.
#' @references
#' Jech, J.M., Horne, J.K., Chu, D., Demer, D.A., Francis, D.T.I., Gorska, N.,
#' Jones, B., Lavery, A.C., Stanton, T.K., and Reeder, D.B. (2015). Comparisons
#' among ten models of acoustic backscattering used in aquatic ecosystem
#' research. \emph{The Journal of the Acoustical Society of America}, 138,
#' 3742-3764. \doi{10.1121/1.4937607}
#'
#' Jech, M., and Macaulay, G. \code{echoSMs}: Acoustic backscattering models
#' used in aquatic ecosystem research. GitHub repository:
#' \url{https://github.com/ices-tools-dev/echoSMs}
#'
#' @usage data(benchmark_ts)
"benchmark_ts"
