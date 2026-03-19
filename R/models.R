#' Model registry comprising target strength models currently supported by the
#' acousticTS R-package.
#' @keywords internal
#' @noRd
.get_models <- function() {
  list(
    DWBA = DWBA,
    DWBA_curved = DWBA_curved,
    SDWBA = SDWBA,
    SDWBA_curved = SDWBA_curved,
    calibration = calibration,
    SOEMS = calibration,
    ESSMS = ESSMS,
    SPHMS = SPHMS,
    KRM = KRM,
    HPA = HPA,
    PSMS = PSMS,
    FCMS = FCMS
  )
}