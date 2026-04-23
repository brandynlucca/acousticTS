source("tools/implementation-figures/helpers/common.R")
impl_load_all()

density_sw <- 1026.8
sound_speed_sw <- 1477.3
theta_body <- pi / 2
phi_body <- pi / 2
theta_scatter <- rep(pi / 2, 181)
phi_scatter <- seq(0, 2 * pi, length.out = 181)
quad <- gauss_legendre(n = 96, a = -1, b = 1)

build_case <- function(boundary, Amn_method) {
  object <- target_strength(
    object = fls_generate(
      shape = prolate_spheroid(length_body = 0.14, radius_body = 0.01, n_segments = 80),
      g_body = 1,
      h_body = 1,
      theta_body = theta_body
    ),
    frequency = 38e3,
    model = "tmm",
    boundary = boundary,
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  acoustics <- object@model_parameters$TMM$parameters$acoustics
  body <- object@model_parameters$TMM$body
  medium <- object@model_parameters$TMM$medium

  exact_f <- complex(length(phi_scatter))

  # Evaluate the exact general-angle spheroidal field at the same angular
  # locations used by the stored TMM post-processing path.
  for (i in seq_along(phi_scatter)) {
    body$theta_body <- theta_body
    body$phi_body <- phi_body
    body$theta_scatter <- theta_scatter[i]
    body$phi_scatter <- phi_scatter[i]

    exact_raw <- prolate_spheroid_fbs(
      acoustics = acoustics,
      body = body,
      medium = medium,
      integration_pts = quad,
      precision = "double",
      Amn_method = Amn_method
    )
    exact_f[i] <- (-2i / acoustics$k_sw) * exact_raw
  }

  stored_sigma_scat_dB <- vapply(
    seq_along(phi_scatter),
    function(i) {
      tmm_scattering(
        object,
        theta_body = theta_body,
        phi_body = phi_body,
        theta_scatter = theta_scatter[i],
        phi_scatter = phi_scatter[i]
      )$sigma_scat_dB
    },
    numeric(1)
  )

  data.frame(
    boundary = boundary,
    phi_scatter = phi_scatter,
    exact_sigma_scat_dB = 10 * log10(Mod(exact_f)^2),
    stored_sigma_scat_dB = stored_sigma_scat_dB
  )
}

df <- rbind(
  build_case("pressure_release", "Amn_pressure_release"),
  build_case("fixed_rigid", "Amn_fixed_rigid")
)

write.csv(
  df,
  "tools/implementation-figures/data/tmm/tmm_prolate_exact_angular_validation.csv",
  row.names = FALSE
)

png(
  filename = "vignettes/tmm/tmm-prolate-exact-angular-validation.png",
  width = 1800,
  height = 800,
  res = 200
)

par(mfrow = c(1, 2), mar = c(4.6, 4.8, 3.2, 1.2), oma = c(0, 0, 0, 0))

plot_case <- function(case_df, title_text) {
  y_range <- range(
    c(case_df$exact_sigma_scat_dB, case_df$stored_sigma_scat_dB),
    finite = TRUE
  )

  plot(
    case_df$phi_scatter,
    case_df$exact_sigma_scat_dB,
    type = "l",
    lwd = 2,
    col = "black",
    xlab = expression(phi[scatter] ~ "(rad)"),
    ylab = expression(10 * log[10](sigma[scat] / (1 ~ m^2)) ~ "(dB)"),
    main = title_text,
    ylim = y_range
  )

  point_index <- seq(1, nrow(case_df), by = 6)
  points(
    case_df$phi_scatter[point_index],
    case_df$stored_sigma_scat_dB[point_index],
    pch = 1,
    cex = 0.8,
    lwd = 1.4,
    col = "#D55E00"
  )

  legend(
    "topright",
    legend = c("Exact spheroidal", "Stored TMM"),
    lty = c(1, NA),
    lwd = c(2, NA),
    pch = c(NA, 1),
    pt.cex = c(NA, 0.8),
    col = c("black", "#D55E00"),
    bty = "n"
  )
}

plot_case(
  subset(df, boundary == "pressure_release"),
  "Pressure-release prolate, 38 kHz"
)
plot_case(
  subset(df, boundary == "fixed_rigid"),
  "Rigid prolate, 38 kHz"
)

dev.off()
