% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_scatterer.R
\name{ess_generate}
\alias{ess_generate}
\title{Generate ESS shape}
\usage{
ess_generate(
  shape = "sphere",
  x_body = NULL,
  y_body = NULL,
  z_body = NULL,
  radius_shell,
  shell_thickness = NULL,
  g_fluid = NULL,
  h_fluid = NULL,
  g_shell,
  h_shell,
  theta_shell = pi/2,
  ID = NULL,
  theta_units = "radians",
  length_units = "m"
)
}
\arguments{
\item{shape}{Optional input argument that dictates shape-type, if desired, for
generalized shapes.}

\item{x_body}{Vector containing x-axis body (m) shape data.}

\item{y_body}{Vector containing y-axis body (m) shape data.}

\item{z_body}{Vector containing z-axis body (m) shape data.}

\item{radius_shell}{Radius of shell (m).}

\item{shell_thickness}{Optional shell thickness (m).}

\item{g_fluid}{Optional density contrast for fluid-like body.}

\item{h_fluid}{Optional sound speed contrast for fluid-like body.}

\item{g_shell}{Density contrast for the shell.}

\item{h_shell}{Sound speed contrast for the shell.}

\item{theta_shell}{Object orientation relative to incident sound wave.}

\item{ID}{Optional metadata entry.}

\item{theta_units}{Units used for orientation. Defaults to "radians".}

\item{length_units}{Units used for position vector. Defaults to "m".}
}
\description{
Generate ESS shape
}
