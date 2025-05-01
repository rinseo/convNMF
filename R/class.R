#-------------------#
# Define S4 Classes #
#-------------------#
#' Wavelets Class
#'
#' Stores spatial wavelets.
#'
#' @slot value A 3D array [neurons × channels × delay] of spatial wavelets.
#' @slot num_neurons Number of neurons.
#' @slot num_channels Number of channels (electrodes/tasks).
#' @slot len_wavelet Length of each spatial wavelet (in samples).
#' @slot weights Weights after re-normalization.
#' @slot channel_order Order of channels.
#' @slot channel_labels List of channel labels.
#'
#' @seealso \code{\link{plot.Wavelets}}
#'
#' @docType class
#' @name Wavelets-class
#' @exportClass Wavelets
setClass("Wavelets",
         slots=c(value="array",
                 num_neurons="numeric",
                 num_channels="numeric",
                 len_wavelet="numeric",
                 weights="array",
                 channel_order="numeric",
                 channel_labels="character"))
#' Amplitudes Class
#'
#' Stores temporal amplitudes.
#'
#' @slot value A 2D array [neurons × time] of temporal amplitudes.
#' @slot num_neurons Number of neurons.
#' @slot num_samples Number of time points.
#'
#' @seealso \code{\link{plot.Amplitudes}}
#'
#' @docType class
#' @name Amplitudes-class
#' @exportClass Amplitudes
setClass("Amplitudes",
         slots=c(value="array",
                 num_neurons="numeric",
                 num_samples="numeric"))
#' Data Class
#'
#' Stores multi-channel time-series data.
#'
#' @slot value A 2D array [channels × time] of input data.
#' @slot num_channels Number of channels (electrodes/tasks).
#' @slot num_samples Number of time points.
#' @slot channel_order Order of channels.
#' @slot channel_labels List of channel labels.
#'
#' @seealso \code{\link{plot.Data}}
#'
#' @docType class
#' @name Data-class
#' @exportClass Data
setClass("Data",
         slots=c(value="array",
                 num_channels="numeric",
                 num_samples="numeric",
                 channel_order="numeric",
                 channel_labels="character"))
#' ConvNMF Class
#'
#' Stores multi-channel time-series data, spatial wavelets, and temporal amplitudes.
#'
#' @slot wavelets An object of \code{\link{Wavelets-class}}.
#' @slot amplitudes An object of \code{\link{Amplitudes-class}}.
#' @slot data An object of \code{\link{Data-class}}.
#' @slot lls An numeric vector of normalized log likelihood values over iterations for fitted objects.
#' @slot fit_args Named arguments passed to \code{\link{convNMF}}.
#'
#' @seealso \code{\link{plot,ConvNMF-method}}, 
#' \code{\link{compute_diss_matrix}},
#' \code{\link{confusion_matrix}},
#' #' \code{\link{dissimilarity}},
#' \code{\link{extract_spikes}},
#' \code{\link{logLik,ConvNMF-method}},
#' \code{\link{power_explained}},
#' \code{\link{predict,ConvNMF-method}},
#' \code{\link{sequenciness}},
#' \code{\link{update,ConvNMF-method}}
#'
#' @docType class
#' @name ConvNMF-class
#' @exportClass ConvNMF
setClass("ConvNMF",
         slots=c(wavelets="Wavelets",
                 amplitudes="Amplitudes",
                 data="Data",
                 lls="numeric",
                 fit_args="list"))
