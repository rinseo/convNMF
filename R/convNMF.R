#' @import dplyr
#' @import torch
#' @import progress
#' @importFrom magrittr "%>%"
#' @importFrom pracma findpeaks
#' @importFrom fido gather_array
#' @importFrom reshape2 melt

NULL

# version 2025.07.15
#-----------------#
# Data Generation #
#-----------------#
#' Generate Random Spatial Wavelets
#'
#' Creates a set of random spatial wavelets for simulated neurons across multiple channels.
#'
#' @param num_neurons Number of neurons to simulate.
#' @param num_channels Number of channels (electrodes/tasks).
#' @param len_wavelet Length of each spatial wavelet (in samples).
#' @param warp Function for time warping.
#'
#' @return A \code{\link{Wavelets-class}} object.
#' @export
generate_wavelets <- function(num_neurons, num_channels, len_wavelet,
                              warp=function(x) -exp(-x)+1){
  # Make (semi) random wavelets
  # num_neurons <- 2; num_channels <- 5; len_wavelet <- 10
  new_wavelets <- array(dim=c(num_neurons, num_channels, len_wavelet))
  times <- array(dim=c(num_neurons))
  for(k in 1:num_neurons){
    # Space warping using a Gaussian spatial kernel
    center <- runif(1, min=0, max=num_channels)
    width <- runif(1, min=1, max=1 + num_channels/10)
    spatial_factor <- exp(-0.5*(seq(0,num_channels-1) - center)**2 / width**2)
    
    # Time warping using a Gaussian temporal kernel
    dt <- seq(0, len_wavelet-1)
    period <- len_wavelet / runif(1, min=1, max=2)
    peak_shift <- runif(1, min=0.25, max=0.75) # 0.5 = peak at midpoint
    spread_factor <- runif(1, min=0.25, max=0.75) # smaller = narrower
    z <- (dt - peak_shift * period) / (spread_factor * period)
    # warp <- function(x) exp(-x) - 1 # positive spike
    # warp <- function(x) -exp(-x) + 1 # negative spike
    window <- exp(-0.5 * z**2)
    shape <- sin(2 * pi * dt / period)
    temporal_factor <- warp(window * shape)
    
    # Outer product of spatial x temporal factor
    wavelet <- outer(spatial_factor, temporal_factor)
    wavelet <- wavelet / norm(wavelet,"2") # L2 norm
    new_wavelets[k,,] <- wavelet
    times[k] <- which(abs(wavelet)==max(abs(wavelet)),arr.ind=T)[2] # get col index
  }
  # Reorder by spiking times
  wavelets_r <- new_wavelets[order(times),,,drop=FALSE]
  wavelets <- new("Wavelets",
                  value=wavelets_r,
                  num_neurons=num_neurons, # K
                  num_channels=num_channels, # N 
                  len_wavelet=len_wavelet, # D
                  channel_order=1:num_channels,
                  channel_labels=as.character(1:num_channels)
  )
  return(wavelets)
}

#' Generate Synthetic Neural Data
#'
#' Simulates multi-channel time-series data with spatial wavelets and temporal amplitudes.
#'
#' @param num_channels Number of channels (electrodes/tasks).
#' @param num_samples Number of time points.
#' @param sample_freq Sampling frequency.
#' @param num_neurons Number of neurons to simulate.
#' @param len_wavelet Length of spatial wavelet (in samples).
#' @param warp Function for time warping.
#' @param mean_spikes Expected number of spikes per unit time.
#' @param mean_amplitude Mean of the Gamma distribution for amplitudes.
#' @param shape_amplitude Shape of the Gamma distribution for amplitudes.
#' @param noise_std Standard deviation of Gaussian noise.
#'
#' @return A synthetic \code{\link{ConvNMF-class}} object.
#' @export
generate <- function(num_channels, # number of channels
                     num_samples, # number of samples
                     sample_freq, # number of samples per unit time
                     num_neurons, # number of neurons
                     len_wavelet, # length of spikes
                     warp=function(x) -exp(-x)+1, # warp function
                     mean_spikes=10, # expected number of spikes per unit time
                     mean_amplitude=15, # alpha/beta of Gamma (alpha, beta)
                     shape_amplitude=3, # alpha of Gamma (alpha, beta)
                     noise_std=1 # MVN(0, noise_std^2)
){
  
  # Make random wavelets
  # num_neurons <- 2; num_channels <- 5; len_wavelet <- 10
  wavelets <- generate_wavelets(num_neurons, num_channels, len_wavelet, warp)@value
  
  # Make random amplitudes
  # num_samples <- 1000; sample_freq <- 1000
  # mean_spikes <- 10; mean_amplitude <- 15; shape_amplitude <- 3
  amplitudes <- array(0, dim=c(num_neurons, num_samples))
  for (k in 1:num_neurons){
    num_spikes <- rpois(n=1, lambda=num_samples / sample_freq * mean_spikes)
    times <- sample(1:num_samples, size=num_spikes, replace=T)
    amps <- rgamma(n=num_spikes, shape=shape_amplitude, rate=shape_amplitude / mean_amplitude)
    amplitudes[k, times] <- amps
    
    # Find the peaks in the cross-correlation
    amps_peaks <- findpeaks(amplitudes[k],
                            minpeakheight = 1e-3,
                            minpeakdistance = len_wavelet)
    heights <- amps_peaks[,1]
    peaks <- amps_peaks[,2]
    
    # Only keep spikes separated by at least D
    amplitudes[k] <- 0
    amplitudes[k, peaks] <- heights
  }
  
  # Compute the model prediction by convolving the amplitude and wavelets
  torch_amplitudes <- torch_tensor(amplitudes)
  torch_wavelets <- ensure3D(torch_tensor(wavelets))
  torch_wavelets_flip <- torch_flip(torch_wavelets, dims=-1)
  torch_data <- torch_conv1d(torch_amplitudes,
                             torch_wavelets_flip$permute(c(2,1,3)),
                             padding=len_wavelet-1)[,1:num_samples]
  data <- as.array(torch_data) # convert back to array
  
  # Make random noises
  data[,] <- data[,] + rnorm(length(data), mean=0, sd=noise_std)
  
  return(new("ConvNMF",
             wavelets=new("Wavelets",
                          value=wavelets,
                          num_neurons=num_neurons, # K
                          num_channels=num_channels, # N 
                          len_wavelet=len_wavelet, # D
                          channel_order=1:num_channels,
                          channel_labels=as.character(1:num_channels)
             ),
             amplitudes=new("Amplitudes",
                            value=amplitudes, 
                            num_neurons=num_neurons, # K
                            num_samples=num_samples), # T
             data=new("Data",
                      value=data,
                      num_channels=num_channels, # N 
                      num_samples=num_samples, # T
                      channel_order=1:num_channels,
                      channel_labels=as.character(1:num_channels)),
             fit_args=list(
               num_channels=num_channels,
               num_samples=num_samples,
               sample_freq=sample_freq,
               num_neurons=num_neurons,
               len_wavelet=len_wavelet,
               mean_spikes=mean_spikes,
               mean_amplitude=mean_amplitude,
               shape_amplitude=shape_amplitude,
               noise_std=noise_std
             )))
}

#------------------------#
# ConvNMF Implementation #
#------------------------#
ensure3D <- function(tensor) {
  # Ensure that input tensor is 3D array
  shape <- tensor$shape
  if (length(shape) == 2) {
    return(tensor$unsqueeze(1))  # Add K dimension
  }
  return(tensor)
}
log_likelihood <- function(wavelets, amplitudes, data, noise_std=1.0){
  # Evaluate the (normalized) log likelihood
  K <- dim(wavelets)[1]
  N <- dim(wavelets)[2]
  D <- dim(wavelets)[3]
  Ts <- dim(data)[2]
  
  # Compute the model prediction by convolving the amplitude and wavelets
  torch_amplitudes <- torch_tensor(amplitudes)
  torch_wavelets <- ensure3D(torch_tensor(wavelets))
  torch_wavelets_flip <- torch_flip(torch_wavelets, dims=-1)
  torch_pred <- torch_conv1d(torch_amplitudes,
                             torch_wavelets_flip$permute(c(2,1,3)),
                             padding=D-1)[,1:Ts]
  pred <- as.array(torch_pred) # convert back to array
  
  # Evaluate the log probability using dist.Normal
  ll <- sum(dnorm(data, mean=pred, sd=noise_std, log=T))
  
  # Return the log probability noramlized by the data size
  return(ll / (N * Ts))
}

update_residual <- function(neuron, footprints, profiles, amplitudes, residual){
  K <- dim(footprints)[1]
  N <- dim(footprints)[2]
  R <- dim(footprints)[3]
  D <- dim(profiles)[3]
  Ts <- dim(residual)[2]
  
  torch_footprints <- ensure3D(torch_tensor(footprints))  # [K, N, R]
  torch_profiles <- ensure3D(torch_tensor(profiles))      # [K, R, D]
  torch_amplitudes <- torch_tensor(amplitudes)            # [K, T]
  
  # Convolve this neuron's amplitude with its weighted temporal factor
  torch_profile_flipped <- torch_flip(torch_profiles[neuron,,], -1)  # [R, D]
  torch_profile <- torch_profile_flipped$unsqueeze(2)                # [R, 1, D]
  torch_amplitude <- torch_amplitudes[neuron,]$reshape(c(1,1,-1))    # [1, 1, T]
  torch_conv <- torch_conv1d(torch_amplitude,
                             torch_profile,
                             padding = D-1)[1,,1:Ts]                 # [R, T]
  
  # Project the result using the neuron's channel factor
  torch_proj <- torch_matmul(torch_footprints[neuron,,], torch_conv)
  proj <- as.array(torch_proj)
  
  # Add that result to the residual
  residual <- residual + proj
  return(residual)
}

downdate_residual <- function(neuron, footprints, profiles, amplitudes, residual){
  K <- dim(footprints)[1]
  N <- dim(footprints)[2]
  R <- dim(footprints)[3]
  D <- dim(profiles)[3]
  Ts <- dim(residual)[2]
  
  torch_footprints <- ensure3D(torch_tensor(footprints))  # [K, N, R]
  torch_profiles <- ensure3D(torch_tensor(profiles))      # [K, R, D]
  torch_amplitudes <- torch_tensor(amplitudes)            # [K, T]
  
  # Convolve this neuron's amplitude with its weighted temporal factor
  torch_profile_flipped <- torch_flip(torch_profiles[neuron,,], -1)  # [R, D]
  torch_profile <- torch_profile_flipped$unsqueeze(2)                # [R, 1, D]
  torch_amplitude <- torch_amplitudes[neuron,]$reshape(c(1,1,-1))    # [1, 1, T]
  torch_conv <- torch_conv1d(torch_amplitude,
                             torch_profile,
                             padding = D-1)[1,,1:Ts]                 # [R, T]
  
  # Project the result using the neuron's channel factor
  torch_proj <- torch_matmul(torch_footprints[neuron,,], torch_conv)
  proj <- as.array(torch_proj)
  
  # Subtract that result to the residual
  residual <- residual - proj
  return(residual)
}

update_amplitude_fast <- function(neuron, footprints, profiles, amplitudes, residual,
                                  noise_std=1.0, amp_rate=5.0){
  K <- dim(footprints)[1]
  N <- dim(footprints)[2]
  R <- dim(footprints)[3]
  D <- dim(profiles)[3]
  Ts <- dim(residual)[2]
  
  new_amplitudes <- array(0, dim=Ts)
  
  # Compute the score by projecting the residual (U_k^T R_k)
  # and cross-correlating with the weighted temporal factor (V_k^T)
  
  # Project the residual onto the channel factor for this neuron
  torch_footprints <- ensure3D(torch_tensor(footprints))  # [K, N, R]
  torch_profiles <- ensure3D(torch_tensor(profiles))      # [K, R, D]
  torch_residual <- torch_tensor(residual)
  torch_proj_residual <- torch_matmul(torch_footprints[neuron]$transpose(1,2),
                                      torch_residual)
  
  # correlate the projected residual with the weighted temporal factor
  torch_score <- torch_conv1d(torch_proj_residual,
                              torch_unsqueeze(torch_profiles[neuron], 1),
                              padding=D-1)[,D:(D+Ts-1)]
  torch_score <- torch_squeeze(torch_score)
  score <- as.array(torch_score)
  
  # Find the peaks in the cross-correlation
  score_peaks <- findpeaks(score,
                           minpeakheight = (noise_std**2)*amp_rate,
                           minpeakdistance = D)
  heights <- score_peaks[,1]
  peaks <- score_peaks[,2]
  
  # Update the amplitudes for this neuron in place
  new_amplitudes[peaks] <- heights
  
  return(new_amplitudes)
}

update_wavelet_factors <- function(neuron, footprints, profiles,
                                   amplitudes, residual, wavelet_rank=1,
                                   warp=function(x) -exp(-x)+1){
  K <- dim(footprints)[1]
  N <- dim(footprints)[2]
  D <- dim(profiles)[3]
  Ts <- dim(residual)[2]
  
  # Check if the factor is used. if not, generate a random target
  if (sum(amplitudes[neuron,]) < 1){
    target <- generate_wavelets(1, N, D, warp)@value               # [1, N, D] 
    torch_target <- torch_tensor(target) %>% torch_squeeze() # [N, D]
  } else{
    # Compute the target (inner product of residual and regressors)
    torch_delay <- torch_eye(D)
    torch_amplitudes <- torch_tensor(amplitudes)
    torch_residual <- torch_tensor(residual)
    torch_regressors <- torch_conv1d(
      torch_reshape(torch_amplitudes[neuron,],shape=c(1,1,Ts)),
      torch_unsqueeze(torch_delay, dim=2) %>% torch_flip(dims=-1),
      padding=D-1)[1,,1:Ts]
    torch_target <- torch_matmul(torch_residual, torch_regressors$transpose(1,2))
  }
  
  # Project the target onto the set of normalized rank-K wavelets
  # and keep only the channel factors (U_n)
  # and the weighted temporal factors (V_n^T)
  USV <- torch_svd(torch_target)
  
  # Truncate the SVD and normalize the singular values
  U <- USV[[1]][ ,1:wavelet_rank] # 2D array
  S <- USV[[2]][  1:wavelet_rank] # 1D array
  V <- USV[[3]][ ,1:wavelet_rank] # 2D array
  
  # Truncate, weight, and transpose the factors as appropriate
  torch_footprint <- U
  torch_profile <- torch_matmul(torch_diag(S) / torch_norm(S),
                                V$transpose(1,2))
  
  # Revert back to array
  new_footprint <- as.array(ensure3D(torch_footprint))
  new_profile <- as.array(ensure3D(torch_profile))
  
  return(list(new_footprint, new_profile))  
}

map_estimate_fast <- function(wavelets, amplitudes, data,
                              noise_std=1.0, amp_rate=5.0, wavelet_rank=1, 
                              warp=function(x) -exp(-x)+1,
                              num_iters=20, tol=1e-06, show_progress=TRUE){
  # Fit the wavelets and amplitudes by maximum a posteriori (MAP) estimation
  K <- dim(wavelets)[1]
  N <- dim(wavelets)[2]
  D <- dim(wavelets)[3]
  Ts <- dim(data)[2]
  
  # Initialize the residual
  torch_amplitudes <- torch_tensor(amplitudes)
  torch_wavelets <- ensure3D(torch_tensor(wavelets))
  torch_wavelets_flip <- torch_flip(torch_wavelets, dims=-1)
  torch_pred <- torch_conv1d(torch_amplitudes,
                             torch_wavelets_flip$permute(c(2,1,3)),
                             padding=D-1)[,1:Ts]
  pred <- as.array(torch_pred) # convert back to array
  residual <- data - pred
  
  # Initialize the wavelet factors
  stopifnot(wavelet_rank >= 1)
  USV <- torch_svd(torch_wavelets)
  U <- USV[[1]][,,1:wavelet_rank] # 3D array
  S <- USV[[2]][ ,1:wavelet_rank] # 2D array
  V <- USV[[3]][,,1:wavelet_rank] # 3D array
  
  # Define footprints/profiles
  torch_footprints <- U
  torch_profiles <- V * torch_unsqueeze(S, 2)
  torch_profiles <- torch_profiles$permute(c(1,3,2))
  
  # Revert back to array
  footprints <- as.array(ensure3D(torch_footprints))
  profiles <- as.array(ensure3D(torch_profiles))
  
  # Track log likelihoods over iterations
  lls <- log_likelihood(wavelets, amplitudes, data, noise_std)
  
  # Coordinate ascent
  if (show_progress){
    pb <- progress_bar$new(
      total = num_iters,
      format = "iter :current/:total [:bar] :percent eta: :eta")
  }
  
  for (iter in 1:num_iters){
    
    # Update neurons one at a time
    for (k in 1:K){
      # Update the residual in place (add a_k \circledast W_k)
      residual <- update_residual(
        k, footprints, profiles, amplitudes, residual)
      
      # Update the wavelet and amplitude with the residual
      amplitudes[k,] <- update_amplitude_fast(
        k, footprints, profiles, amplitudes, residual, noise_std, amp_rate)
      
      new_wavelets <- update_wavelet_factors(
        k, footprints, profiles, amplitudes, residual, wavelet_rank, warp)
      footprints[k,,] <- new_wavelets[[1]]
      profiles[k,,] <- new_wavelets[[2]]
      
      # Downdate the residual in place (subtract a_k \circledast W_k)
      residual <- downdate_residual(
        k, footprints, profiles, amplitudes, residual)
    }
    
    # Reconstruct the wavelets in place
    torch_footprints <- ensure3D(torch_tensor(footprints))
    torch_profiles <- ensure3D(torch_tensor(profiles))
    torch_wavelets <- torch_matmul(torch_footprints, torch_profiles)
    wavelets <- as.array(ensure3D(torch_wavelets))
    
    # Compute the log likelihood
    lls <- append(lls, log_likelihood(wavelets, amplitudes, data, noise_std))
    
    # Check for convergence
    if (abs(lls[iter+1] - lls[iter]) < tol){
      break  
    }
    
    # Step the progress bar
    if (show_progress){
      pb$tick() 
    }
  }
  
  # Check for convergence and warn if necessary
  if (abs(lls[iter+1] - lls[iter]) < tol){
    cat(sprintf("Iteration [%d/%d] Convergence detected!\n",iter,num_iters))
  } else{
    warning(sprintf("Increase num_iters > %d",num_iters))
  }
  return(list(wavelets=wavelets, amplitudes=amplitudes, lls=lls))
}

#---------------#
# Main Function #
#---------------#
#' Fit ConvNMF Model
#'
#' Fits a convolutional NMF model to multi-channel time-series data using MAP estimation.
#'
#' @param data Matrix of dimensions \code{num_channels} × \code{num_samples}.
#' @param sample_freq Sampling frequency.
#' @param num_neurons Number of neurons to be fitted.
#' @param len_wavelet Length of spatial wavelet (in samples).
#' @param warp Function for time warping.
#' @param mean_spikes Expected number of spikes per unit time.
#' @param mean_amplitude Mean of the Gamma distribution for amplitudes.
#' @param shape_amplitude Shape of the Gamma distribution for amplitudes.
#' @param noise_std Standard deviation of Gaussian noise.
#' @param amp_rate Detection threshold in SD units.
#' @param wavelet_rank SVD rank for wavelet decomposition.
#' @param num_iters Maximum number of MAP iterations.
#' @param tol Convergence tolerance.
#' @param weight_wavelets Renormalize and reorder wavelets
#'
#' @return A fitted \code{\link{ConvNMF-class}} object.
#' @export
convNMF <- function(data, # multi-channel time-series array (NxT)
                    sample_freq, # number of samples per unit time
                    num_neurons, # number of neurons
                    len_wavelet, # length of spikes
                    warp=function(x) -exp(-x)+1, # warp function
                    mean_spikes=10, # expected number of spikes per unit time
                    mean_amplitude=15, # alpha/beta of Gamma (alpha, beta)
                    shape_amplitude=3, # alpha of Gamma (alpha, beta)
                    noise_std=1, # MNV(0, noise_std^2)
                    amp_rate=3, # detect amplitudes (amp_rate) times higher than SD
                    wavelet_rank=1, # Dimension reduction in SVD
                    num_iters=50, # number of iterations
                    tol=1e-06, # tolerance level
                    weight_wavelets=FALSE, # renormalize and reorder wavelets
                    show_progress=TRUE # show progress bar
){
  
  # Get data dimensions
  num_channels <- dim(data)[1] # number of channels, N
  num_samples <- dim(data)[2] # number of samples, T
  cat(paste0("Using (N x T) = (",num_channels," x ",num_samples,") data.\n"))
  if(!is.null(rownames(data))) {
    channel_labels <- rownames(data)
  } else{
    channel_labels <- as.character(1:num_channels)
  }
  data <- array(data, dim=c(num_channels, num_samples))
  
  # Make random wavelets
  cat(paste0("Generating (K x N x D) = (",num_neurons," x ",num_channels," x ",len_wavelet,") wavelets.\n"))
  wavelets <- generate_wavelets(num_neurons, num_channels, len_wavelet, warp)@value
  
  # Make zero amplitudes
  cat(paste0("Initializing (K x T) = (",num_neurons," x ",num_samples,") amplitudes.\n"))
  amplitudes <- array(0, dim=c(num_neurons, num_samples))
  
  # Fit the method
  est <- map_estimate_fast(wavelets, amplitudes, data,
                           noise_std, amp_rate, wavelet_rank, warp,
                           num_iters, tol, show_progress)

  # Reorder by spiking times
  times <- apply(est$wavelets, 1, function(x) which(abs(x)==max(abs(x)),arr.ind=T)[2])
  wavelets_r <- array(est$wavelets[order(times),,,drop=FALSE],
                      dim=c(num_neurons, num_channels, len_wavelet))
  amplitudes_r <- array(est$amplitudes[order(times),],
                        dim=c(num_neurons, num_samples))
  
  if(weight_wavelets){
    # Re-normalize wavelets
    w_kn  <- apply(apply(wavelets_r, c(1,2), max), 2, sum)     + 1e-06
    w_k   <- apply(wavelets_r, 1, function(x) quantile(x,.95)) + 1e-06
    w_knd <- array(w_k %*% t(w_kn), dim(wavelets_r))
    wavelets_r_knd <- wavelets_r / w_knd
    weights <- w_knd * apply(wavelets_r_knd, 1, norm, "2")
    wavelets_r <- wavelets_r / weights
    
    # Re-order channels
    chMax <- apply(wavelets_r,1,function(x) apply(x,1,function(x) max(abs(x))))
    # chMaxIdx <- which(chMax > .1, arr.ind=TRUE) # abs threshold
    chMaxIdx <- which(chMax == apply(chMax,1, max), arr.ind=TRUE) # neuron max
    new_channels <- as.data.frame(chMaxIdx) %>%
      mutate(val=chMax[chMaxIdx]) %>%
      arrange(col,-val) %>%
      distinct(row) %>%
      pull(row)
    # new_channels <- unlist(apply(wavelets_r, 1, # per neuron
    #                              function(x) which(apply(x,1,max) > .1))) # channel max > .1
    # new_channels <- unique(c(new_channels, 1:dim(wavelets_r)[2]))
    new_channels <- rev(new_channels) # for plotting
    wavelets_r <- wavelets_r[,new_channels,]
    weights <- weights[,new_channels,]
    data <- data[new_channels,]
  } else{
    weights <- array(1, dim=dim(wavelets_r))
    new_channels <- 1:num_channels
  }
  
  # Create ConvNMF object
  est.ConvNMF <- new("ConvNMF",
                     data=new("Data",
                              value=data,
                              num_channels=num_channels, # N 
                              num_samples=num_samples, # T
                              channel_order=new_channels,
                              channel_labels=channel_labels), 
                     wavelets=new("Wavelets",
                                  value=wavelets_r,
                                  num_neurons=num_neurons, # K
                                  num_channels=num_channels, # N 
                                  len_wavelet=len_wavelet, # D
                                  weights=weights,
                                  channel_order=new_channels,
                                  channel_labels=channel_labels),
                     amplitudes=new("Amplitudes",
                                    value=amplitudes_r, 
                                    num_neurons=num_neurons, # K
                                    num_samples=num_samples), # T
                     lls=est$lls,
                     fit_args = list(
                       data=data,
                       sample_freq=sample_freq,
                       num_neurons=num_neurons,
                       len_wavelet=len_wavelet,
                       warp=warp,
                       mean_spikes=mean_spikes,
                       mean_amplitude=mean_amplitude,
                       shape_amplitude=shape_amplitude,
                       noise_std=noise_std,
                       amp_rate=amp_rate,
                       wavelet_rank=wavelet_rank,
                       num_iters=num_iters,
                       tol=tol,
                       weight_wavelets=weight_wavelets,
                       show_progress=show_progress
                     )
  )
  return(est.ConvNMF)
}
