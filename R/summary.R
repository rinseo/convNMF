#' @importFrom tidyr pivot_longer expand_grid
#' @importFrom stringr str_extract

NULL

# Set generics
setGeneric("extract_spikes", function(object, ...) standardGeneric("extract_spikes"))
setGeneric("power_explained", function(object) standardGeneric("power_explained"))
setGeneric("sequenciness", function(object, ...) standardGeneric("sequenciness"))
setGeneric("diss", function(object1, object2) standardGeneric("diss"))
setGeneric("compute_diss_matrix", function(object, ...) standardGeneric("compute_diss_matrix"))
setGeneric("confusion_matrix", function(object1, object2) standardGeneric("confusion_matrix"))

#---------------#
# Model Summary #
#---------------#
#' Extract Spikes from ConvNMF Model
#'
#' Identifies spike occurrences based on the fitted wavelets and amplitudes.
#'
#' @param object A fitted \code{\link{ConvNMF-class}} object.
#' @param spike_threshold Quantile threshold for spike sorting.
#'
#' @return A tibble of detected spikes (\code{k}=neuron, \code{n}=channel, \code{t}=time).
#' @name extract_spikes
#' @export
setMethod("extract_spikes",
          signature(object="ConvNMF"),
          function(object, spike_threshold=0.95){
            # Find spiking channels
            if(length(object@wavelets@weights)>0){
              wavelets <- object@wavelets@value * object@wavelets@weights
            } else{
              wavelets <- object@wavelets@value
            }
            spike_channels <- melt(wavelets) %>% # Var1 (k) x Var2 (n) x Var3 (d) x value
              rename(k=Var1, n=Var2, d=Var3) %>% 
              group_by(n) %>% # per channel 
              filter(abs(value)>=quantile(abs(value),spike_threshold)) %>% 
              distinct(k,n)
            
            # Identify spiking times
            amplitudesT <- data.frame(t(object@amplitudes@value))
            names(amplitudesT) <- paste0("A",1:object@amplitudes@num_neurons)
            dataT <- data.frame(t(object@data@value)) %>% 
              cbind(amplitudesT)
            
            # Add predictions per channel
            spike_times <- dataT %>% 
              mutate(t=row_number() %>% as.integer()) %>% 
              pivot_longer(cols=matches("A\\d{1,}"),names_to="k",values_to="value_A") %>%
              pivot_longer(cols=matches("X\\d{1,}"),names_to="n",values_to="value_X") %>%
              filter(value_A>0) %>% 
              # Add durations of a spike
              expand_grid(data.frame(d=0:(object@wavelets@len_wavelet-1))) %>% 
              mutate(t=t+d,
                     k=str_extract(k,"\\d{1,}") %>% as.integer(),
                     n=str_extract(n,"\\d{1,}") %>% as.integer()) %>%
              inner_join(spike_channels,by=c("k","n")) %>% 
              distinct(k,n,t)
            
            # Return Var1 (k) x Var2 (n) x Var3 (d)
            return(spike_times)
          })

#' Log Likelihood of ConvNMF Model
#'
#' Returns the (un-normalized) log likelihood of the fitted ConvNMF model.
#'
#' @param object A fitted \code{\link{ConvNMF-class}} object.
#'
#' @return A numeric scalar representing the log likelihood.
#' @rdname logLik.ConvNMF
#' @aliases logLik.ConvNMF
#' @export
setMethod("logLik",
          signature(object="ConvNMF"),
          function(object){
            if(length(object@lls)==0){
              stop("ConvNMF object is not fitted. Fit the model!")
            }
            N <- object@data@num_channels
            Ts <- object@data@num_samples
            ll <- object@lls[length(object@lls)]*N*Ts
            return(ll)
          })

#' Reconstruct Data from ConvNMF Model
#'
#' Computes the predicted data from the fitted wavelets and amplitudes.
#'
#' @param object A fitted \code{\link{ConvNMF-class}} object.
#'
#' @return A matrix representing the predicted data.
#' @rdname predict.ConvNMF
#' @aliases predict.ConvNMF
#' @export
setMethod("predict",
          signature(object="ConvNMF"),
          function(object){
            N <- object@data@num_channels
            D <- object@wavelets@len_wavelet
            Ts <- object@data@num_samples
            # Compute the model prediction by convolving the amplitude and wavelets
            torch_amplitudes <- object@amplitudes@value
            if(length(object@wavelets@weights)>0){
              wavelets <- object@wavelets@value * object@wavelets@weights
            } else{
              wavelets <- object@wavelets@value
            }
            torch_wavelets <- torch_tensor(wavelets)
            torch_wavelets_flip <- torch_flip(torch_wavelets, dims=-1)
            torch_pred <- torch_conv1d(torch_amplitudes,
                                       torch_wavelets_flip$permute(c(2,1,3)),
                                       padding=D-1)[,1:Ts]
            pred <- as.array(torch_pred)
            data_pred <- new("Data",
                             value=pred,
                             num_channels=N,
                             num_samples=Ts,
                             channel_order=object@data@channel_order,
                             channel_labels=object@data@channel_labels)
            return(data_pred)
          })

#' Update ConvNMF Model with New Parameters
#'
#' Refits the ConvNMF model using the updated model parameters.
#'
#' @param object A \code{\link{ConvNMF-class}} object.
#' @param ... Named arguments passed to \code{\link{convNMF}} to update the model formula.
#'
#' @return A newly fitted \code{ConvNMF} object.
#' @rdname update.ConvNMF
#' @aliases update.ConvNMF
#' @export
setMethod("update",
          signature(object = "ConvNMF"),
          function(object, ...) {
            new_args <- list(...)
            args <- modifyList(object@fit_args, new_args)
            do.call(convNMF, args)
          })

#------------#
# Evaluation #
#------------#
#' Power Explained by ConvNMF Model
#'
#' Calculates the percent power explained by a factorization.
#'
#' @param object A fitted \code{\link{ConvNMF-class}} object.
#'
#'#' @references
#' Mackevicius, E. L., Bahle, A. H., Williams, A. H., Gu, S., Denisenko, N. I., Goldman, M. S., & Fee, M. S. (2019).
#' Unsupervised discovery of temporal sequences in high-dimensional datasets, with applications to neuroscience.
#' \emph{eLife, 8}, e38471. \doi{10.7554/eLife.38471}
#' 
#' @return A numeric scalar representing the \eqn{1 - \frac{\Sigma(\bold{X} - \hat{\bold{X}})^2}{\bold{X}^2}} value.
#' @name power_explained
#' @export
setMethod("power_explained",
          signature(object="ConvNMF"),
          function(object){
            data <- object@data@value
            pred <- predict(object)@value
            residual <- data - pred
            power <- 1 - sum(residual^2) / sum(data^2)
            return(power)
          })

#' Compute Sequenciness Score for ConvNMF Model
#'
#' Estimates the "sequenciness" of a fitted \code{\link{ConvNMF-class}} object by comparing the power explained
#' by the model on the actual data versus shuffled data.
#'
#' @param object A fitted \code{\link{ConvNMF-class}} object.
#' @param num_shuffle Integer. Number of shuffled datasets to generate. Defaults to 20.
#' @param num_iters Integer. Number of model fits to run per shuffle. Defaults to 20.
#'
#' @references
#' Mackevicius, E. L., Bahle, A. H., Williams, A. H., Gu, S., Denisenko, N. I., Goldman, M. S., & Fee, M. S. (2019).
#' Unsupervised discovery of temporal sequences in high-dimensional datasets, with applications to neuroscience.
#' \emph{eLife, 8}, e38471. \doi{10.7554/eLife.38471}
#' 
#' @return A numeric value representing the sequenciness score. Values closer to 1 indicate stronger sequential patterns.
#' @name sequenciness
#' @export
setMethod("sequenciness",
          signature(object="ConvNMF"),
          function(object, num_shuffle=20, num_iters=20){
            if(num_shuffle < 20) warning("num_shuffle < 20")
            if(num_iters < 20) warning("num_iters < 20")
            
            Ts <- object@data@num_samples
            D <- object@wavelets@len_wavelet
            
            # Calculate the median power among multiple shuffles
            # num_shuffle <- 20
            power_full_max <- rep(0, num_shuffle)
            power_col_max <- rep(0, num_shuffle)
            pb <- progress_bar$new(
              total = num_shuffle,
              format = "shuffle :current/:total [:bar] :percent eta: :eta")
            pb$tick(0)
            
            for (j in 1:num_shuffle){
              shuffle_full <- sample(1:Ts, Ts)
              shuffle_col <- c(1:D, sample((D+1):Ts, Ts-D))
              
              # For each shuffling
              # Find the best performing model among multiple iterations
              # num_iters <- 20
              power_full <- rep(0, num_iters)
              power_col <- rep(0, num_iters)
              capture.output({
                for(i in 1:num_iters){
                  # Fully shuffled
                  object_shuffled_full <- update(object,
                                                 data=object@data@value[,shuffle_full],
                                                 show_progress=F)                  
                  power_full[i] <- power_explained(object_shuffled_full)
                  # Column shuffled
                  object_shuffled_col <- update(object,
                                                data=object@data@value[,shuffle_col],
                                                show_progress=F)
                  power_col[i] <- power_explained(object_shuffled_col)
                }
              })
              power_full_max[j] <- max(power_full)
              power_col_max[j] <- max(power_col)
              
              # Step the progress bar
              pb$tick()
            }
            
            # Compare with the current power
            power <- power_explained(object)
            seq_score <- (power - median(power_col_max)) / (power - median(power_full_max))
            return(seq_score)
          })

# Calculate the dissimilarity score
reconstruct <- function(object, flatten=TRUE) {
  if(length(object@wavelets@weights)>0){
    wavelets <- object@wavelets@value * object@wavelets@weights
  } else{
    wavelets <- object@wavelets@value
  }
  # adding Dim 1 when num_neurons = 1
  if(length(dim(wavelets))==2) wavelets <- array(wavelets, dim=c(1,dim(wavelets)))
  amplitudes <- object@amplitudes@value
  K <- dim(wavelets)[1]
  D <- dim(wavelets)[3]
  Ts <- dim(amplitudes)[2]
  torch_conv_list <- lapply(1:K, function(k){
    W_k <- torch_tensor(wavelets[k,,,drop=F])$permute(c(2,1,3)) # [N, 1, D]
    A_k <- torch_tensor(amplitudes[k,])$reshape(c(1,1,-1))      # [1, 1, T]
    torch_conv1d(A_k, W_k$flip(dims=-1), padding=D-1)[1,,1:Ts]  # [N, T]
  })                                                                 
  torch_conv <- torch_stack(torch_conv_list)                    # [K, N, T]   
  if (flatten){
    torch_conv <- torch_conv$reshape(c(K, -1))                  # [K, N*T] 
  }
  conv <- as.array(torch_conv)
  return(conv)
}

#' Compute Dissimilarity Between Two ConvNMF Models
#'
#' Calculates the dissimilarity score between two \code{\link{ConvNMF-class}} models
#' by comparing their reconstructed factor matrices.
#'
#' @param object1 A fitted \code{\link{ConvNMF-class}} object (first model).
#' @param object2 A fitted \code{\link{ConvNMF-class}} object (second model).
#'
#' @references
#' Mackevicius, E. L., Bahle, A. H., Williams, A. H., Gu, S., Denisenko, N. I., Goldman, M. S., & Fee, M. S. (2019).
#' Unsupervised discovery of temporal sequences in high-dimensional datasets, with applications to neuroscience.
#' \emph{eLife, 8}, e38471. \doi{10.7554/eLife.38471}
#' 
#' @return A numeric scalar indicating the dissimilarity score between the two models.
#' Values closer to 0 indicate high similarity or stability of factorizations.
#' @name diss
#' @export
setMethod("diss",
          signature(object1="ConvNMF", object2="ConvNMF"),
          function(object1, object2) {
            # Prepare 2D reconstructions per object
            X1 <- reconstruct(object1, flatten=T) # [K, T]
            X2 <- reconstruct(object2, flatten=T) # [K, T]
            if(!all(dim(X1),dim(X2))) {warning("Two shapes should match!")}
            
            # Compute Frobenius norm (sqrt sum of squares)
            norm1 <- sqrt(apply(X1^2, 1, sum))    # [K]
            norm2 <- sqrt(apply(X2^2, 1, sum))    # [K]  
            
            # Denominator: outer product of norms
            denom <- norm1 %*% t(norm2)           # [K × K]
            # Numerator: inner products between all rows
            numer <- X1 %*% t(X2)                 # [K × K]
            
            # Final similarity matrix
            sim_matrix <- numer / (denom + 1e-06)
            sim_max <- apply(sim_matrix, 1, max, na.rm=T) # per neuron
            diss_score <- 1 - mean(sim_max, na.rm=T) # 1 - avg best similarity
            return(diss_score)
          })

#' Estimate ConvNMF Stability via Repeated Dissimilarity Scores
#'
#' Re-fits the given \code{ConvNMF} model multiple times and computes dissimilarity scores
#' comparing the original model to each new fit.
#'
#' @param object A fitted \code{\link{ConvNMF-class}} object. Used to define the data and hyperparameters.
#' @param num_iters Integer. Number of re-fits to perform. Defaults to 20.
#' @param ... Additional arguments passed to \code{update} function.
#'
#' @references
#' Mackevicius, E. L., Bahle, A. H., Williams, A. H., Gu, S., Denisenko, N. I., Goldman, M. S., & Fee, M. S. (2019).
#' Unsupervised discovery of temporal sequences in high-dimensional datasets, with applications to neuroscience.
#' @seealso \code{\link{diss}}, \code{\link{update}}
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{diss_matrix}}{A symmetric matrix of pairwise dissimilarity scores.}
#'     \item{\code{mean_diss}}{The mean dissimilarity across all model pairs.}
#'   }
#' @name compute_diss_matrix
#' @export
setMethod("compute_diss_matrix",
          signature(object="ConvNMF"),
          function(object, num_iters=20, ...){
            
            # Calculate the dissimilarity matrix
            pb <- progress_bar$new(
              total = num_iters,
              format = "iter :current/:total [:bar] :percent eta: :eta")
            pb$tick(0)
            
            # Step 1: Fit models with different seeds
            models <- vector("list", num_iters)
            for (i in 1:num_iters) {
              capture.output({
                models[[i]] <- update(object, show_progress=F, weight_wavelets=F, ...)
              })
              pb$tick()
            }
            
            # Step 2: Compute pairwise dissimilarities
            diss_matrix <- matrix(NA_real_, num_iters, num_iters)
            for (i in 1:(num_iters - 1)) {
              for (j in (i + 1):num_iters) {
                diss_val <- diss(models[[i]], models[[j]])
                diss_matrix[i, j] <- diss_val
                diss_matrix[j, i] <- diss_val
              }
            }
            
            # Step 3: Compute the average diss score
            mean_diss <- mean(diss_matrix[upper.tri(diss_matrix)], na.rm = TRUE)
            
            return(list(
              diss_matrix = diss_matrix,
              mean_diss = mean_diss
            ))
          })

#' Compute Confusion Matrix for ConvNMF Model
#'
#' Compares ground-truth and predicted spike labels to compute a confusion matrix.
#'
#' @param object1 A synthetic \code{\link{ConvNMF-class}} object with true labels using \code{\link{generate}} function.
#' @param object2 A fitted \code{\link{ConvNMF-class}} object with predicted labels after fitting \code{\link{convNMF}} model.
#' @param spikes Default to \code{"auto"} that uses \code{\link{extract_spikes}} function.
#'
#' @references
#' Buccino, A. P., Hurwitz, C. L., Garcia, S., Magland, J., Siegle, J. H., Hurwitz, R., & Hennig, M. H. (2020).
#' SpikeInterface, a unified framework for spike sorting.
#' \emph{eLife, 9}, e61834. \doi{10.7554/eLife.61834}
#' 
#' @return A 2x2 confusion matrix.
#' @name confusion_matrix
#' @export
setMethod("confusion_matrix",
          signature(object1="ConvNMF",object2="ConvNMF"),
          function(object1, object2){
            N <- object2@data@num_channels
            Ts <- object2@data@num_samples
            
            # Object 1 = reference
            ref <- extract_spikes(object1) %>% mutate(n_t=paste0(n,"_",t)) %>% pull(n_t)
            # Object 2 = prediction
            pred <- extract_spikes(object2) %>% mutate(n_t=paste0(n,"_",t)) %>% pull(n_t)
            
            # True/False Positive/Negative
            TP <- length(intersect(ref, pred))
            FN <- length(ref) - TP
            FP <- length(pred) - TP
            TN <- N*Ts - TP - FN - FP
            confusion_matrix <- matrix(c(TP,FN,FP,TN),nrow=2)
            rownames(confusion_matrix) <- c("(Reference) Positive","            Negative")
            colnames(confusion_matrix) <- c("  Positive","  Negative")
            
            # Calculate scores
            accuracy <- TP/(TP+FN+FP)
            recall <- TP/(TP+FN)
            precision <- TP/(TP+FP)
            miss <- FN/(TP+FN)
            false_discovery <- FP/(TP+FP)
            
            # Print outputs
            cat("                     (Prediction)\n")
            print(confusion_matrix)
            cat("------------------------------------------\n")
            cat(sprintf("            Accuracy: TP/(TP+FN+FP) = %.2f",TP/(TP+FN+FP)),"\n")
            cat(sprintf("              Recall: TP/(TP+FN)    = %.2f",TP/(TP+FN)),"\n")
            cat(sprintf("           Precision: TP/(TP+FP)    = %.2f",TP/(TP+FP)),"\n")
            cat(sprintf("           Miss rate: FN/(TP+FN)    = %.2f",FN/(TP+FN)),"\n")
            cat(sprintf("False discovery rate: FP/(TP+FP)    = %.2f",FP/(TP+FP)),"\n")
            invisible(confusion_matrix)
          })
