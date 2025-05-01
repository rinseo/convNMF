#' @import ggplot2
#' @import patchwork

NULL

#--------------------#
# Data Visualization #
#--------------------#
#' Plot Method for Data Objects
#'
#' @param x A \code{\link{Data-class}} object.
#'
#' @return A ggplot2 object.
#' @rdname plot.Data
#' @aliases plot.Data
#' @export
setMethod("plot",
          signature(x="Data"),
          function(x){
            # Retrieve labels
            channel_labels <- x@channel_labels[x@channel_order]
            melt(x@value) %>% # Var1 (n) x Var2 (t) x value
              mutate(value=value/(max(abs(value),na.rm=T)*1.05), # rescale values by 5%
                     value=Var1-1 + value) %>% # add n-1 to value (0,N)
              ggplot(aes(x=Var2,y=value,group=Var1))+
              geom_hline(yintercept=0:(x@num_channels-1),color="grey",linewidth=0.2)+
              geom_path(linewidth=0.2)+ # thinner lines
              # Add plot title
              annotate("rect",
                       xmin=-Inf,xmax=Inf,
                       ymin=x@num_channels,ymax=x@num_channels+1, # (N,N+1)
                       alpha=.5,linewidth=0)+
              annotate("text",
                       x=(1+x@num_samples)/2,
                       y=x@num_channels+0.5,
                       label=paste("'Multi-Channel Time-Series'~bold(X)[italic(N)%*%italic(T)]~'('*italic(N)==",
                                   x@num_channels,"*','~italic(T)==",x@num_samples,"*')'"),
                       parse=TRUE)+
              # Adjust x/y scales
              scale_x_continuous(expand=expansion(mult=c(0,.02)))+ # add 2% at the end
              scale_y_continuous(breaks=0:(x@num_channels-1), # break every integer (0,N)
                                 labels=channel_labels, 
                                 limits=c(-0.5,x@num_channels+1), # (0,N+1)
                                 expand=expansion(mult=c(0,0)))+
              guides(color="none")+
              theme_classic()+
              labs(x=parse(text="Timestamp~italic(t)"),
                   y=parse(text="Channels~italic(n)"))
          })

#' Plot Method for Wavelets Objects
#'
#' @param x A \code{\link{Wavelets-class}} object.
#'
#' @return A ggplot2 object.
#' @rdname plot.Wavelets
#' @aliases plot.Wavelets
#' @export
setMethod("plot",
          signature(x="Wavelets"),
          function(x){
            # Retrieve labels
            channel_labels <- x@channel_labels[x@channel_order]
            melt(x@value) %>% # Var1 (k) x Var2 (n) x Var3 (d) x value
              mutate(value=Var2-1 + value) %>% # add n-1 to value (0,N)
              ggplot(aes(x=Var3,y=value,color=as.factor(Var1), # colors specific to k
                         group=interaction(Var1,Var2)))+ # lines specific to k x n
              geom_hline(yintercept=0:(x@num_channels-1),color="grey",linewidth=0.2)+
              geom_path()+
              # geom_smooth()+
              # Add plot title per neuron
              geom_rect(aes(fill=as.factor(Var1)),
                        xmin=-Inf,xmax=Inf,
                        ymin=x@num_channels,ymax=x@num_channels+1, # (N,N+1)
                        alpha=.5,linewidth=0)+
              geom_text(aes(label=paste0("italic(W)[",Var1,"]"),
                            x=(1+x@len_wavelet)/2,y=x@num_channels+0.5),
                        color="black",parse=TRUE,
                        data=melt(x@value) %>% filter(Var2==1,Var3==1))+
              facet_wrap(~Var1,nrow=1)+
              geom_vline(xintercept=1,color="black", linewidth=1)+ # add border lines
              # Adjust x/y scales
              scale_x_continuous(limits=c(1,x@len_wavelet),
                                 breaks=x@len_wavelet,
                                 expand=c(0,0))+
              scale_y_continuous(breaks=0:(x@num_channels-1), # break every integer (0,N)
                                 labels=channel_labels, 
                                 limits=c(-0.5,x@num_channels+1), # (0,N+1)
                                 expand=expansion(mult=c(0,0)))+
              guides(color="none",fill="none")+
              theme_classic()+
              theme(axis.line.y=element_blank(), # added y border lines manually
                    strip.background = element_blank(), 
                    strip.text = element_blank(),
                    panel.spacing.x = unit(0.2, "lines"))+ # spacing between neurons
              labs(x=parse(text="Delay~italic(d)"),
                   y=parse(text="Channels~italic(n)"))
          })

#' Plot Method for Amplitudes Objects
#'
#' @param x An \code{\link{Amplitudes-class}} object.
#'
#' @return A ggplot2 object.
#' @rdname plot.Amplitudes
#' @aliases plot.Amplitudes
#' @export
setMethod("plot",
          signature(x="Amplitudes"),
          function(x){
            melt(x@value) %>% # Var1 (k) x Var2 (t) x value
              mutate(value=value/(max(value,na.rm=T)*1.05), # rescale values by 5%
                     value=(x@num_neurons-Var1) + value) %>% # change k to K-k
              ggplot(aes(x=Var2, y=value, color=as.factor(Var1)))+
              geom_path()+
              # Adjust x/y scales
              scale_x_continuous(expand=expansion(mult=c(.0,.02)))+ # add 2% at the end
              scale_y_continuous(breaks=0:x@num_neurons, # break every integer
                                 label=function(i) parse(text=paste0("italic(A)[",
                                                                    x@num_neurons-as.integer(i),
                                                                    "]")), # revert back to k
                                 expand=expansion(mult=c(.05,0)))+ # add 5% at zero
              guides(color="none")+
              theme_classic()+
              labs(x=parse(text="Timestamp~italic(t)"),
                   y=parse(text="Amplitudes~italic(A[k])"))
          })

#' Plot Method for ConvNMF Objects
#'
#' @param x A \code{\link{ConvNMF-class}} object.
#' @param spikes Mode of spike overlay. Options are \code{"auto"} or \code{"none"}.
#' @param spike_threshold Threshold for spike sorting when \code{spikes = "auto"}
#'
#' @return A ggplot2 object.
#' @rdname plot.ConvNMF
#' @aliases plot.ConvNMF
#' @export
setMethod("plot",
          signature(x="ConvNMF"),
          function(x, spikes="auto", spike_threshold=0.95){
            K <- x@wavelets@num_neurons
            N <- x@data@num_channels
            if(!is.null(x@fit_args)){
              sample_freq <- x@fit_args$sample_freq
            } else{
              sample_freq <- 1
            }
            # Amplitudes
            pA <- plot(x@amplitudes)+
              theme(axis.title.x = element_blank())
            # Combine separate plots into one
            # Data overlaid with spikes
            if(spikes == "auto"){
              spike_times <- extract_spikes(x, spike_threshold)
              pD <- plot(x@data)+
                geom_path(data=. %>% left_join(spike_times, by=c("Var1"="n","Var2"="t")),
                          aes(x=Var2,y=value,color=as.factor(k)))
            }
            if(spikes == "none"){
              pD <- plot(x@data)
            }
            # Wavelets
            pT <- plot(x@wavelets)+
              theme(axis.title.y = element_blank(),
                    axis.text.y = element_blank())
            # Combine separate plots into one
            pA+plot_spacer()+pD+pT+
              plot_layout(ncol=2,
                          heights=c(K,N+2),
                          width  =c(max(11-K,2),1)) #min=2, max=10
          })
