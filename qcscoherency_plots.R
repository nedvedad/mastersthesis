library(quantspec)
library(ggplot2)

set.seed(123)

# simulated time series
# x_t = epsilon_t, y_t = epsilon_{t-1}
nObs <- 1000
epsilon <- rnorm(nObs+1)
XSimul <- matrix(c(epsilon[2:(nObs+1)], epsilon[1:nObs]^2), ncol = 2)

# takes a matrix, tau and number of frequencies and returns a list of quantile coherency, its sd and frequencies
QCSCoherency <- function(X, tau, nFreq, B = 250, l = 32){
    CR <- quantilePG(X, levels.1 = tau, levels.2 = tau, type = "clipped", type.boot = "mbb", B = B, l = l)
    sPG <- smoothedPG(CR, weight = kernelWeight(W = W1, bw = 0.07))
    
    # get an equally spaced subset of length nFreq from [0, pi)
    freq <- getFrequencies(sPG)
    freq <- freq[which(freq < pi)]
    freqSubset <- c(1, floor((1:(nFreq-1))*length(freq)/nFreq))
    freq <- freq[freqSubset]
    
    coh <- getCoherency(sPG, frequencies = freq)
    sdNaive <- getCoherencySdNaive(sPG, frequencies = freq)
    
    return(list(coh = coh, sdNaive = sdNaive, freq = freq))
}

# parameters
tau<- c(0.05, 0.25, 0.5, 0.75, 0.95) # quantile levels
nFreq <- 32 # number of equally spaced Fourier frequencies from [0, pi)
QSCH <- QCSCoherency(XSimul, tau, nFreq)


###

coh <- QSCH$coh
cohSd <- QSCH$sdNaive
cohFreq <- QSCH$freq

j1 <- 1 # first component of X
j2 <- 2 # second component of X
d <- dim(coh)[2] # number of X_j's
J <- length(cohFreq) # number of frequencies
k1 <- dim(coh)[3] # number of tau_1's
k2 <- dim(coh)[5] # number of tau_2's

plot.x <- rep(cohFreq/(2*pi), d^2*k1*k2)
plot.y <- numeric()
plot.tau1 <- numeric()
plot.tau2 <- numeric()
plot.ciLower <- numeric()
plot.ciUpper <- numeric()

confBandAlpha <- 0.05
plotType <- "Re"

# helper function for decomposition of complex valued vectors
decomposeComplex <- function(z, type){
    switch(type,
           "Re" = return(Re(z)),
           "Im" = return(Im(z)))
}

# construct plot data data frame for ggplot
for(tau1 in 1:k1){
    for(tau2 in 1:k2){
        cohBSMean <- apply(coh[, j1, tau1, j2, tau2, ], 1, mean) # mean of all bootstrap runs
        plot.y <- c(plot.y, cohBSMean)
        plot.tau1 <- c(plot.tau1, rep(paste("tau[1]*'='*", tau[tau1]), J))
        plot.tau2 <- c(plot.tau2, rep(paste0("tau[2]*'='*", tau[tau2]), J))
        plot.ciLower <- c(plot.ciLower, cohBSMean + qnorm(confBandAlpha/2) * cohSd[, j1, tau1, j2, tau2])
        plot.ciUpper <- c(plot.ciUpper, cohBSMean + qnorm(1 - confBandAlpha/2) * cohSd[, j1, tau1, j2, tau2])
    }
}

plotData <- data.frame(x = plot.x, y = decomposeComplex(plot.y, plotType), 
                       tau1 = plot.tau1, tau2 = plot.tau2, 
                       ciLower = decomposeComplex(plot.ciLower, plotType), 
                       ciUpper = decomposeComplex(plot.ciUpper, plotType))

# plot the plot data
ggp <- ggplot(plotData, aes(x = x, y = y, ymin = ciLower, ymax = ciUpper)) +
    geom_ribbon(fill = "grey95") +
    geom_line(aes(y = ciUpper), lty = 2, colour = "grey60", size = 0.25) +
    geom_line(aes(y = ciLower), lty = 2, colour = "grey60", size = 0.25) +
    geom_line() +
    facet_grid(tau2 ~ tau1, labeller = label_parsed) +
    ylim(-1, 1) + 
    xlab(expression(omega*"/"*2*pi)) +
    ylab(paste0("Quantile Coherency (", plotType ,")")) + 
    theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "grey85", size = 0.5, fill = NA),
        strip.background = element_blank()
    )
print(ggp)
