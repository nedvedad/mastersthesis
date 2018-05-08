# load data
if(!exists('termSt')){
  termSt <- fread(file=file.path('..', 'data', 'realised_volatility.csv'))
  termSt[, timestamp := as.Date(timestamp)]
}

# vectors of column indices of each type
cYieldLobs <- c(1, grep('yield_lobs', names(termSt)))
cBetaLobs <- c(1, grep('beta\\d{1}_lobs', names(termSt)))
cCloseLobs <- c(1, grep('_close_lobs', names(termSt)))
cYieldRV <- c(1, grep('yield_rv', names(termSt)))
cBetaRV <- c(1, grep('beta\\d{1}_rv', names(termSt)))
cNObsSum <- c(1, grep('_nobs_sum', names(termSt)))

# compute first differences of beta coefficients
cBetaDiffNames <- gsub('_lobs', '_diff', names(termSt)[cBetaLobs[-1]])
termSt[, (cBetaDiffNames) := .SD - shift(.SD), .SDcols=cBetaLobs[-1]]
termSt <- termSt[2:nrow(termSt), ]
cBetaDiff <- c(1, grep('beta\\d{1}_diff', names(termSt)))

# load custom font
if(!('Lato' %in% fonts())){
  font_import(paths=file.path('font'), prompt=F)
  loadfonts() 
}

# custom ggplot2 theme
thesisPlotTheme <- theme(
  text=element_text(family='Lato', color='grey20'),
  plot.margin=unit(c(0, 0, 0, 0), 'cm'),
  plot.title=element_text(size=18, face='bold'),
  plot.subtitle=element_text(face='italic', size=12, colour='grey40'),
  plot.caption=element_text(size=9, colour='grey60', hjust=0.5),
  panel.background=element_rect(fill='grey97', colour='grey97'),
  panel.spacing=unit(0.5, 'cm'),
  panel.grid.major.y=element_line(size=0.1, colour='grey40'),
  panel.grid.major.x=element_blank(),
  panel.grid.minor=element_blank(),
  strip.background=element_blank(),
  axis.ticks=element_blank(),
  axis.title=element_text(size=9),
  legend.position='top',
  legend.title=element_blank()
)

# statistics to be included in summary tables
thesisRequiredSats <- c(
  'nobs', 'Minimum', 'Maximum', '1. Quartile', 
  '3. Quartile', 'Mean', 'Median', 'Variance', 
  'Stdev', 'Skewness', 'Kurtosis'
)

# quantiles required in the quantile cross-spectral plots
thesisRequiredQunatiles <- c(0.5, 0.25, 0.75, 0.05, 0.95)

# compute the fitted values from N-S model
fitYields <- function(betas, date, maturities=seq(1, 30, by=0.5), lambda=0.0609){
  res <- numeric()
  for(m in maturities){
    factors <- t(c(
      1,
      (1 - exp(-lambda * m)) / (lambda * m), 
      (1 - exp(-lambda * m)) / (lambda * m) - exp(-lambda * m)
    ))
    res <- c(res, (sum(betas * factors)))
  }
  data.frame(maturity=maturities, yield=res, type='fitted', date=date)
}
fitYield <- function(beta0, beta1, beta2, m, lambda=0.0609){
  return(
    beta0 + 
      beta1 * ((1 - exp(-lambda * m)) / (lambda * m)) +
      beta2 * ((1 - exp(-lambda * m)) / (lambda * m) - exp(-lambda * m))
  )
}

# spectrum plots
getSpectrumPlotData <- function(x, title, cl=0.95, k.type='modified.daniell', k.m=c(7, 11, 35)) {
  k <- kernel(k.type, m=k.m)
  rawSpec <- spec.pgram(x, kernel=k, plot=FALSE)
  
  specDb <-  10 * log10(rawSpec$spec)
  freq <- rawSpec$freq * 2
  df <- rawSpec$df
  bandwidth <- rawSpec$bandwidth
  
  upper.quantile <- 1 - (1 - cl) * pchisq(df, df, lower.tail = FALSE)
  lower.quantile <- (1 - cl) * pchisq(df, df)
  conf.lim <- 10 * log10(1/(qchisq(c(upper.quantile, lower.quantile), df)/df))
  
  resDF <- data.frame(freq = freq, spec = specDb, ci_lower = specDb + conf.lim[1], ci_upper = specDb + conf.lim[2])
  list(data = resDF, df = df, bandwidth = bandwidth, title = title, kernel = k.type, kernel.dim = k.m, cl = cl)
}

plotSpectrum <- function(plotData) {
  ggplot(plotData$data, aes(x=freq, y=spec, ymin=ci_lower, ymax=ci_upper)) +
    geom_ribbon(alpha=0.15) +
    geom_line(size=0.3) +
    labs(
      title='Power Spectrum',
      subtitle=paste0('of ', plotData$title, ' series'),
      caption=paste0('bandwidth: ', round(plotData$bandwidth, 5), 
                     ', smooth. kernel: ', plotData$kernel, 
                     ' (', paste(plotData$kernel.dim, collapse=', '), ')',
                     ', conf. int.: ', plotData$cl*100, '%'),
      x=expression(paste('frequency (', pi, ')')),
      y='spectrum (dB)'
    ) +
    thesisPlotTheme
}

# sliding-window spectrum plots
getSWSpectrumPlotData <- function(x, title, cl=0.95, k.type='modified.daniell', k.m=c(11), win.length=120){
  x <- as.numeric(unlist(x))
  k <- kernel(k.type, m=k.m)
  
  freqAll <- numeric()
  specDbAll <- numeric()
  indexAll <- integer()
  pValAll <- numeric()
  
  for(i in 1:(pmax(1, length(x)-win.length+1))){
    xSub <- x[i:(i + win.length - 1)]
    # remove linear trend
    fit <- lm(xSub ~ seq_along(xSub))
    xSubDT <- fit$residuals
    rawSpec <- spec.pgram(xSubDT, kernel=k, plot=FALSE, detrend = FALSE)
    specDbAll <- c(specDbAll, 10 * log10(rawSpec$spec))
    freqAll <- c(freqAll, rawSpec$freq * 2)
    indexAll <- c(indexAll, rep(i, win.length))
    pVal <- suppressWarnings(adf.test(xSubDT)$p.value)
    pValAll <- c(pValAll, pVal)
  }
  
  resDT <- data.table(freq = freqAll, spec = specDbAll, index = indexAll)
  resQuant <- resDT[, .(spec_025 = quantile(spec, 0.25),
                        spec_050 = quantile(spec, 0.50),
                        spec_075 = quantile(spec, 0.75)), by = 'freq']
  list(data = resDT, dataQuant = resQuant, df = df, title = title, 
       kernel = k.type, kernel.dim = k.m, win.length = win.length, p.vals = pValAll)
}

plotSWSpectrum <- function(plotData, type){
  if(type == 'spectogram'){
    # spectogram
    ggplot(plotData$data, aes(x=freq, y=index, fill=spec)) + 
      geom_raster() + 
      scale_fill_gradient(low="black", high="white") +
      labs(
        title='Spectrogram',
        subtitle=paste0('for sliding window samples of ', plotData$title, ' series'),
        caption=paste0('smooth. kernel: ', plotData$kernel, 
                       ' (', paste(plotData$kernel.dim, collapse=', '), ')',
                       ', window length: ', plotData$win.length),
        x=expression(paste('frequency (', pi, ')')),
        y='index',
        fill='spectrum (dB)'
      ) +
      thesisPlotTheme
  } else if(type == 'samples-quants'){
    # quantiles of estimated spectra
    ggplot(plotData$dataQuant, aes(x=freq, y=spec_050, ymin=spec_025, ymax=spec_075)) +
      geom_ribbon(alpha=0.15) +
      geom_line(size=0.3) +
      labs(
        title='Median Power Spectrum',
        subtitle=paste0('of ', plotData$title, ' series'),
        caption=paste0('smooth. kernel: ', plotData$kernel, 
                       ' (', paste(plotData$kernel.dim, collapse=', '), ')',
                       ', bands: lower and upper quartiles'),
        x=expression(paste('frequency (', pi, ')')),
        y='spectrum (dB)'
      ) +
      thesisPlotTheme
  } else if(type == 'samples-spectra'){
    # random set of sliding window samples
    set.seed(100)
    randomObs <- sample(1:max(plotData$data$index), size = 16)
    
    # plot of spectrum of the random set of sliding window samples
    ggplot(plotData$data[index %in% randomObs, ], aes(x=freq, y=spec)) +
      geom_line(size=0.25) +
      facet_wrap(~ index) +
      labs(
        title='Power Spectrum',
        subtitle=paste0('for selected samples of series', plotData$title),
        caption=paste0('smooth. kernel: ', plotData$kernel, 
                       ' (', paste(plotData$kernel.dim, collapse=', '), ')',
                       ', window length: ', plotData$win.length),
        x=expression(paste('frequency (', pi, ')')),
        y='spectrum (dB)'
      ) +
      thesisPlotTheme
  } else if(type == 'pvals'){
    # plot of p-values from the ADF test of sliding window sampels
    ggplot(mapping=aes(x=1:length(plotData$p.vals), y=plotData$p.vals)) +
      geom_line(size=0.25) +
      labs(
        title='P-Values',
        subtitle=paste0('of Augmented Dickey-Fuller test for sliding window samples of ', plotData$title),
        caption=paste0('window length: ', plotData$win.length), 
        x='index',
        y='p-value'
      ) +
      thesisPlotTheme
  }
}

# cross-spectrum plots
getCrossSpectrumPlotData <- function(x, titles, cl=0.95, k.type='modified.daniell', k.m=c(3, 9, 21)) {
  k <- kernel(k.type, m=k.m)
  rawSpec <- spec.pgram(x, kernel=k, plot=FALSE)
  
  coh <-  rawSpec$coh
  phase <- rawSpec$phase
  freq <- rawSpec$freq
  df <- rawSpec$df
  bandwidth <- rawSpec$bandwidth
  
  gg <- 2/df
  se <- sqrt(gg/2)
  z <- -qnorm((1 - cl)/2)
  
  cohLim <- pmin(0.99999, sqrt(coh))
  cohCIUpper <- (tanh(atanh(cohLim) + z * se))^2
  cohCILower <- (pmax(0, tanh(atanh(cohLim) - z * se)))^2
  
  phaseLim <- asin(pmin(0.9999, qt(cl, 2/gg - 2) * sqrt(gg * (sqrt(coh)^{-2} - 1)/(2 * (1 - gg)))))
  phaseCIUpper <- phase + phaseLim
  phaseCILower <- phase - phaseLim
  
  resDF <- data.frame(freq = freq * 2, coh = coh, phase = phase,
                      coh_ci_lower = cohCILower, coh_ci_upper = cohCIUpper,
                      phase_ci_lower = phaseCILower, phase_ci_upper = phaseCIUpper)
  list(data = resDF, df = df, bandwidth = bandwidth, titles = titles, kernel = k.type, kernel.dim = k.m)
}

plotCrossSpectrum <- function(plotData){
  ggplot(plotData$data, aes(x=freq, y=coh, ymin=coh_ci_lower, ymax=coh_ci_upper)) +
    geom_ribbon(alpha=0.15) +
    geom_line(size=0.3) +
    coord_cartesian(ylim=c(0, 1)) +
    labs(
      title='Squared Coherency',
      subtitle=paste0('between ', paste(plotData$title, collapse=' and '), ' series'),
      caption=paste0('bandwidth: ', round(plotData$bandwidth, 5), 
                     ', smoothing kernel: ', plotData$kernel, 
                     ' (', paste(plotData$kernel.dim, collapse=', '), ')'),
      x=expression(paste('frequency (', pi, ')')),
      y='squared coherency'
    ) +
    thesisPlotTheme +
    theme(plot.margin = unit(c(0, 0.1, 0, 0), 'cm')) ->
    cohPlot
  
  ggplot(plotData$data, aes(x=freq, y=phase, ymin=phase_ci_lower, ymax=phase_ci_upper)) +
    geom_ribbon(alpha=0.15) +
    geom_line(size=0.3) +
    coord_cartesian(ylim=c(-3, 3)) +
    labs(
      title='Phase',
      subtitle=paste0('between ', paste(plotData$title, collapse=' and '), ' series'),
      caption=paste0('bandwidth: ', round(plotData$bandwidth, 5), 
                     ', smoothing kernel: ', plotData$kernel, 
                     ' (', paste(plotData$kernel.dim, collapse=', '), ')'),
      x=expression(paste('frequency (', pi, ')')),
      y='phase'
    ) +
    thesisPlotTheme +
    theme(plot.margin = unit(c(0, 0, 0, 0.1), 'cm')) ->
    phasePlot
  
  grid.arrange(cohPlot, phasePlot, nrow=1, padding = unit(1, 'cm'))
}

# quantile cross-spectrum coherency plots
getQuantileCoherencyPlotData <- function(x, titles, freqCut=(0:127)/128, cl=0.95) {
  w = kernelWeight(W=W1, b=0.5*nrow(series)^(-1/4))
  sPG <- smoothedPG(series, levels.1=thesisRequiredQunatiles, weight=w)
  coh <- getCoherency(sPG, frequencies=(2 * pi * freqCut))
  cohCI <- getPointwiseCIs(sPG, quantity='coherency', frequencies=(2 * pi * freqCut), alpha=1-cl)
  
  tau1 <- numeric()
  tau2 <- numeric()
  cohRe <- numeric()
  cohReCILower <- numeric()
  cohReCIUpper <- numeric()
  freq <- numeric()
  
  for(i in 1:length(thesisRequiredQunatiles)){
    for(j in 1:length(thesisRequiredQunatiles)){
      cohRe <- c(cohRe, Re(coh[1:(length(freqCut)/2), 1, i, 2, j, 1]))
      cohReCILower <- c(cohReCILower, Re(cohCI$lowerCIs[1:(length(freqCut)/2), 1, i, 2, j]))
      cohReCIUpper<- c(cohReCIUpper, Re(cohCI$upperCIs[1:(length(freqCut)/2), 1, i, 2, j]))
      tau1 <- c(tau1, rep(thesisRequiredQunatiles[i], length(freqCut)/2))
      tau2 <- c(tau2, rep(thesisRequiredQunatiles[j], length(freqCut)/2))
      freq <- c(freq, freqCut[1:(length(freqCut)/2)] * 2)
    }
  }
  
  res <- data.frame(coh = cohRe, coh_ci_lower = cohReCILower, coh_ci_upper = cohReCIUpper, freq = freq, tau1 = tau1, tau2 = tau2)
  list(data=res, titles=titles, cl=cl)
}

plotQuantileCoherency <- function(plotData){
  ggplot(plotData$data, aes(x=freq, y=coh, ymin=coh_ci_lower, ymax=coh_ci_upper)) +
    geom_ribbon(alpha=0.15) +
    geom_line(size=0.25) +
    facet_grid(tau1 ~ tau2, labeller=label_both) +
    coord_cartesian(ylim=c(-1, 1)) +
    labs(
      title='Quantile Coherency',
      subtitle=paste0('between series ', plotData$titles[1], ' and ', plotData$titles[2]),
      caption=paste0('conf. int.: ', plotData$cl*100, '%'),
      x=expression(paste('frequency (', pi, ')')),
      y='quantile coherency (Re)'
    ) +
    thesisPlotTheme -> cohPlot
  print(cohPlot)
}
