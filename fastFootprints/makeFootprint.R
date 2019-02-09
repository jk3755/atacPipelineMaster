



## Calc cov constant
libSize <- length(chr21_counts[,1])
coverageSize <- length(uq)
libFactor <- libSize / coverageSize


## Make profile
Profile <- colMeans(sites_signals, na.rm = TRUE)/libFactor
Profile <- c(Profile, rev(Profile))

## Plot FP
PWMin <- pwm2pfm(PWM)
plotFootprints(Profile = Profile,
               Mlen = motif_size,
               motif = PWMin,
               newpage = TRUE)
#segmentation = Profile.seg)