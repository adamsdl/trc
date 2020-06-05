# trc_code.R
#
# performs transform-roughness correlation using maximal overlap discrete wavelet transform (MODWT) and
# roughness correlation from Forooghi et al. (2017)
#
# continuous wavelet transform (CWT) is also applied for comparison
#
# 2020-FEB DLA
#
################################################################################################################################

# set working directory

setwd("C:/Users/David/Desktop/University/Sweave/Papers/ESURF_Drag_size_distribution/code_data/standalone")

#setwd("yourworkingdirectory")


library(data.table) # for fread function, superior for importing large datasets
library(pracma) # for detrend function

res = 0.001 # resolution of DEM in m

elev <- fread("elev.txt", header = TRUE) # import thalweg elevation profiles for each DEM

# define matrix to list function
mat2lst <- function(m,rowcol=1) {
  if (rowcol == 1) m <- t(m)
  dm <- as.data.frame(m)
  as.list(dm)
}

elev <- mat2lst(elev) # convert to list
elev <- lapply(elev, function(x) x[x != 999]) # remove 999s used to pad profiles of different lengths

# define despike function to remove values greater or lesser than than 0.5* the previous
despike <- function(x){
  
  x <- x[!is.na(x)] # ensure vector has no NA values
  x.diff = abs(c(diff(x, lag = 1), 0)) # make lagged difference
  x.new <- as.numeric(x[x.diff < (0.5 * sd(x))]) # delete values where the successive difference is greater than 0.5*SD of x
  
  return(x.new)
  
}

elev <- lapply(elev, despike) # use this function a few times in case of consecutive erroneous values
elev <- lapply(elev, despike)
elev <- lapply(elev, despike) 

# normalize elev profiles by the mean
elev.ref <- mean(elev[[1]])
elev <- lapply(elev, function(x) x - elev.ref)


################################################################################################################################


# define function to perform MODWT, where y is a numeric vector and ps is the point-spacing
dwt.transform.thalweg <- function(y, ps = res){
  
  library(waveslim) # for MODWT function
  
  # remove NAs in case they are at the end of vector, and ensure that vector is numeric
  y <- as.numeric(y[!is.na(y)])
  
  # determine the maximum number of wavelengths that can be extracted using MODWT
  n.levels.y <- floor(log(length(y), 2))
  
  # perform MODWT using Daubechies 4 wavelet
  wc = mra(y, wf = "d4", J = n.levels.y, boundary = "reflection", method = "modwt")
  
  # determine which spatial scales the wavelengths correspond to
  start.wl <- ps * 2 # starting wavelength
  scale <- start.wl*2^(0:n.levels.y)
  
  list(scale = scale, wc = wc)
  
}

# apply MODWT function and extract wavelengths and wavelet coefficients
dwt <- lapply(elev, dwt.transform.thalweg, ps = res)
dwt.scale <- lapply(dwt, `[[`, 1)
dwt.coef <- lapply(dwt, `[[`, 2) 


################################################################################################################################


# define function to perform continuous wavelet transform (CWT)
cwt.transform.thalweg <- function(y, ps = res){
  
  require(biwavelet) # for CWT function
  
  # remove NAs - in the case that they are at the end of the vector
  y <- y[!is.na(y)]
  
  # continuous wavelet transform require x coordinates (i.e. distance along profile)
  len = length(y) # find number of data points
  x <- seq(from = ps, to = len * ps, by = ps)
  
  # perform CWT using morlet
  w <- wt(cbind(x,y), pad = FALSE, mother = "morlet", do.sig = FALSE)
  
  # extract wavelengths and coefficients
  scale <- w[[7]]
  wc <- Re(w[["wave"]])
  
  list(scale = scale, wc = wc)
  
}

# apply CWT function and extract wavelengths and wavelet coefficients
cwt <- lapply(elev, cwt.transform.thalweg, ps = res)
cwt.scale <- lapply(cwt, `[[`, 1)
cwt.coef <- lapply(cwt, `[[`, 2)


################################################################################################################################


# define function to perform roughness correlation
roughness.correlation <- function(y, ps = res){
  
  library(fBasics) # for skewness function
  
  # define effective slope (Napoli et al. 2008) function 
  effective.slope <- function(y, dx = ps){
    n = length(y)
    y1 = y[1:(n-1)]
    y2 = y[2:n]
    dy = y2-y1
    m = abs(dy)/dx
    return(sum(m)/(n-1))
  }
  
  # calculate standard deviation, skewness, 4.4 standard deviations above the mean, and effective slope
  sd <- sd(y)
  sk <- skewness(y)
  sd4.4 = sd * 4.4 # calculate 4.4 standard deviations above the mean, used as kz in Forooghi et al. (2017)
  es <- effective.slope(y, dx = ps)
  
  # perform Forooghi et al. (2017) roughness correlation
  F1 <- 0.67 * sk^2 + 0.93 * sk + 1.3
  F2 <- 1.05 * (1 - exp(-3.8 * es))
  ksk <- F1 * F2
  
  # calculate Nikuradse equivalent sand roughness ks
  ks <- ksk * sd
  
  stats <- cbind(sd, sk, sd4.4, es, ks, ksk)
  
  colnames(stats) <- c("sd", "sk", "4.4.sd", "es", "ks", "ksk")
  
  return(stats)
  
  
}

# apply roughness correlation function to elevation profiles
# the profiles must be initially detrended, unlike the wavelet transform
elev.roughness.correlation <- as.data.frame(do.call(rbind, lapply(lapply(elev, detrend, tt = 'linear'), roughness.correlation)))

# define nested function to apply roughness correlation to CWT and MODWT coefficents
nested.function <- function(x) {
  
  # if x is a list of vectors, convert to matrix
  if(is.list(x) == TRUE){ 
    
    x <- do.call(rbind, x) 
    
    } 
  
  # apply roughness correlation to each wavelength
  stats <- t(apply(x, 1, roughness.correlation, ps = res)) 
  
}

# apply roughness correlation to CWT coefficients
cwt.roughness.correlation <- lapply(cwt.coef, nested.function) 

# apply roughness correlation to MODWT coefficients
dwt.roughness.correlation <- lapply(dwt.coef, nested.function) 











































