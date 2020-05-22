# trc_plot.R
#
# produces plots from trc_code.R
#
# 2020-FEB DLA
#
################################################################################################################################


library(magicaxis) # for plot axes
library(RColorBrewer) # for plot colour scheme

# define maximum grain size from Experiment 1 for use in the plots
dmax = 0.0064

# define point in the experiment that each elevation profile corresponds to
time.int = c(0, 5, 5, 5, 5, 10, 10, 10, 10, 15, 15, 15, 15, 30, 30, 30, 30, 60, 60, 60, 60, 120, 120, 120, 120)
time.cum = cumsum(time.int) # cumulative time

# let's use the DWT results, rather than the CWT results
roughness = dwt.roughness.correlation
scale = dwt.scale

n = length(roughness) # determine necessarily length of colour ramp
getPalette = colorRampPalette(brewer.pal(11, "RdYlBu")) # red to yellow to blue colour scheme
cl = getPalette(n) # use getPalette function to interpolate colours

# define powers-of-ten rounding functions for log plots
roundup <- function(x) 10^ceiling(log10(x))
rounddown <- function(x) 10^floor(log10(x))

# import flow resistance data from other experiments
exp1.fr <- read.csv("str_exp1_summary.txt")
exp2.fr <- read.csv("str_exp2_summary.txt")
exp3.fr <- read.csv("str_exp3_summary.txt")
exp4.fr <- read.csv("nar_exp1_summary.txt")
exp5.fr <- read.csv("nar_exp3_summary.txt")
sp.fr <- data.frame(read.csv("hohermuth2018_summary.txt")) # laboraory step-pool experiments from Hohermuth and Weitbrecht (2018)

# combine original experiments conducted in the A-BES together
pr.fr <- data.frame(rbind(exp1.fr, exp2.fr, exp3.fr, exp4.fr, exp5.fr))


################################################################################################################################


# FORM SIZE DISTRIBUTION (FSD)

png("fsd.png", width = 15, height = 15, units = 'cm', res = 600)

#find y axis limits
ymax = numeric(n)
ymin = numeric(n)
for (i in 1:n){
  
  ymax[i] <- max(roughness[[i]][, 1])
  ymin[i] <- min(roughness[[i]][, 1])
}

# make initial plot, without axes
plot(scale[[1]], roughness[[1]][, 1],
     type = "l",
     log = "xy",
     ylim = c(rounddown(min(ymin)), roundup(max(ymax))),
     col = cl[1],
     xaxt = 'n',
     yaxt = 'n',
     ann = FALSE)

# add all other lines
for(i in 2:n){
  
  lines(scale[[i]], roughness[[i]][, 1], 
        type = "l",
        col = cl[i])
  
}

# plot line representing dmax
abline(v = dmax, lty = 2)

# draw axes
magaxis(xlab = "Wavelength (m)",
        ylab = expression(sigma[z]))

# define legend
legend("bottomright", bty = "n", inset = c(0.05, 0.05), ncol = 1, lty = 1, 
       col = c(cl[1], cl[length(cl)/2], cl[length(cl)]), 
       legend = c(paste(time.cum[1], "mins"), paste(time.cum[length(time.cum)/2], "mins"), paste(time.cum[n], "mins")))

legend("topleft", expression(bold("(a)")), bty="n")


dev.off()


# cumulative FORM SIZE DISTRIBUTION (cFSD)

png("cfsd.png", width = 15, height = 15, units = 'cm', res = 600)

# for cumulative plot, let's limit the number of wavelengths to the minimum across the dataset
wavelengths = numeric(n) 
for (i in 1:n){
  
  wavelengths[i] <- length(roughness[[i]][, 1]) 
  
}

max.wavelength <- min(wavelengths)

# define cumsum function
cum.sum.perc <- function(x){
  
  cs <- cumsum(x)
  csp <- (cs / max(cs)) * 100
  
  return(csp)
  
}

plot(scale[[1]][1:max.wavelength], cum.sum.perc(roughness[[1]][1:max.wavelength, 1]),
     type = "l",
     log = "x",
     ylim = c(0, 100),
     col = cl[1],
     xaxt = 'n',
     yaxt = 'n',
     ann = FALSE)

for(i in 2:n){
  
  lines(scale[[i]][1:max.wavelength], cum.sum.perc(roughness[[i]][1:max.wavelength, 1]), 
        type = "l",
        col = cl[i])
  
}

abline(v = dmax, lty = 2)

magaxis(xlab = "Wavelength (m)",
        ylab = expression("Cumulative %" ~  sigma[z]))

legend("bottomright", bty = "n", inset = c(0.05, 0.05), ncol = 1, lty = 1, 
       col = c(cl[1], cl[length(cl)/2], cl[length(cl)]), 
       legend = c(paste(time.cum[1], "mins"), paste(time.cum[length(time.cum)/2], "mins"), paste(time.cum[n], "mins")))

legend("topleft", expression(bold("(b)")), bty="n")

dev.off()



# DRAG SIZE DISTRIBUTION (DSD)

png("dsd.png", width = 15, height = 15, units = 'cm', res = 600)

ymax = numeric(n)
ymin = numeric(n)
for (i in 1:n){
  
  ymax[i] <- max(roughness[[i]][, 5])
  ymin[i] <- min(roughness[[i]][, 5])
}

plot(scale[[1]], roughness[[1]][, 5],
     type = "l",
     log = "xy",
     ylim = c(rounddown(min(ymin)), roundup(max(ymax))),
     col = cl[1],
     xaxt = 'n',
     yaxt = 'n',
     ann = FALSE)

for(i in 2:n){
  
  lines(scale[[i]], roughness[[i]][, 5], 
        type = "l",
        col = cl[i])
  
}

abline(v = dmax, lty = 2)

magaxis(xlab = "Wavelength (m)",
        ylab = expression(k["s,pred"]))

legend("bottom", bty = "n", inset = c(0.05, 0.05), ncol = 1, lty = 1, 
       col = c(cl[1], cl[length(cl)/2], cl[length(cl)]), 
       legend = c(paste(time.cum[1], "mins"), paste(time.cum[length(time.cum)/2], "mins"), paste(time.cum[n], "mins")))

legend("topleft", expression(bold("(a)")), bty="n")

dev.off()



# cumulative DRAG SIZE DISTRIBUTION (cDSD)

png("cdsd.png", width = 15, height = 15, units = 'cm', res = 600)

plot(scale[[1]][1:max.wavelength], cum.sum.perc(roughness[[1]][1:max.wavelength, 5]),
     type = "l",
     log = "x", 
     ylim = c(0, 100),
     col = cl[1],
     xaxt = 'n',
     yaxt = 'n',
     ann = FALSE)

for(i in 2:n){
  
  lines(scale[[i]][1:max.wavelength], cum.sum.perc(roughness[[i]][1:max.wavelength, 5]), 
        type = "l",
        col = cl[i])
  
}

abline(v = dmax, lty = 2)

magaxis(xlab = "Wavelength (m)",
        ylab = expression("Cumulative %" ~ k["s,pred"]))

legend("bottomright", bty = "n", inset = c(0.05, 0.05), ncol = 1, lty = 1, 
       col = c(cl[1], cl[length(cl)/2], cl[length(cl)]), 
       legend = c(paste(time.cum[1], "mins"), paste(time.cum[length(time.cum)/2], "mins"), paste(time.cum[n], "mins")))

legend("topleft", expression(bold("(b)")), bty="n")

dev.off()


# EFFECTIVE SLOPE

png("es.png", width = 15, height = 15, units = 'cm', res = 600)

ymax = numeric(n)
ymin = numeric(n)
for (i in 1:n){
  
  ymax[i] <- max(roughness[[i]][, 4])
  ymin[i] <- min(roughness[[i]][, 4])
}

plot(scale[[1]], roughness[[1]][, 4],
     type = "l",
     log = "xy",
     ylim = c(rounddown(min(ymin)), roundup(max(ymax))),
     col = cl[1],
     xaxt = 'n',
     yaxt = 'n',
     ann = FALSE)

# plot range of ES values used to develop Forooghi's (2017) roughness correlation
es.high <- rep(0.89, times = length(scale[[1]]))
es.low <- rep(0.07, times = length(scale[[1]]))

polygon(c(scale[[1]], rev(scale[[1]])), 
        c(es.high, rev(es.low)),
        col = "grey90", border = NA)

for(i in 1:n){
  
  lines(scale[[i]], roughness[[i]][, 4], 
        type = "l",
        col = cl[i])
  
}

abline(v = dmax, lty = 2)

magaxis(xlab = "Wavelength (m)",
        ylab = "ES")

legend("bottom", bty = "n", inset = c(0.05, 0.05), ncol = 1, lty = 1, 
       col = c(cl[1], cl[length(cl)/2], cl[length(cl)]), 
       legend = c(paste(time.cum[1], "mins"), paste(time.cum[length(time.cum)/2], "mins"), paste(time.cum[n], "mins")))

dev.off()



# ks/k

png("ksk.png", width = 15, height = 15, units = 'cm', res = 600)

ymax = numeric(n)
ymin = numeric(n)
for (i in 1:n){
  
  ymax[i] <- max(roughness[[i]][, 6])
  ymin[i] <- min(roughness[[i]][, 6])
}

plot(scale[[1]], roughness[[1]][, 6],
     type = "l",
     log = "xy",
     ylim = c(rounddown(min(ymin)), roundup(11.78)),
     col = cl[1],
     xaxt = 'n',
     yaxt = 'n',
     ann = FALSE)

# plot range of ks/k values used to develop Forooghi's (2017) roughness correlation
ksk.high <- rep(11.78, times = length(scale[[1]]))
ksk.low <- rep(1.67, times = length(scale[[1]]))

polygon(c(scale[[1]], rev(scale[[1]])), 
        c(ksk.high, rev(ksk.low)),
        col = "grey90", border = NA)

for(i in 1:n){
  
  lines(scale[[i]], roughness[[i]][, 6], 
        type = "l",
        col = cl[i])
  
}

abline(v = dmax, lty = 2)

magaxis(xlab = "Wavelength (m)",
        ylab = expression(k[s] / k))

legend("bottom", bty = "n", inset = c(0.05, 0.05), ncol = 1, lty = 1, 
       col = c(cl[1], cl[length(cl)/2], cl[length(cl)]), 
       legend = c(paste(time.cum[1], "mins"), paste(time.cum[length(time.cum)/2], "mins"), paste(time.cum[n], "mins")))

dev.off()


################################################################################################################################


# SENSITIVITY ANALYSIS OF KSpred USING MOTHER WAVELET AND K VALUE

png("sensitivity.png", width = 15, height = 15, units = 'cm', res = 600)

sensitivity <- read.csv("sensitivity.txt") #import data from sensitivity analysis

ymax <- as.numeric(max(sensitivity))
ymin <- as.numeric(min(sensitivity))

plot(time.cum[2:n], sensitivity[, 1], 
     type = "o",
     pch = 15,
     ylim = c(ymin, ymax), 
     log = "y",
     xaxt = 'n',
     yaxt = 'n',
     ann = FALSE)

magaxis(xlab = "Time (min)",
        ylab = expression(Sigma ~ k["s,pred"]))

points(time.cum[2:n], sensitivity[, 2], pch = 16, type = "o")
points(time.cum[2:n], sensitivity[, 3], pch = 17, type = "o")
points(time.cum[2:n], sensitivity[, 4], pch = 1, type = "o")
points(time.cum[2:n], sensitivity[, 5], pch = 2, type = "o")

items <- c(expression(k ~ "=" ~ sigma[z] ~ "," ~  Psi ~ "=" ~ "db4"),
           expression(k ~ "=" ~ sigma[z] ~ "," ~  Psi ~ "=" ~ "db8"),
           expression(k ~ "=" ~ sigma[z] ~ "," ~  Psi ~ "=" ~ "db16"),
           expression(k ~ "=" ~ 2*sigma[z] ~ "," ~  Psi ~ "=" ~ "db4"),
           expression(k ~ "=" ~ z[max] - z[min] ~ "," ~  Psi ~ "=" ~ "db4"))
                      
legend("topright", bty = "n", inset = c(0.05, 0.05), ncol = 1,
       lty = 1, pch = c(15, 16, 17, 1, 2), 
       legend = items)

dev.off()


# plot two wavelengths of interest from the MPDWT and CWT transforms

png("comparison.png", width = 15, height = 10, units = 'cm', res = 600)

woi <- c(0.004, 1.024) # wavelengths of interest
doi <- 25 # MODWT of interst

#find wavelengths corresponding to these DEM wavelengths
dwt.one <- which(dwt.scale[[doi]] == woi[1])
dwt.two <- which(dwt.scale[[doi]] == woi[2])

cwt.one <- which(cwt.scale[[doi]] == woi[1])
cwt.two <- which(cwt.scale[[doi]] == woi[2])

recon <- colSums(do.call(rbind, dwt.coef[[doi]]))

xvals <- seq(from = res, to = (length(elev[[doi]])*res), by = res)

par(mfrow=c(4,1), mar = c(0.5,0.5,0.5,0.5))

par(oma=c(1,1,0,0)) # all sides have 3 lines of space



plot(xvals, elev[[doi]]*100, 
     type = "l",
     xaxt = 'n',
     yaxt = 'n',
     ann = FALSE)

magaxis(side = 2, xlab = "", minorn = FALSE, majorn = 4)

legend("bottomleft", expression(bold("(a) Original")), bty="n")



plot(xvals, cwt.coef[[doi]][cwt.two, ],
     type = "l",
     xaxt = 'n',
     yaxt = 'n',
     ann = FALSE)
lines(xvals, cwt.coef[[doi]][cwt.one, ])

legend("bottomleft", expression(bold("(b) CWT")), bty="n")


plot(xvals, dwt.coef[[doi]][[dwt.two]],
     type = "l",
     xaxt = 'n',
     yaxt = 'n',
     ann = FALSE)
lines(xvals, dwt.coef[[doi]][[dwt.one]])

legend("bottomleft", expression(bold("(c) MODWT")), bty="n")

plot(xvals, recon*100,
     type = "l",
     xaxt = 'n',
     yaxt = 'n',
     ann = FALSE)

legend("bottomleft", expression(bold("(d) MODWT (Recon.)")), bty="n")

magaxis(side = 2, xlab = "", minorn = FALSE, majorn = 4)

magaxis(side = 1, xlab = "", minorn = FALSE, majorn = 5)

dev.off()



# SUM Kspred vs Ks*pred

png("ksvssumks.png", width = 15, height = 15, units = 'cm', res = 600)

pr.ks <- log(pr.fr$ks)
pr.sumks <- log(pr.fr$sumks)
sp.ks <- log(sp.fr$ks)
sp.sumks <- log(sp.fr$sumks)

xmin <- min(c(pr.ks, sp.ks))
xmax <- max(c(pr.ks, sp.ks))
ymin <- min(c(pr.sumks, sp.sumks))
ymax <- max(c(pr.sumks, sp.sumks))

plot(pr.ks, pr.sumks, 
     pch = 16,
     xlim = c(xmin, xmax),
     ylim = c(ymin, ymax),
     xaxt = 'n',
     yaxt = 'n',
     ann = FALSE)

points(sp.ks, sp.sumks,
       pch = 17)

abline(0,1)

pr.reg = lm(pr.sumks ~ pr.ks)
abline(pr.reg, untf = T, lty = 3)

sp.reg = lm(sp.sumks ~ sp.ks)
abline(sp.reg, untf = T, lty = 5)

magaxis(xlab = expression("log("*k["s,pred"]^"*"*")"),
        ylab = expression("log("*Sigma ~ k["s,pred"]*")"))



legend("topleft", bty = "n", inset = c(0.05, 0.05), ncol = 1, 
       pch = c(26, 16, 17),
       lty = c(1, 3, 5),
       legend = c("1:1",
                  expression("lm PR (n = 164)"),
                  expression("lm SP (n = 83)")))

dev.off()

# FLOW RESISTANCE VS RELATIVE ROUGHNESS

png("fvsk.png", width = 15, height = 15, units = 'cm', res = 600)

# extract variables of interest, for simplicity
pr.x1 = pr.fr$h/pr.fr$k
pr.x2 = pr.fr$h/pr.fr$ks
pr.y = pr.fr$UUstar

sp.x1 = sp.fr$h/sp.fr$k
sp.x2 = sp.fr$h/sp.fr$ks
sp.y = sp.fr$UUstar

par(mar=c(0.5, 0.5, 0.2, 0.2), mfrow=c(1,2),
    oma = c(3, 3, 0.5, 0.5))

#bottom left top right

# pool riffle k = ks
plot(pr.x1, pr.y,
     log = "xy",
     ylim = c(1, 10),
     xlim = c(1, 100),
     pch = 16,
     xaxt = 'n',
     yaxt = 'n',
     ann = FALSE)

# step pool k = ks
points(sp.x1, sp.y,
       pch = 17)

# PLOT VPE EQUATIONS
x.seq <- seq(from = 0, to = 100, by = 0.1)

#fits from chen meta-analysis
c1 = 5.77
c2 = 1.24 

vpe.pred <- (c1 * c2 * (x.seq)) / (c1^2 + c2^2 * (x.seq)^(5/3))^(1/2)

lines(x.seq, vpe.pred, lty = 1)

magaxis(xlab = expression(h/sigma[z]),
        ylab = expression("U/U* =" ~ (8/f)^{1/2}))

legend("topleft", expression(bold("(a)")), bty="n")

legend("bottomright", bty = "n", inset = c(0.05, 0.05), ncol = 1,
       pch = c(16, 17),
       legend = c("PR", "SP"))

# pool riffle k = sd
plot(pr.x2, pr.y,
     log = "xy",
     ylim = c(1, 10),
     xlim = c(1, 100),
     pch = 16,
     xaxt = 'n',
     yaxt = 'n',
     ann = FALSE)

# step pool k = sd
points(sp.x2, sp.y,
       pch = 17)

lines(x.seq, vpe.pred, lty = 1)

#fits 
c1 = 4.82
c2 = 2.17

vpe.pred <- (c1 * c2 * (x.seq)) / (c1^2 + c2^2 * (x.seq)^(5/3))^(1/2)

lines(x.seq, vpe.pred, lty = 2)

magaxis(side = 1, xlab = expression(h/k["s,pred"]^"*"))

legend("bottomright", bty = "n", inset = c(0.05, 0.05), ncol = 1,
       lty = c(1, 2),
       legend = c("VPE (Chen)", "VPE (fitted)"))

legend("topleft", expression(bold("(b)")), bty="n")
       
dev.off()
