#!/usr/bin/env Rscript

# close all devices
graphics.off()

suppressMessages(library(GA))
source("fit_functions.R", local=T)

# set max number of GA iterations and population size
maxiter <- 4
popSize <- 10

### FUNCTIONS #########################################

plot.initial <- function() {
  par(mar=c(6, 15, 2, 2))
  d$phase <- (d$mjd - init.t0) %% init.p / init.p
  
  xlab <- "Orbital phase"
  ylab <- "PD [%]"
  
  
  plot.t <- seq(min(d$mjd), min(d$mjd) + init.p, 0.1)
  plot.pd <- fit.pd.e(plot.t, init.t0, init.p, init.i, init.e, init.o, init.p.max)
  plot.df <- data.frame(t=plot.t, pd=plot.pd)
  plot.df$phase <- (plot.df$t - init.t0) %% init.p / init.p
  
  plot.df <- plot.df[order(plot.df$phase), ]
  
  ylim <- c(min(d$pd - d$dpd, plot.df$pd), max(d$pd + d$dpd, plot.df$pd))
  
  with(d, plot(phase, pd, xlim=c(0, 1), ylim=ylim, xlab=xlab, ylab=ylab, pch=19, main="Initial parameters plot"))
  with(d, segments(x0=phase, y0=pd - dpd, y1=pd + dpd))
  with(plot.df, lines(phase, pd))
  mtext(text=paste("t0 =", t0), side=2, las=1, padj=0, adj=0)
}

pdf(file="orb_fit.pdf", title="orb fit", width=9, height=5)

### READ DATA #########################################

d <- read.csv("v0332.csv")
d$mjd <- d$jd - 2400000.5

d <- d[, c("mjd", "pd", "dpd")]

### SET INITIAL PARAMETERS ############################

init.t0 <- 57150.068
init.p <- 33.83
init.i <- 10.3
init.e <- 0.37
init.o <- 280
init.p.max <- 3.84

lower <- list(t0=init.t0 - 200, p=init.p - 10, i=0,  e=0.0, o=0, p.max=3.5)
upper <- list(t0=init.t0 + 200, p=init.p + 10, i=90, e=0.95, o=360, p.max=4.0)


### PLOT INITIAL CURVE ################################

plot.initial()

### FIT TEST DATA #####################################
test = F
if (test) {
  n <- 32
  mjd <- sort(runif(n, 55000, 55300))
  t0 <- 57234
  p <- 26.78
  i <- 38.1
  e <- 0.12
  o <- 45.9
  p.max <- 2.8
  dpd <- abs(rnorm(m=0.0, sd=0.02, n=n))
  pd <- fit.pd.e(mjd, t0, p, i, e, o, p.max) + dpd
  
  phase <- (mjd - t0) %% p / p
  
  plot(phase, pd)
  segments(x0=phase, y0=pd + dpd, y1=pd - dpd)

  lower <- list(t0=t0 - 101.16, p=p - 3, i=0, e=0, o=0, p.max=1.0)
  upper <- list(t0=t0 + 101.16, p=p + 3, i=90, e=1, o=180, p.max=5.0)

  test.df <- data.frame(mjd=mjd, pd=pd, dpd=rep(0.1, n))

  fit <- nls(pd ~ fit.pd.e(mjd, t0, p, i, e, o, p.max), 
           start=list(t0=t0, p=p, i=i, e=e, o=o, p.max=p.max), 
           lower=lower, 
           upper=upper, 
           algorithm="port")
  fit <- coef(fit)
  cat("parameter (initial - fitted) -> difference\n")
  cat(paste("t0", t0, round(fit[[1]], 2), round(t0 - fit[[1]], 2), "\n"))
  cat(paste("p", p, round(fit[[2]], 2), round(p - fit[[2]], 2), "\n"))
  cat(paste("i", i, round(fit[[3]], 2), round(i - fit[[3]], 2), "\n"))
  cat(paste("e", e, round(fit[[4]], 2), round(e - fit[[4]], 2), "\n"))
  cat(paste("o", o, round(fit[[5]], 2), round(o - fit[[5]], 2), "\n"))
  cat(paste("p.max", p.max, round(fit[[6]], 2), round(p.max - fit[[6]], 2), "\n"))
}

### FIT ###############################################

x <- d$mjd
y <- d$pd

fit <- nls(y ~ fit.pd.e(x, t0, p, i, e, o, p.max), 
           start=list(t0=init.t0, p=init.p, i=init.i, e=init.e, o=init.o, p.max=init.p.max), 
           lower=lower, 
           upper=upper, 
           algorithm="port")


dev.off()






