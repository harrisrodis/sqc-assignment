###################
##  Preparation  ##
###################

# load bias correlation constants
bcc <- read.csv("bcc.csv")

# set seed for replicating simulated results
set.seed(10001)

#####################################
##  1. Control Chart for the Mean  ##
#####################################

# setup function that creates the Xbar control chart for a given type I error
cc.Xbar <- function(dt, a) {
  
  # function to find the numeric value of a range
  range.sn <- function(x) { max(x) - min(x) }
  # find k based on type I error probability
  k <- abs(qnorm(a/2))
  # number of observations in each subgroup (sample size)
  n <- nrow(dt)
  # determine the value of d2 constant based on the sample size
  d2 <- bcc[n - 1, 2]
  
  # function to calculate the central line, lower and upper limits
  xbar.lims <- function(means, ranges) {
    # calculate the mean range of all subgroups
    rbar <- mean(ranges, na.rm = T)
    # find the central line of the chart (mean of means)
    cl <- mean(xbars, na.rm = T)
    # calculate the upper control limit
    ucl <- cl + (k/(d2 * sqrt(nrow(dt)))) * rbar
    # calculate the lower control limit
    lcl <- cl - (k/(d2 * sqrt(nrow(dt)))) * rbar
    # group and return them
    return(list(cl, lcl, ucl))
  }
  
  # ----------------------------------------------------------------------------- #
  
  # Phase 1
  
  # number of subgroups in phase 1
  m <- ncol(dt) / 2
  # calculate the mean of each subgroup
  xbars <- apply(dt[, 1:m], 2, mean)
  # calculate the range of each subgroup
  Rs <- apply(dt[, 1:m], 2, range.sn)
  # calculate initial limits
  lims <- xbar.lims(xbars, Rs)
  # check and correct for outliers
  while(any(xbars < lims[[2]] | xbars > lims[[3]], na.rm = T)) {
    xbars[which(xbars < lims[[2]] | xbars > lims[[3]])] <- NA
    Rs[which(xbars < lims[[2]] | xbars > lims[[3]])] <- NA
    lims <- xbar.lims(xbars, Rs)
  }
  
  # ----------------------------------------------------------------------------- #
  
  # Phase 2
  
  xbars.all <- c(na.omit(xbars), apply(dt[, (m + 1):ncol(dt) ], 2, mean))
  
  # ----------------------------------------------------------------------------- #
  
  # generate chart
  par(
    family = "Source Sans Pro",
    cex = 0.8,
    las = 1,
    cex.lab = 1.2,
    col.axis = "#333333",
    fg = "#333333"
  )
  plot(
    xbars.all, 
    type = "b", 
    pch = 16, 
    lwd = 2, 
    col = "#f6b26b",
    frame = F,
    ylim = c(2 * lims[[2]] - lims[[1]], 2 * lims[[3]] - lims[[1]]),
    xlab = "Δειγματοληψία",
    ylab = "Μέση τιμή"
  )
  abline(h = lims[[1]], lwd = 2, col = "#9fc5e8")
  abline(h = lims[[2]], lwd = 2, lty = 2, col = "#e06666")
  abline(h = lims[[3]], lwd = 2, lty = 2, col = "#e06666")
  abline(v = length(xbars[which(!is.na(xbars))]), lwd = 1, lty = 3)
  text(
    60, 
    lims[[3]] + (lims[[3]] * 0.01),
    pos = 2,
    labels = paste0("UCL: ", round(lims[[3]],3)), 
    col = "#e06666",
    font = 2
  )
  text(
    60, 
    lims[[2]] - (lims[[2]] * 0.01),
    pos = 2,
    labels = paste0("LCL: ", round(lims[[2]],3)), 
    col = "#e06666",
    font = 2
  )
  text(
    length(xbars[which(!is.na(xbars))]) - 0.2, 
    lims[[3]] + (lims[[3]] * 0.01),
    pos = 2,
    labels = "Τέλος Φάσης 1",
    font = 2,
    cex = 0.9
  )
}

# generate random data for scenario 1.1 (under control)
data.1.1 <- matrix(rep(rnorm(600, 1, .05)), ncol = 60, nrow = 10)
cc.Xbar(data.1.1, .0027)

# generate random data for scenario 1.2 (out of control)
data.1.2 <- cbind(
  matrix(rep(rnorm(300, 1, .05)), ncol = 30, nrow = 10),
  matrix(rep(rnorm(300, 1.05, .05)), ncol = 30, nrow = 10)
)
cc.Xbar(data.1.2, .0027)

######################################
##  2. Control Chart for the Range  ##
######################################

# setup function that creates the R control chart for a given type I error
cc.R <- function(dt, a) {
  
  # function to find the numeric value of a range
  range.sn <- function(x) { max(x) - min(x) }
  # find k based on type I error probability
  k <- abs(qnorm(a/2))
  # number of observations in each subgroup (sample size)
  n <- nrow(dt)
  # determine the value of d2 & d3 constant based on the sample size
  d2 <- bcc[n - 1, 2]
  d3 <- bcc[n - 1, 3]
  
  # function to calculate the central line, lower and upper limits
  rchart.lims <- function(ranges) {
    # find the central line of the chart (mean range of all subgroups)
    cl <- mean(ranges, na.rm = T)
    # calculate the upper control limit
    ucl <- cl + ((k * d3) / d2) * cl
    # calculate the lower control limit
    lcl <- cl - ((k * d3) / d2) * cl
    if (lcl < 0) { lcl <- 0 }
    # group and return them
    return(list(cl, lcl, ucl))
  }
  
  # ----------------------------------------------------------------------------- #
  
  # Phase 1
  
  # number of subgroups in phase 1
  m <- ncol(dt) / 2
  # calculate the range of each subgroup
  Rs <- apply(dt[, 1:m], 2, range.sn)
  # calculate initial limits
  lims <- rchart.lims(Rs)
  # check and correct for outliers
  while(any(Rs < lims[[2]] | Rs > lims[[3]], na.rm = T)) {
    Rs[which(Rs < lims[[2]] | Rs > lims[[3]])] <- NA
    lims <- rchart.lims(Rs)
  }
  
  # ----------------------------------------------------------------------------- #
  
  # Phase 2
  
  Rs.all <- c(na.omit(Rs), apply(dt[, (m + 1):ncol(dt) ], 2, range.sn))
  
  # ----------------------------------------------------------------------------- #
  
  # generate chart
  par(
    family = "Source Sans Pro",
    cex = 0.8,
    las = 1,
    cex.lab = 1.2,
    col.axis = "#333333",
    fg = "#333333"
  )
  plot(
    Rs.all, 
    type = "b", 
    pch = 16, 
    lwd = 2, 
    col = "#f6b26b",
    frame = F,
    ylim = c(2 * lims[[2]] - lims[[1]], 2 * lims[[3]] - lims[[1]]),
    xlab = "Δειγματοληψία",
    ylab = "Μέσο εύρος"
  )
  abline(h = lims[[1]], lwd = 2, col = "#9fc5e8")
  abline(h = lims[[2]], lwd = 2, lty = 2, col = "#e06666")
  abline(h = lims[[3]], lwd = 2, lty = 2, col = "#e06666")
  abline(v = length(Rs[which(!is.na(Rs))]), lwd = 1, lty = 3)
  text(
    60, 
    lims[[3]] + (lims[[3]] * 0.05),
    pos = 2,
    labels = paste0("UCL: ", round(lims[[3]],3)), 
    col = "#e06666",
    font = 2
  )
  text(
    60, 
    lims[[2]] - (lims[[2]] * 0.4),
    pos = 2,
    labels = paste0("LCL: ", round(lims[[2]],3)), 
    col = "#e06666",
    font = 2
  )
  text(
    length(Rs[which(!is.na(Rs))]) - 0.2, 
    lims[[3]] + (lims[[3]] * 0.05),
    pos = 2,
    labels = "Τέλος Φάσης 1",
    font = 2,
    cex = 0.9
  )
}

# generate random data for scenario 2.1 (under control)
data.2.1 <- matrix(rep(rnorm(600, 1, .05)), ncol = 60, nrow = 10)
cc.R(data.2.1, .0027)

# generate random data for scenario 2.2 (out of control)
data.2.2 <- cbind(
  matrix(rep(rnorm(300, 1, .05)), ncol = 30, nrow = 10),
  matrix(rep(rnorm(300, 1, .08)), ncol = 30, nrow = 10)
)
cc.R(data.2.2, .0027)

###################################################
##  3. Control Chart for the Standard Deviation  ##
###################################################

# setup function that creates the S control chart for a given type I error
cc.S <- function(dt, a) {

  # find k based on type I error probability
  k <- abs(qnorm(a/2))
  # number of observations in each subgroup (sample size)
  n <- nrow(dt)
  # determine the value of c4 constant based on the sample size
  c4 <- bcc[n - 1, 4]
  
  # function to calculate the central line, lower and upper limits
  schart.lims <- function(errors) {
    # find the central line of the chart (mean range of all subgroups)
    cl <- mean(errors, na.rm = T)
    # calculate the upper control limit
    ucl <- cl + k * (cl / c4) * sqrt(1 - c4^2)
    # calculate the lower control limit
    lcl <- cl - k * (cl / c4) * sqrt(1 - c4^2)
    if (lcl < 0) { lcl <- 0 }
    # group and return them
    return(list(cl, lcl, ucl))
  }
  
  # ----------------------------------------------------------------------------- #
  
  # Phase 1
  
  # number of subgroups in phase 1
  m <- ncol(dt) / 2
  # calculate the standard deviation of each subgroup
  SDs <- apply(dt[, 1:m], 2, sd)
  # calculate initial limits
  lims <- schart.lims(SDs)
  # check and correct for outliers
  while(any(SDs < lims[[2]] | SDs > lims[[3]], na.rm = T)) {
    SDs[which(SDs < lims[[2]] | SDs > lims[[3]])] <- NA
    lims <- schart.lims(SDs)
  }
  
  # ----------------------------------------------------------------------------- #
  
  # Phase 2
  
  SDs.all <- c(na.omit(SDs), apply(dt[, (m + 1):ncol(dt) ], 2, sd))
  
  # ----------------------------------------------------------------------------- #
  
  # generate chart
  par(
    family = "Source Sans Pro",
    cex = 0.8,
    las = 1,
    cex.lab = 1.2,
    col.axis = "#333333",
    fg = "#333333"
  )
  plot(
    SDs.all, 
    type = "b", 
    pch = 16, 
    lwd = 2, 
    col = "#f6b26b",
    frame = F,
    ylim = c(2 * lims[[2]] - lims[[1]], 2 * lims[[3]] - lims[[1]]),
    xlab = "Δειγματοληψία",
    ylab = "Τυπική απόκλιση"
  )
  abline(h = lims[[1]], lwd = 2, col = "#9fc5e8")
  abline(h = lims[[2]], lwd = 2, lty = 2, col = "#e06666")
  abline(h = lims[[3]], lwd = 2, lty = 2, col = "#e06666")
  abline(v = length(SDs[which(!is.na(SDs))]), lwd = 1, lty = 3)
  text(
    60, 
    lims[[3]] + (lims[[3]] * 0.05),
    pos = 2,
    labels = paste0("UCL: ", round(lims[[3]],3)), 
    col = "#e06666",
    font = 2
  )
  text(
    60, 
    lims[[2]] - (lims[[2]] * 0.4),
    pos = 2,
    labels = paste0("LCL: ", round(lims[[2]],3)), 
    col = "#e06666",
    font = 2
  )
  text(
    length(SDs[which(!is.na(SDs))]) - 0.2, 
    lims[[3]] + (lims[[3]] * 0.05),
    pos = 2,
    labels = "Τέλος Φάσης 1",
    font = 2,
    cex = 0.9
  )
}

# generate random data for scenario 3.1 (under control)
data.3.1 <- matrix(rep(rnorm(1200, 1, .05)), ncol = 60, nrow = 20)
cc.S(data.3.1, .0027)

# generate random data for scenario 3.2 (out of control)
data.3.2 <- cbind(
  matrix(rep(rnorm(600, 1, .05)), ncol = 30, nrow = 20),
  matrix(rep(rnorm(600, 1, .065)), ncol = 30, nrow = 20)
)
cc.S(data.3.2, .0027)

#########################################################
##  4. Control Chart for the percentage of defectives  ##
#########################################################

# setup function that creates the P control chart for a given type I error
cc.P <- function(dt, a) {
  
  # find k based on type I error probability
  k <- abs(qnorm(a/2))
  # number of observations in each subgroup (sample size)
  n <- nrow(dt)
  
  # function to calculate the central line, lower and upper limits
  pchart.lims <- function(defectives) {
    # find the central line of the chart (probability of defective)
    cl <- mean(defectives, na.rm = T)
    # calculate the upper control limit
    ucl <- cl + k * sqrt((cl * (1 - cl)) / n)
    # calculate the lower control limit
    lcl <- cl - k * sqrt((cl * (1 - cl)) / n)
    if (lcl < 0) { lcl <- 0 }
    # group and return them
    return(list(cl, lcl, ucl))
  }
  
  # ----------------------------------------------------------------------------- #
  
  # Phase 1
  
  # number of subgroups in phase 1
  m <- ncol(dt) / 2
  # calculate the number of defectives of each subgroup
  defs <- apply(dt[, 1:m], 2, sum) / n
  # calculate initial limits
  lims <- pchart.lims(defs)
  # check and correct for outliers
  while(any(defs < lims[[2]] | defs > lims[[3]], na.rm = T)) {
    defs[which(defs < lims[[2]] | defs > lims[[3]])] <- NA
    lims <- pchart.lims(defs)
  }
  
  # ----------------------------------------------------------------------------- #
  
  # Phase 2
  
  defs.all <- c(na.omit(defs), apply(dt[, (m + 1):ncol(dt) ], 2, sum) / n)
  
  # ----------------------------------------------------------------------------- #
  
  # generate chart
  par(
    family = "Source Sans Pro",
    cex = 0.8,
    las = 1,
    cex.lab = 1.2,
    col.axis = "#333333",
    fg = "#333333"
  )
  plot(
    defs.all, 
    type = "b", 
    pch = 16, 
    lwd = 2, 
    col = "#f6b26b",
    frame = F,
    ylim = c(2 * lims[[2]] - lims[[1]], 2 * lims[[3]] - lims[[1]]),
    xlab = "Δειγματοληψία",
    ylab = "Ποσοστό ελαττωματικών"
  )
  abline(h = lims[[1]], lwd = 2, col = "#9fc5e8")
  abline(h = lims[[2]], lwd = 2, lty = 2, col = "#e06666")
  abline(h = lims[[3]], lwd = 2, lty = 2, col = "#e06666")
  abline(v = length(defs[which(!is.na(defs))]), lwd = 1, lty = 3)
  text(
    40, 
    lims[[3]] + (lims[[3]] * 0.1),
    pos = 2,
    labels = paste0("UCL: ", round(lims[[3]],3)), 
    col = "#e06666",
    font = 2
  )
  text(
    40, 
    lims[[2]] - 0.005,
    pos = 2,
    labels = paste0("LCL: ", round(lims[[2]],3)), 
    col = "#e06666",
    font = 2
  )
  text(
    length(defs[which(!is.na(defs))]) - 0.2, 
    lims[[3]] + (lims[[3]] * 0.1),
    pos = 2,
    labels = "Τέλος Φάσης 1",
    font = 2,
    cex = 0.9
  )
}

# generate random data for scenario 4.1 (under control)
data.4.1 <- matrix(rep(rbinom(10000, 1, .02)), ncol = 40, nrow = 250)
cc.P(data.4.1, .0027)

# generate random data for scenario 4.2 (out of control)
data.4.2 <- cbind(
  matrix(rep(rbinom(5000, 1, .02)), ncol = 20, nrow = 250),
  matrix(rep(rbinom(5000, 1, .04)), ncol = 20, nrow = 250)
)
cc.P(data.4.2, .0027)

##################################################
##  5. Control Chart for the number of defects  ##
##################################################

# setup function that creates the C control chart for a given type I error
cc.C <- function(dt, a) {
  
  # find k based on type I error probability
  k <- abs(qnorm(a/2))
  
  # function to calculate the central line, lower and upper limits
  cchart.lims <- function(defects) {
    # find the central line of the chart (defects per subgroup)
    cl <- mean(defects, na.rm = T)
    # calculate the upper control limit
    ucl <- cl + k * sqrt(cl)
    # calculate the lower control limit
    lcl <- cl - k * sqrt(cl)
    if (lcl < 0) { lcl <- 0 }
    # group and return them
    return(list(cl, lcl, ucl))
  }
  
  # ----------------------------------------------------------------------------- #
  
  # Phase 1

  # calculate the number of defects of each subgroup
  defs <- dt[, 1]
  # calculate initial limits
  lims <- cchart.lims(defs)
  # check and correct for outliers
  while(any(defs < lims[[2]] | defs > lims[[3]], na.rm = T)) {
    defs[which(defs < lims[[2]] | defs > lims[[3]])] <- NA
    lims <- cchart.lims(defs)
  }
  
  # ----------------------------------------------------------------------------- #
  
  # Phase 2
  
  defs.all <- c(na.omit(defs), dt[, 2])
  
  # ----------------------------------------------------------------------------- #
  
  # generate chart
  par(
    family = "Source Sans Pro",
    cex = 0.8,
    las = 1,
    cex.lab = 1.2,
    col.axis = "#333333",
    fg = "#333333"
  )
  plot(
    defs.all, 
    type = "b", 
    pch = 16, 
    lwd = 2, 
    col = "#f6b26b",
    frame = F,
    ylim = c(2 * lims[[2]] - lims[[1]], 2 * lims[[3]] - lims[[1]]),
    xlab = "Δειγματοληψία",
    ylab = "Αριθμός ελαττωμάτων"
  )
  abline(h = lims[[1]], lwd = 2, col = "#9fc5e8")
  abline(h = lims[[2]], lwd = 2, lty = 2, col = "#e06666")
  abline(h = lims[[3]], lwd = 2, lty = 2, col = "#e06666")
  abline(v = length(defs[which(!is.na(defs))]), lwd = 1, lty = 3)
  text(
    60, 
    lims[[3]] + 2,
    pos = 2,
    labels = paste0("UCL: ", round(lims[[3]],3)), 
    col = "#e06666",
    font = 2
  )
  text(
    60, 
    lims[[2]] - 2,
    pos = 2,
    labels = paste0("LCL: ", round(lims[[2]],3)), 
    col = "#e06666",
    font = 2
  )
  text(
    length(defs[which(!is.na(defs))]) - 0.2, 
    lims[[3]] + 2,
    pos = 2,
    labels = "Τέλος Φάσης 1",
    font = 2,
    cex = 0.9
  )
}

# generate random data for scenario 5.1 (under control)
data.5.1 <- matrix(rep(rpois(60, 10)), ncol = 2, nrow = 30)
cc.C(data.5.1, .0027)

# generate random data for scenario 5.2 (out of control)
data.5.2 <- cbind(
  matrix(rep(rpois(30, 10)), ncol = 1, nrow = 30),
  matrix(rep(rpois(30, 15)), ncol = 1, nrow = 30)
)
cc.C(data.5.2, .0027)

#####################
##  6. EWMA chart  ##
#####################

# setup function that creates the EWMA chart for a given type I error & lambda parameter
cc.EWMA <- function(dt, a, lambda) {
  
  # find k based on type I error probability
  k <- abs(qnorm(a/2))
  
  # function to calculate the central line, lower and upper limits
  ewma.lims <- function(obs, mu, sigma) {
    # estimate the central line of the chart (mean of observations)
    cl <- mu
    # calculate the upper control limits
    ucl <- vector(length = length(obs))
    for (i in 1:length(obs)) {
      ucl[i] <- cl + k * sigma * sqrt(lambda * (1 - (1 - lambda)^(2 * i)) / (2 - lambda)) 
    }
    # calculate the lower control limits
    lcl <- vector(length = length(obs))
    for (i in 1:length(obs)) {
      lcl[i] <- cl - k * sigma * sqrt(lambda * (1 - (1 - lambda)^(2 * i)) / (2 - lambda)) 
    }
    # group and return them
    return(list(cl, lcl, ucl))
  }
  
  # ----------------------------------------------------------------------------- #
  
  # Phase 1
  
  # Calculate the mean of phase 1 observations
  mu.0 <- mean(dt[1:(length(dt) / 2)])
  # Estimate the standard deviation of phase 1 observations
  s.hat <- sd(dt[1:(length(dt) / 2)])
  
  # ----------------------------------------------------------------------------- #
  
  # Phase 2
  
  # calculate the ewma statistic for each observation
  z <- mu.0
  ewma.obs <- vector(length = length(dt))
  for (i in 1:length(dt)) {
    ewma.obs[i] <- lambda * dt[i] + (1 - lambda) * z
    z <- ewma.obs[i]
  }
  # calculate ewma chart limits
  lims <- ewma.lims(ewma.obs, mu.0, s.hat)
  
  # ----------------------------------------------------------------------------- #
  
  # generate chart
  par(
    family = "Source Sans Pro",
    cex = 0.8,
    las = 1,
    cex.lab = 1.2,
    col.axis = "#333333",
    fg = "#333333"
  )
  plot(
    ewma.obs, 
    type = "b", 
    pch = 16, 
    lwd = 2, 
    col = "#f6b26b",
    frame = F,
    ylim = c(lims[[2]][length(ewma.obs)] - 0.05 * lims[[2]][length(ewma.obs)], 
             lims[[3]][length(ewma.obs)] + 0.05 * lims[[3]][length(ewma.obs)]),
    xlab = "Παρατήρηση",
    ylab = "Εκθετικά σταθμισμένη τιμή παρατήρησης"
  )
  abline(h = lims[[1]], lwd = 2, col = "#9fc5e8")
  lines(lims[[2]], lwd = 2, type = "s", col = "#e06666")
  lines(lims[[3]], lwd = 2, type = "s", col = "#e06666")
  text(
    length(ewma.obs), 
    lims[[3]][length(ewma.obs)] + 0.01 * lims[[3]][length(ewma.obs)],
    pos = 2,
    labels = paste0("UCL: ", round(lims[[3]][length(ewma.obs)],3)), 
    col = "#e06666",
    font = 2
  )
  text(
    length(ewma.obs), 
    lims[[2]][length(ewma.obs)] - 0.01 * lims[[2]][length(ewma.obs)],
    pos = 2,
    labels = paste0("LCL: ", round(lims[[2]][length(ewma.obs)],3)), 
    col = "#e06666",
    font = 2
  )
}

# generate random data for scenario 6.1 (under control)
data.6.1 <- rnorm(100, 1, .05)
cc.EWMA(data.6.1, .0027, .2)

# generate random data for scenario 6.2 (out of control)
data.6.2 <- c(rnorm(50, 1, .05), rnorm(50, 1.05, .05))
cc.EWMA(data.6.2, .0027, .2)

######################
##  7. CUSUM chart  ##
######################

# setup function that creates the CUSUM chart for a given minimum shift from target value
cc.CUSUM <- function(dt, minshift) {
  # calculate the recommended value of the decision interval
  H <- 5 * sd(dt)
  # calculate K for the minimum shift from target value
  K <- abs(minshift) / 2
  
  # ----------------------------------------------------------------------------- #
  
  # Phase 1

  # find the central line of the chart (mean of phase 1 observations)
  cl <- mean(dt[1:(length(dt) / 2)])
  # calculate the upper control limit
  ucl <- H
  # calculate the lower control limit
  lcl <- -H
  
  # ----------------------------------------------------------------------------- #
  
  # Phase 2
  
  # calculate positive cumulative sums of observations
  c <- 0
  c.plus <- vector(length = length(dt))
  for(i in 1:length(dt)) {
    c.plus[i] <- max(0, dt[i] - (cl + K) + c)
    c <- c.plus[i]
  }
  
  # calculate negative cumulative sums of observations
  c <- 0
  c.minus <- vector(length = length(dt))
  for(i in 1:length(dt)) {
    c.minus[i] <- max(0, (cl - K) - dt[i] + c)
    c <- c.minus[i]
  }
  c.minus <- -c.minus
  
  # ----------------------------------------------------------------------------- #
  
  # generate chart
  par(
    family = "Source Sans Pro",
    cex = 0.8,
    las = 1,
    cex.lab = 1.2,
    col.axis = "#333333",
    fg = "#333333"
  )
  plot(
    1, 
    type = "n",
    frame = F,
    ylim = c(min(lcl * 2, min(c.minus)), max(ucl * 2, max(c.plus))),
    xlim = c(0, length(dt)),
    xlab = "Παρατήρηση",
    ylab = "Συσσωρευμένο άθροισμα"
  )
  abline(h = 0, lwd = 2, col = "#9fc5e8")
  abline(h = lcl, lwd = 2, lty = 2, col = "#e06666")
  abline(h = ucl, lwd = 2, lty = 2, col = "#e06666")
  lines(c.minus, lwd = 2, type = "l", col = "#76a5af")
  points(seq_along(c.minus), c.minus, pch = 21, col = "#76a5af", bg = "#ffffff")
  lines(c.plus, lwd = 2, type = "l", col = "#f6b26b")
  points(seq_along(c.plus), c.plus, pch = 16, col = "#f6b26b")
  text(
    length(dt), 
    ucl * 1.2,
    pos = 2,
    labels = paste0("Upper CUSUM: ", round(ucl, 3)), 
    col = "#e06666",
    font = 2
  )
  text(
    length(dt),
    lcl * 1.2,
    pos = 2,
    labels = paste0("Lower CUSUM: ", round(lcl, 3)), 
    col = "#e06666",
    font = 2
  )
}

# generate random data for scenario 7.1 (under control)
data.7.1 <- rnorm(100, 1, .05)
cc.CUSUM(data.7.1, .05)

# generate random data for scenario 7.2 (out of control)
data.7.2 <- c(rnorm(50, 1, .05), rnorm(50, 1.05, .05))
cc.CUSUM(data.7.2, .05)
