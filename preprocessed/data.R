library(RSQLite)

## Fetch aggregate tables from the SQLite database
con <- dbConnect(SQLite(), "data.db")

rs <- dbSendQuery(con, "SELECT COUNT(*),age,year FROM data GROUP BY age,year ORDER BY year")
total_births <- dbFetch(rs, -1)

rs2 <- dbSendQuery(con, "SELECT COUNT(*),age,year FROM data WHERE children > 1 GROUP BY age,year ORDER BY year")
twin_births <- dbFetch(rs2, -1)

dbDisconnect(con)

## Data treatment
library(plyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(tikzDevice)
path <- ""

names(total_births) <- names(twin_births) <- c("count", "age", "year")

## Some extreme values for mother's age are only present in certain years. This
## code snipped harmonizes age to be between 12 and 49 while discarding the rest
## and filling up zeroes if necessary. It's much faster to do it here rather
## than directly via SQL.

total_births <- subset(total_births, age >= 12 & age <= 49)
twin_births <- subset(twin_births, age>=12 & age <= 49)

unifData <- function(x) {
  year_span <- 1971:2015
  for (i in 12:49) {
    ind <- which(!(year_span %in% subset(x, age == i)$year))
    if (length(ind) > 0) {
      x <- rbind(x, cbind(year=year_span[ind], age=i, count=0))
    }
  }
  return(arrange(x, year, age))
}

total_births <- unifData(total_births)
twin_births <- unifData(twin_births)

## Create the density curves per year
total_births <- ddply(total_births, .(year), mutate, density=count/sum(count))
twin_births <- ddply(twin_births, .(year), mutate, density=count/sum(count))

## Plot the fertility density curves
plotf <- function(x, title, legend, ylab) {
  ifelse(legend, position <- "right", position <- "none")
  p <- ggplot(x, aes(x = age, y = density, group = year, col = year)) +
    geom_line() +
    xlab("Mother's age") + ylab(ylab) +
    theme_bw() +
    theme(text = element_text(size = 11),
          plot.title = element_text(size=10),
          legend.position = position,
          legend.title = element_text(size=10)) +
    scale_colour_continuous(name = "Year", breaks=c(1971, 1980, 1990, 2000, 2010),
                            low = "green", high = "blue")
  ggtitle(title)
  return(p)
}


tikz(file = paste0(path, "fertility_plots.tex"), width = 7, height = 4)
plot1 <- plotf(total_births, "Total fertility density", FALSE, "Density")
plot2 <- plotf(twin_births, "Fertility density for multiple births", TRUE, "")
grid.arrange(plot1, plot2, ncol = 2)
dev.off()
## There's an issue with raster legend for this plot in TeX compilation.
## For a workaround, generate .tex and .png files from above, then replace
## \pgfimage with \includegraphics and add code/ to the path in the generated
## TeX file (~ line 3891).


## Get GDP growth data
library(quantmod)
getSymbols("A939RC0A052NBEA", src="FRED")
gdp <- A939RC0A052NBEA
gdpg <- (gdp/lag(gdp,1) - 1)*100
y <- as.numeric(gdpg[time(gdpg) > "1970-01-01" & time(gdpg) < "2016-01-01"])
## Apply functional instrumental regression estimators
library(functionalIV)
years <- 1971:2015
age <- 12:49
x <- matrix(total_births$density, nrow = length(years),
            ncol = length(age), byrow = TRUE)
w <- matrix(twin_births$density, nrow = length(years),
            ncol = length(age), byrow = TRUE)
## Wavelets need domain length to be a power of two
age_wav <- age[4:35]
wav_x <- x[, 4:35]
wav_w <- w[, 4:35]

fit <- endogenousEstimation(x = t(x)*1000, w = t(w)*1000, y = y,
                            basis = "fourier", m = 7, domain = c(0,1),
                            PCA = TRUE, pca.dim = 6,
                            tikhonov = TRUE, tikhonov.alpha = 0.05)

fitw <- waveletEndogenous(x = t(wav_x)*10000, w = t(wav_w)*10000, y = y,
                          minScale = 3, level = 3, family = "DaubLeAsymm",
                          filter.number = 4)

## Plot Fourier, Tikhonov and wavelet estimators
age_bounds <- c(15, 35)
age_span <- which(age >= age_bounds[1] & age <= age_bounds[2])
age_span_wav <- which(age_wav >= age_bounds[1] & age_wav <= age_bounds[2])
allEstimators <- as.data.frame(cbind(age=age[age_span], galerkin.fourier=fit$galest_endoSA[age_span],
                                     tikhonov=fit$tikhonov_est[age_span], # pca=fit$pca_est[age_span],
                                     wavelet.thresh=fitw$wav_thresh_est[age_span_wav]))
allEstimators <- melt(allEstimators, id.var="age")
allEstimators$variable <- factor(allEstimators$variable, levels=c("galerkin.fourier", "tikhonov", "wavelet.thresh"),
                                 labels=c("Fourier basis estimator",
                                          "Tikhonov estimator",
                                          "Wavelet thresholded estimator"))

tikz(file = paste0(path, "estimator_plots.tex"), width = 6, height = 5)
p <- ggplot(allEstimators, aes(x = age, y = value)) +
  geom_line(aes(color = variable, linetype = variable)) +
  geom_point(aes(x = age, y = value, color = variable, shape = variable), size=3) +
  xlab("Mother's age") + ylab("") +
  theme_bw() +
  theme(text = element_text(size = 11),
        plot.title = element_text(size = 10),
        legend.justification = c(1.05, 0),
        legend.position = c(1, 0),
        legend.title = element_blank(),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.key.width=unit(2.5, "cm")) +
  guides(col = guide_legend(ncol = 1))+
  scale_color_manual(values = c("blue", "darkgreen", "darkred", "black")) +
  ggtitle("Estimates of regressand function on growth of GDP/capita")
print(p)
dev.off()

## Fourier confidence interval
getPrediction <- function(fit, xx, ww) {

  step <- 1/(nrow(xx)-1)
  n <- ncol(xx)

  sigma.f <- as.numeric(var(y - (1000*x %*% fit$galest_endoSA)*step))
  That <- step*(ww %*% t(xx))/n; Sigma <- step*(ww %*% t(ww))/n; TST <- That %*% Sigma %*% t(That)
  project <- t(fit$basis_values) %*% TST %*% fit$basis_values * step^3
  gm <- diag(1, 38) %*% (fit$basis_values)*step
  aa <- sigma.f * solve(fit$gal_matrix) %*% project %*% solve(fit$gal_matrix) / n
  bound <- 1.96*sqrt(diag(gm %*% aa %*% t(gm)))

  pred <- as.data.frame(cbind("age" = age[age_span],
                              "galerkin.fourier" = fit$galest_endoSA[age_span]))
  pred <- melt(pred, id.var = "age")
  pred$variable <- factor(pred$variable,
                          levels=c("galerkin.fourier"),
                          labels=c("Fourier basis estimator"))
  pred$upper <- (fit$galest_endoSA + bound)[age_span]
  pred$lower <- (fit$galest_endoSA - bound)[age_span]
  pred
}

pred <- getPrediction(fit, t(x) * 1000, t(w) * 1000)

tikz(file = paste0(path, "confidence_bands.tex"), width = 5, height = 3)
p <- ggplot(pred, aes(x = age, y = value)) +
  geom_line(data = pred) +
  geom_smooth(aes(ymin = lower, ymax = upper), stat = "identity") +
  geom_point(aes(x = age, y = value, color = variable, shape = variable), size = 3) +
  xlab("Mothers' age") + ylab("") +
  theme_bw() +
  theme(text = element_text(size = 11),
        plot.title = element_text(size = 10),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key.width=unit(2.5, "cm")) +
  scale_color_manual(values=c("blue", "darkgreen", "darkred", "black")) +
  scale_fill_manual("", values = "grey75") +
  scale_y_continuous(limits=c(-4.8, 3.5)) +
  ggtitle("Pointwise confidence interval for Fourier basis estimate")
print(p)
dev.off()
