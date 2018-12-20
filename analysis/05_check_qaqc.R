# 2018-12-18
# For the ABIDE 'tools matter' analysis of ANTS and FS, see whether the comparison of the volumes from each tool are dependent in some way on the QA/QC data.
# Possible proxies for the poor corr among methods besides QA/QC:
# DX_GROUP
# site (same as SITE_ID?)
# AGE_AT_SCAN (should it be mean centered?)

# QA/QC measures:
# anat_cnr 
# anat_efc
# anat_fber
# anat_fwhm
# anat_qi1
# anat_snr
# NOTE: of these, anat_cnr and anat_snr are priority for now.

# qc_rater_1 | qc_anat_rater_2 | qc_anat_rater_3
# ----------------

# ----------------
# Load the required packages

# ----------------

# ----------------
# set the working dir to the comp-surf-tools directory
setwd("~/abide/compare-surf-tools/.")
# ----------------

# ----------------
# Load the data
load("analysis/abide_ct.RData")
# 'abide
# ----------------

# ----------------
# Look at the correlations among the qc measures
# Define the plot functions
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "orange", ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = "pairwise.complete.obs")
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    # text(0.5, 0.5, txt, cex = cex.cor * r)
    text(0.5, 0.5, txt, cex = 2)
}
panel.smooth.ab <- function(x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) {
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), col = col.smooth, ...)
        abline(a = 0, b = 1, lty = 3, col = "gray50")
}
# ----------------

# ----------------
# Get a collection of qc data from one roi and method ... and a second dataset to check the first (should be the same data)
abide1 <- subset(abide, method %in% "ANTS" & label_abbrev %in% "L.CN")
abide2 <- subset(abide, method %in% "FS53" & label_abbrev %in% "R.SM")

dim(abide1)
# [1] 1101  114
dim(abide2)
# [1] 976 114
# ... the datasets are of differing length
# This is because not all of the data were run through FS:
head(table(abide$label_abbrev, abide$method))
            
# #              ANTS FS51 FS53
  # L.cauantCG 1101 1112  976
  # L.caumidF  1101 1112  976
  # L.CN       1101 1112  976
  # L.ENT      1101 1112  976
  # L.FUS      1101 1112  976
  # L.infP     1101 1112  976

# ----------------

# ----------------
# qc meas corrs
qa <- c('anat_cnr', 'anat_efc', 'anat_fber', 'anat_fwhm', 'anat_qi1', 'anat_snr')
pairs(
  abide1[, qa],
  gap = 0,
  lower.panel = panel.smooth.ab,
  upper.panel = panel.cor,
  diag.panel = panel.hist,
  pch = 21,
  bg = as.factor(abide1$site),
  main = "Title"
)
dev.new()
pairs(
  abide2[, qa],
  gap = 0,
  lower.panel = panel.smooth.ab,
  upper.panel = panel.cor,
  diag.panel = panel.hist,
  pch = 21,
  bg = as.factor(abide1$site),
  main = "Title"
)
# hmm ... histograms show measures that are skewed ... and the datasets are not the same (though the summary values are not too different)

summary(abide1$anat_cnr)
summary(abide1$anat_snr)
summary(abide1$anat_fwhm)
summary(abide1$anat_qi1)

dev.new()
hist(abide1$anat_cnr)
hist(abide1$anat_snr)
hist(abide1$anat_fwhm)
hist(abide1$anat_qi1)

# these work
hist(sqrt(abide1$anat_cnr))
hist(sqrt(abide1$anat_snr)) # <- but not for this one
hist(sqrt(abide1$anat_fwhm))
hist(sqrt(abide1$anat_qi1))

# This is okay:
hist(abide1$anat_snr[abide1$anat_snr < 30])
hist(abide1$anat_cnr[abide1$anat_cnr < 20])

# how many data points do these exclude?
sum(!is.na(abide1$anat_snr)) - sum(abide1$anat_snr < 30, na.rm = TRUE)
# [1] 124
sum(!is.na(abide1$anat_snr)) - sum(abide1$anat_cnr < 20, na.rm = TRUE)
# [1] 70
# SO ... not too many of the n = 1098 available (and non-missing) for ANTS
# ----------------

# ----------------
# Get a list of the regions of interst
roi <- sort(unique(abide$label_abbrev))

# Get a subset of the ANTS and FS53 data, including only some of the columns of interest:
voi <- c("SubjID", "thickness",  "method", "label_abbrev", "site", "anat_cnr", "anat_snr")
abides <- subset(abide, method %in% c("ANTS", "FS53"), select = voi)

# reshape into wide format by method
# NOT WORKING AS EXPECTED ... MOVE ALONG
# abidesw <- reshape(abides, direction = "wide", idvar = c("SubjID", "label_abbrev", "site", "anat_cnr", "anat_snr"), timevar = "method")
# ----------------

# ----------------
# create loop of plots and correlation results:
# For each roi:
# ANTS ~ FS53, abs(resid) ~ anat_cnr, abs(resid ~ anat_snr)
# for pairwise complete cases and where anat_snr < 30 and abide1$anat_cnr < 20
# text output: roi, n, coef, p.value, cnr coef, p.value, snr coef, p.value
# graphic output: y ~ x, abs(resid) ~ cnr, abs(resid ~ snr) ... with fit lines, pch color by site?

# a long example
roi1 <- subset(abides, label_abbrev %in% roi[1])
roi1w <- reshape(roi1, direction = "wide", idvar = c("SubjID", "label_abbrev", "site", "anat_cnr", "anat_snr"), timevar = "method")
# pare down to complete cases AND roi1w$anat_cnr < 20 & roi1w$anat_snr < 30
cc <- complete.cases(roi1w[, c("thickness.ANTS", "thickness.FS53")])
roi1ws <- subset(roi1w, complete.cases(roi1w[, c("thickness.ANTS", "thickness.FS53")]) & roi1w$anat_cnr < 20 & roi1w$anat_snr < 30)
lm1 <- lm(thickness.ANTS ~ thickness.FS53, data = roi1ws)
rlmcnr <- lm(abs(resid(lm1)) ~ roi1ws$anat_cnr)
rlmsnr <- lm(abs(resid(lm1)) ~ roi1ws$anat_snr)

mylim <- range(c(roi1ws$thickness.ANTS, roi1ws$thickness.FS5))

pdf(paste0("results/qa/", roi[1], ".pdf"), height = 3.75, width = 10)
# dev.new(height = 3.75, width = 10)
par(mfrow = c(1, 3))
plot(
	thickness.ANTS ~ thickness.FS53,
	data = roi1ws,
	pch = 21,
	col = "gray50",
	bg = as.factor(roi1ws$site),
	xlim = mylim,
	ylim = mylim,
	main = roi[1]
)
abline(lm1, col = "red", lwd = 2)
abline(a = 0, b = 1, lty = 3)

plot(
	abs(resid(lm1)) ~ roi1ws$anat_cnr,
	pch = 21,
	col = "gray50",
	bg = as.factor(roi1ws$site),
	xlab = "CNR",
	ylab = "abs(Residuals)",
	main = roi[1]
)
abline(rlmcnr, col = "red", lwd = 2)

plot(
	abs(resid(lm1)) ~ roi1ws$anat_snr,
	pch = 21,
	col = "gray50",
	bg = as.factor(roi1ws$site),
	xlab = "SNR",
	ylab = "abs(Residuals)",
	main = roi[1]
)
abline(rlmsnr, col = "red", lwd = 2)
dev.off()

# ----------------
# For PDF printing:
# pdf(paste0("results/qa/", roi[i], ".pdf"), height = 7, width = 7)
# [...]
# dev.off()


# Create a results table for later:
# roi, lm1.coef, lm1p, cnr.coef, cnr.p, snr.coef, snr.p
out_list <- matrix(NA, nrow = length(roi), ncol = 7)
out_list[1, ] <- c(
	roi[1], 
	summary(lm1)$coef[2, 1],
	summary(lm1)$coef[2, 4],
	summary(rlmcnr)$coef[2, 1],
	summary(rlmcnr)$coef[2, 4],
	summary(rlmsnr)$coef[2, 1],
	summary(rlmsnr)$coef[2, 4]
)
# save it to "results/qa/."
# ----------------


# ----------------
# Now do it in a loop
# Create a results table for later:
# roi, lm1.coef, lm1p, cnr.coef, cnr.p, snr.coef, snr.p
out_list <- matrix(NA, nrow = length(roi), ncol = 7)

for(i in seq_along(roi)){
	cat(roi[i], "\n")

	# extract and prep the data
	df <- subset(abides, label_abbrev %in% roi[i])
	dfw <- reshape(df, direction = "wide", idvar = c("SubjID", "label_abbrev", "site", "anat_cnr", "anat_snr"), timevar = "method")
	# pare down to complete cases AND roi1w$anat_cnr < 20 & roi1w$anat_snr < 30
	dfws <- subset(dfw, complete.cases(dfw[, c("thickness.ANTS", "thickness.FS53")]) & dfw$anat_cnr < 20 & dfw$anat_snr < 30)

	# run the models
	lm1 <- lm(thickness.ANTS ~ thickness.FS53, data = dfws)
	rlmcnr <- lm(abs(resid(lm1)) ~ dfws$anat_cnr)
	rlmsnr <- lm(abs(resid(lm1)) ~ dfws$anat_snr)

	# write the model results to the output file
	out_list[i, ] <- c(
		roi[i], 
		summary(lm1)$r.square,
		summary(lm1)$coef[2, 4],
		summary(rlmcnr)$r.square,
		summary(rlmcnr)$coef[2, 4],
		summary(rlmsnr)$r.square,
		summary(rlmsnr)$coef[2, 4]
	)

	# # plot the results
	mylim <- range(c(dfws$thickness.ANTS, dfws$thickness.FS5))

	pdf(paste0("results/qa/", roi[i], ".pdf"), height = 3.75, width = 10)
	# dev.new(height = 3.75, width = 10)
	par(mfrow = c(1, 3))
	plot(
		thickness.ANTS ~ thickness.FS53,
		data = dfws,
		pch = 21,
		col = "gray50",
		bg = as.factor(dfws$site),
		xlim = mylim,
		ylim = mylim,
		main = roi[i]
	)
	abline(lm1, col = "red", lwd = 2)
	abline(a = 0, b = 1, lty = 3)

	plot(
		abs(resid(lm1)) ~ dfws$anat_cnr,
		pch = 21,
		col = "gray50",
		bg = as.factor(dfws$site),
		xlab = "CNR",
		ylab = "abs(Residuals)",
		main = roi[i]
	)
	abline(rlmcnr, col = "red", lwd = 2)

	plot(
		abs(resid(lm1)) ~ dfws$anat_snr,
		pch = 21,
		col = "gray50",
		bg = as.factor(dfws$site),
		xlab = "SNR",
		ylab = "abs(Residuals)",
		main = roi[i]
	)
	abline(rlmsnr, col = "red", lwd = 2)
	dev.off()
}

# adjust the output results
out_df <- data.frame(out_list, stringsAsFactors = FALSE)
out_df[,-1] <- sapply(out_df[, -1], as.numeric)
# save it to "results/qa/."
names(out_df) <- c("roi", "method.rsqr", "method.p", "cnr.rsqr", "cnr.p", "snr.rsqr", "snr.p")
write.csv(out_df, file = "results/qa/results2.csv", row.names = FALSE)

# ----------------


# ----------------
# Try using non-abs(residual)
# Get a list of the regions of interst
roi <- sort(unique(abide$label_abbrev))

# Get a subset of the ANTS and FS53 data, including only some of the columns of interest:
voi <- c("SubjID", "thickness",  "method", "label_abbrev", "site", "anat_cnr", "anat_snr")
abides <- subset(abide, method %in% c("ANTS", "FS53"), select = voi)

# Create a results table for later:
# roi, lm1.coef, lm1p, cnr.coef, cnr.p, snr.coef, snr.p
out_list <- matrix(NA, nrow = length(roi), ncol = 7)

for(i in seq_along(roi)){
	cat(roi[i], "\n")

	# extract and prep the data
	df <- subset(abides, label_abbrev %in% roi[i])
	dfw <- reshape(df, direction = "wide", idvar = c("SubjID", "label_abbrev", "site", "anat_cnr", "anat_snr"), timevar = "method")
	# pare down to complete cases AND roi1w$anat_cnr < 20 & roi1w$anat_snr < 30
	dfws <- subset(dfw, complete.cases(dfw[, c("thickness.ANTS", "thickness.FS53")]) & dfw$anat_cnr < 20 & dfw$anat_snr < 30)

	# run the models
	lm1 <- lm(thickness.ANTS ~ thickness.FS53, data = dfws)
	rlmcnr <- lm(resid(lm1) ~ dfws$anat_cnr)
	rlmsnr <- lm(resid(lm1) ~ dfws$anat_snr)

	# write the model results to the output file
	out_list[i, ] <- c(
		roi[i], 
		summary(lm1)$r.square,
		summary(lm1)$coef[2, 4],
		summary(rlmcnr)$r.square,
		summary(rlmcnr)$coef[2, 4],
		summary(rlmsnr)$r.square,
		summary(rlmsnr)$coef[2, 4]
	)

	# # plot the results
	mylim <- range(c(dfws$thickness.ANTS, dfws$thickness.FS5))

	pdf(paste0("results/qa/", roi[i], "_resid.pdf"), height = 3.75, width = 10)
	# dev.new(height = 3.75, width = 10)
	par(mfrow = c(1, 3))
	plot(
		thickness.ANTS ~ thickness.FS53,
		data = dfws,
		pch = 21,
		col = "gray50",
		bg = as.factor(dfws$site),
		xlim = mylim,
		ylim = mylim,
		main = roi[i]
	)
	abline(lm1, col = "red", lwd = 2)
	abline(a = 0, b = 1, lty = 3)

	plot(
		resid(lm1) ~ dfws$anat_cnr,
		pch = 21,
		col = "gray50",
		bg = as.factor(dfws$site),
		xlab = "CNR",
		ylab = "Residuals",
		main = roi[i]
	)
	abline(rlmcnr, col = "red", lwd = 2)

	plot(
		resid(lm1) ~ dfws$anat_snr,
		pch = 21,
		col = "gray50",
		bg = as.factor(dfws$site),
		xlab = "SNR",
		ylab = "Residuals",
		main = roi[i]
	)
	abline(rlmsnr, col = "red", lwd = 2)
	dev.off()
}

# adjust the output results
out_df <- data.frame(out_list, stringsAsFactors = FALSE)
out_df[,-1] <- sapply(out_df[, -1], as.numeric)
# save it to "results/qa/."
names(out_df) <- c("roi", "method.rsqr", "method.p", "cnr.rsqr", "cnr.p", "snr.rsqr", "snr.p")
write.csv(out_df, file = "results/qa/results3.csv", row.names = FALSE)

# ----------------