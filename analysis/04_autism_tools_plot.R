# 2018-05-04
# Compare-surf-tools: recreate plot of CT ~ Dx for each toolkit (ANTS, FS5.1, FS5.3)
# see ~/abide/compare-surf-tools/results/AutismTools.png
# The orig plot is from Erin Dickie
# ? Was -log10(p) intended to be plotted, instead of adjust the scale to log10?
# ? NOTE that the result is not signed (the direction of the coefficient estimate is not assigned to the p-values).
# --------------

# --------------
# Load required software
library(ggplot2)
# --------------

# --------------
# set working directory and load the data
# setwd("/Volumes/hodges/abide/compare-surf-tools/")
setwd("~/abide/compare-surf-tools/")
load("analysis/abide_ct.RData")
# --------------

# --------------
# Add level names to the DX_GROUP factor
abide$DX_GROUP <- factor(abide$DX_GROUP, levels = c(1, 2), labels = c("ASD", "TD"))
# NOTE: these labels are consistent with the file 'dx_effect_testing.Rmd" on the github repo
# --------------

# --------------
# run the models
# thickness ~ DX_GROUP ... for each level of "method" for each PU in "label_abbrev"
# extract the p-value from the DX_GROUP coefficient.

tools <- unique(abide$method)
pus <- unique(abide$label_abbrev)

# start a results matrix: tool, pu, p-value
reslist <- list()

# run the models
for(j in seq_along(tools)){
	for(i in seq_along(pus)){
		
		dat <- subset(abide, method == tools[j] & label_abbrev == pus[i], select = c("SubjID", "method", "label_abbrev", "thickness", "DX_GROUP"))
		lm_i <- lm(thickness ~ DX_GROUP, data = dat)
		reslist[[j * length(pus) + i - length(pus)]] <- cbind(tools[j], pus[i], summary(lm_i)$coef[2, 1], summary(lm_i)$coef[2, 4])
	}
}
# --------------


# --------------
# Convert the results to a data.frame
row.names(reslist) <- NULL
results_df <- data.frame(do.call(rbind, reslist), stringsAsFactors = FALSE)
names(results_df) <- c("tool", "pus", "coef", "pval")
results_df$coef <- as.numeric(results_df$coef)
results_df$pval <- as.numeric(results_df$pval)
results_df$pvallog <- -1 * log10(results_df$pval)
# --------------


# --------------
# plot the results
# For DK: just ANTS and FS53
gg <- ggplot(results_df[results_df$tool %in% c("ANTS", "FS53"),], aes(x = pus, y = pvallog, color = tool)) + 
geom_point() + 
scale_color_brewer(palette = "Set1") + 
scale_y_reverse() + 
scale_y_log10() + 
coord_flip() + 
labs(title = "Effect of Diagnosis, no covariates", x = NULL, y = "Inverse Log p-value (uncorrected)") + 
facet_grid(. ~ tool) +
guides(color = FALSE)

gg
ggsave(gg, filename = "results/AutismTools_redux_redux.png", height = 8, width = 7)
# --------------

# --------------
# Variation of the plot:
# plot inverse log p-value, rather than log scale. Add reference line at -log10(0.05) = 1.3
gg <- ggplot(results_df[results_df$tool %in% c("ANTS", "FS53"),], aes(x = pus, y = pvallog, color = tool)) + 
  geom_point() + 
  scale_color_brewer(palette = "Set1") + 
# scale_y_reverse() + 
# scale_y_log10() + 
  coord_flip() + 
  labs(title = "Effect of Diagnosis, no covariates", x = NULL, y = "Inverse Log p-value (uncorrected)", subtitle = "Model: Volume ~ Dx ... run on data from each tool in [ANTS, FS5.1, FS5.3]") + 
  facet_grid(. ~ tool) +
  guides(color = FALSE) +
# add line at p = 0.05 ... and adjust the number of tick marks and axis labels.
  geom_hline(yintercept = -1*log10(0.05), lty = 3) 
# scale_y_log10(breaks = pretty(results_df$pval[results_df$tool %in% c("ANTS", "FS53")], n = 15))

ggsave(gg, filename = "results/AutismTools_redux_pline3.png", height = 8, width = 7)
# --------------