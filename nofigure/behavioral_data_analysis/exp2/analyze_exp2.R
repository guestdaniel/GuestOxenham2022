# Libraries
library(ggplot2)
library(dplyr)
library(lme4)
library(merTools)
library(effects)
library(phia)
options(contrasts=c("contr.sum","contr.poly"))

# Define location to store results
model_dir = 'nofigure/behavioral_data_analysis/exp2/'

# Define function to jointly correct p-values in a list of model/contrast objects
correct_p <- function(objs) {
	# Identify number of rows and the column containing p-values in each object
	n_rows = sapply(objs, nrow)
	col_names = sapply(objs, colnames)
	col_arrays = sapply(sapply(col_names, FUN=stringr::str_detect, pattern="Pr\\(>"), as.numeric)
	col_idxs = sapply(col_arrays, which.max)
	# Extract p values
	p_values = numeric()
	for (ii in 1:length(objs)) {
			p_values = c(p_values, objs[[ii]][, col_idxs[ii]])
	}
	# Correct p values
	p_values = p.adjust(p_values, method="holm")
	# Replace p values
	idx = 1
	for (ii in 1:length(objs)) {
			objs[[ii]][, col_idxs[ii]] = p_values[idx:(idx+n_rows[ii]-1)]
			idx = idx + n_rows[[ii]]
	}
	return(objs)
}

# Load data
load('data/exp2.RData')
data = data_exp2

# Fit mixed effects model to data
model = lmer(threshold ~ F0*interval + (F0*interval|subj), data=data,
		control=lmerControl(optimizer="optimx",
			optCtrl=list(method="nlminb",
			starttests=FALSE,
			kkt=FALSE)))

# Plot various diagnostic plots
# Residual plot
png(paste0(model_dir, "residuals_vs_fitted_values.png"))
	plot(model) # Residuals
dev.off()
# Scale location
png(paste0(model_dir, "scale_location.png"))
plot(model, sqrt(abs(resid(.)))~fitted(.),
	 type=c("p", "smooth"),
	 ylab=expression(sqrt(abs(resid))),
	 main="Scale Location Plot (sqrt(resid) vs fitted)")
dev.off()
# Normalized QQ plot
png(paste0(model_dir, "qq.png"))
qqnorm(scale(resid(model)),
	 main="QQ Plot of Normalized Resiudals")
abline(c(0,1))
dev.off()
# Residual vs condition
png(paste0(model_dir, "residuals_vs_fitted_values_by_condition.png"))
	boxplot(residuals(model) ~ interaction(data$F0, data$interval)) # Residuals
dev.off()

# Analyze model with ANOVA
model_anova = Anova(model, type=3, test="F")
model_anova[, "Pr(>F)"] = p.adjust(model_anova[, "Pr(>F)"], method="holm")

# Calculate various contrast tests
contrast_F0 = testInteractions(model, pairwise="F0", adjustment="none", test="F")
contrast_interval = testInteractions(model, pairwise="interval", adjustment="none", test="F")
contrast_interval_by_F0 = testInteractions(model, pairwise="interval", fixed="F0", adjustment="none", test="F")
contrast_interval_and_F0 = testInteractions(model, pairwise=c("interval", "F0"), adjustment="none", test="F")

# Jointly correct contrast test p-values
contrasts = list(contrast_F0, contrast_interval, contrast_interval_by_F0, contrast_interval_and_F0)
unadj_values = numeric()
for (ii in 1:4) {
	unadj_values = c(unadj_values, contrasts[[ii]][, "Pr(>F)"])
}
adj_values = p.adjust(unadj_values, method="holm")
contrast_F0[, "Pr(>F)"] = adj_values[1:2]
contrast_interval[, "Pr(>F)"] = adj_values[3:4]
contrast_interval_by_F0[, "Pr(>F)"] = adj_values[5:7]
contrast_interval_and_F0[, "Pr(>F)"] = adj_values[8:9]

# Sink model results to disk
sink(paste0(model_dir, "model_output.txt"))
"ANOVA // Satterthwaite's Estimation of denominator DF // Joint correction by Holm-Bonferroni"
model_anova
"Contrasts to test F0"
contrast_F0
"Contrasts to test interval"
contrast_interval
"Contrasts to test interval by F0"
contrast_interval_by_F0
"Contrasts to test interval F0"
contrast_interval_and_F0
sink()