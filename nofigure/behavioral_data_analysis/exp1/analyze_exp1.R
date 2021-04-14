# Libraries
library(ggplot2)
library(dplyr)
library(lme4)
library(merTools)
library(effects)
library(phia)
options(contrasts=c("contr.sum","contr.poly"))

# Define location to store results
model_dir = 'nofigure/behavioral_data_analysis/exp1/'

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
load('data/exp1.RData')
data = data_exp1

# Model
model = lmer(threshold ~ F0*masker*experiment + (F0*masker|subj), data=data)

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
	boxplot(residuals(model) ~ interaction(data$F0, data$masker)) # Residuals
dev.off()
# Interaction means
png(paste0(model_dir, "interaction_means.png"))
	plot(interactionMeans(model))
dev.off()

# Analyze model with ANOVA
model_anova = Anova(model, type=3, test="F")
model_anova[, "Pr(>F)"] = p.adjust(model_anova[, "Pr(>F)"], method="holm")

# Compute contrast tests
test_F0 = testInteractions(model, pairwise="F0", fixed="masker", adjustment="none", test="F")
test_masker = testInteractions(model, pairwise="masker", fixed="F0", adjustment="none", test="F")
test_F0_and_masker = testInteractions(model, pairwise=c("F0", "masker"), adjustment="none", test="F")
#test_experiment_pairwise = testInteractions(model, pairwise="experiment", fixed=c("masker"), adjustment="none",
#											test="F")

# Combine contrast tests and correct their p-values
contrasts = list(test_F0, test_masker, test_F0_and_masker, test_experiment_pairwise)
unadj_values = numeric()
for (ii in 1:3) {
	unadj_values = c(unadj_values, contrasts[[ii]][, "Pr(>F)"])
}
adj_values = p.adjust(unadj_values, method="holm")
test_F0[, "Pr(>F)"] = adj_values[1:3]
test_masker[, "Pr(>F)"] = adj_values[4:6]
test_F0_and_masker[, "Pr(>F)"] = adj_values[7:8]
#test_experiment_pairwise[, "Pr(>F)"] = adj_values[9:13]

# Change all the difference measures to ratio measures in the contrast tests
test_F0[, "Value"] = 10^(abs(test_F0[, "Value"])/10)
test_masker[, "Value"] = 10^(abs(test_masker[, "Value"])/10)
test_F0_and_masker[, "Value"] = 10^(abs(test_F0_and_masker[, "Value"])/10)
#test_experiment_pairwise[, "Value"] = 10^(abs(test_experiment_pairwise[, "Value"])/10)

# Sink the results to file
sink(paste0(model_dir, "model_output.txt"))
"ANOVA // KR approximation of denominator DF // Joint correction by Holm-Bonferroni"
model_anova
"Contrasts to test F0 by masker // KR approximation of denominator DF // Joint correction by Holm-Bonferroni"
test_F0
"Contrasts to test masker by F0 // KR approximation of denominator DF // Joint correction by Holm-Bonferroni"
test_masker
"Contrasts to test F0-masker interaction // KR approximation of denominator DF // Joint correction by Holm-Bonferroni"
test_F0_and_masker
#"Contrasts to test experiment pairwise // KR approximation of denominator DF // Joint correction by Holm-Bonferroni"
#test_experiment_pairwise
sink()