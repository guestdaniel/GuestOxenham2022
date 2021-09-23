# Libraries
library(ggplot2)
library(dplyr)
library(lme4)
library(merTools)
library(effects)
library(phia)
options(contrasts=c("contr.sum","contr.poly"))

# Load data
load('data/exp1.RData')
data = data_exp1

# Model
model = lmer(threshold ~ F0*masker*experiment + (F0*masker|subj), data=data)

# Preprocess data for visualization
temp = data
temp_ind = temp %>%
		# Data preprocessing
		group_by(F0, masker, subj, experiment) %>%
		summarize(threshold_low = 10^((mean(threshold)-1.00*sd(threshold)/sqrt(n()))/10),
				  threshold_high = 10^((mean(threshold)+1.00*sd(threshold)/sqrt(n()))/10),
				  threshold_mean = 10^(mean(threshold)/10))
temp = temp %>%
		# Data preprocessing
		group_by(F0, masker, subj, experiment) %>%
		summarize(threshold = mean(threshold, na.rm=TRUE)) %>%
		group_by(F0, masker, experiment) %>%
		summarize(threshold_low = 10^((mean(threshold)-1.00*sd(threshold)/sqrt(n()))/10),
				  threshold_high = 10^((mean(threshold)+1.00*sd(threshold)/sqrt(n()))/10),
				  threshold_mean = 10^(mean(threshold)/10))

# Combine preprocessed data with model predictions
temp$source = "Data"
model_predictions = as.data.frame(allEffects(model))[[1]]
model_predictions$threshold_mean = 10^(model_predictions$fit/10)
model_predictions$source = "Model"
temp = rbind(temp, model_predictions)

# Plot
font_scale=8
size_point = 3
size_line = 2
size_point_sub = 1
size_line_sub = 0.5
size_error = 1
width_error = 0.06

# Define breaks for plotting logarithmic axes
breaks_ST = 1/(1/2^(-7:2 + 0.46))
n_breaks = 10
breaks = c(seq(2, 10, 10/n_breaks) %o% 10^(-3:3))
labels = as.character(breaks)
labels[!(log10(breaks)%%1==0)] = ''
ticksizes = rep(.25, length(breaks))
ticksizes[log10(breaks)%%1==0] = 1

# Aesthetics
temp %>%
	ggplot(aes(x=masker, y=threshold_mean, color=F0, group=F0, shape=source)) +
	# Geoms
	geom_hline(yintercept=c(100*2^(breaks_ST/12) - 100), linetype="dashed", color="gray") +
	geom_point(size=size_point) +
	geom_errorbar(size=size_error, width=width_error, aes(ymin=threshold_low, ymax=threshold_high)) +
	# Modify axes
	scale_y_log10(breaks=breaks, labels=labels, limits=c(1.5/10, 30)) +
	# Set theme and guide settings
	theme_bw() +
	guides(shape=guide_legend(
					keywidth=0.1,
					keyheight=0.4,
					default.unit="inch")) +
	# Configure panel grids, ticks, and font sizes
	theme(axis.text.y=element_text(size=1*font_scale),           # axis tick label font size
			  axis.text.x=element_text(size=1*font_scale),
			  axis.title.y=element_text(size=1.2*font_scale),    # axis label font size
			  axis.title.x=element_text(size=1.2*font_scale),
			  legend.text=element_text(size=1*font_scale),       # legend text font size
			  legend.title=element_text(size=1.2*font_scale),    # legend title font size
			  strip.text.x=element_text(size=1*font_scale),      # facet label font size
			  plot.title=element_text(size=1.5*font_scale),      # figure title font size
			  panel.grid.major=element_blank(),
			  panel.grid.minor = element_blank(),
			  axis.ticks.x=element_line(size=ticksizes)) +
	# Configure scales
	scale_color_manual(values=c("#7fc97f", "#a28ac1"), guide=guide_legend(title="F0 (Hz)")) +
	scale_shape_manual(values=c(16, 24)) +
	# Titles and labels
	xlab("Masker Type") +
	ylab("F0DL (%)") +
	# Facets
	facet_grid(. ~ experiment)
