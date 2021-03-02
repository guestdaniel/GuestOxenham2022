# Libraries
library(ggplot2)
library(dplyr)
library(lme4)
library(car)
library(phia)
library(effects)
library(lmerTest)

# Global parameters
font_scale=8
size_point = 3
size_line = 2
size_point_sub = 1
size_line_sub = 0.5
size_error = 1
width_error = 0.06

# Load data
load('data/exp2.RData')
data = data_exp2

# Calculate means and standard errors
data$interval = factor(data$interval)
temp_ind = data %>%
	group_by(subj, F0, interval) %>%
	summarize(threshold_sd = sd(threshold)/sqrt(n()), threshold_mean = mean(threshold))
temp_avg = temp_ind %>%
	group_by(F0, interval) %>%
	summarize(threshold_sd = sd(threshold_mean)/sqrt(n()), threshold_mean = mean(threshold_mean))

# Plot
temp_avg %>%
    # Aesthetics
	ggplot(aes(x=interval, y=threshold_mean, ymin=threshold_mean-threshold_sd, ymax=threshold_mean+threshold_sd, color=F0, group=F0)) +
	# Geoms
	geom_hline(yintercept=seq(0, 16, by=4), linetype="dashed", color="gray") +
	geom_point(size=size_point) +
	geom_line(size=size_line) +
	geom_errorbar(size=size_error, width=width_error) +
	geom_pointrange(size=size_line_sub, fatten=size_point_sub, position=position_dodge(width=0.10), data=temp_ind, aes(x=as.numeric(interval)+0.15, group=interaction(subj, F0), ymin=threshold_mean-threshold_sd, ymax=threshold_mean+threshold_sd)) +
	# Configure theme, axis text, and panels
	theme_bw() +
	theme(axis.text.y=element_text(size=1*font_scale),   # axis tick label font size
	  axis.text.x=element_text(size=1*font_scale),
	  axis.title.y=element_text(size=1.2*font_scale),    # axis label font size
	  axis.title.x=element_text(size=1.2*font_scale),
	  legend.text=element_text(size=1*font_scale),     # legend text font size
	  legend.title=element_text(size=1.2*font_scale),  # legend title font size
	  strip.text.x=element_text(size=1*font_scale),    # facet label font size
	  plot.title=element_text(size=1.4*font_scale),      # figure title font size
	  panel.grid.major=element_blank(),
	  panel.grid.minor = element_blank()) +
	# Configure scales
	scale_color_manual(values=c("#7fc97f", "#a28ac1"), guide=guide_legend(title="F0 (Hz)")) +
	# Titles and labels
	xlab("Interval Size (multiple of F0DL)") +
	ylab("TMR at threshold")

# Save to disk
ggsave('../plots/fig2.png')