source('config.R')
library(RcppCNPy)

# Load simulations
sims = list.files('figure7', pattern='.csv')
f0dls = data.frame()
for (sim in 1:length(sims)) {
	# Import each simulation CSV
	temp = read.csv(file.path('figure7', sims[sim]))
	# If level is numeric, that means it's a phase roving simulation --- change level to str
	if (class(temp$level) == 'numeric') {
		temp$level = as.character(temp$level)
	}
	f0dls = bind_rows(f0dls, temp)
}
f0dls$threshold = f0dls$result
f0dls$nominal_level = factor(f0dls$nominal_level)
f0dls$roving_type = factor(f0dls$roving_type, levels=c("none", "level", "phase"),
						   labels=c("None", "Level Roved", "Phase Randomized"))
f0dls$model = factor(f0dls$model, levels=c("Heinz2001", "Zilany2014", "Verhulst2018"),
					labels=c("Heinz et al. (2001)", "Zilany et al. (2014)", "Verhulst et al. (2018)"))

# Load in vector strength data
freqs = npyLoad('nofigure/vector_strength_curves/freqs.npy')
freqs = as.data.frame(freqs)
colnames(freqs) = c('freq')
vs = data.frame()
for (model in c('Heinz2001', 'Zilany2014', 'Verhulst2018')) {
	temp = npyLoad(file.path('nofigure/vector_strength_curves', paste0(model, '_means.npy')))
	temp = as.data.frame(temp)
	colnames(temp) = c('vs')
	temp = cbind(freqs, temp)
	temp$model = model
	vs = rbind(vs, temp)
}
vs$model = factor(as.factor(vs$model), levels=c('Heinz2001', 'Zilany2014', 'Verhulst2018'), c('Heinz et al. (2001)', 'Zilany et al. (2014)', 'Verhulst et al. (2018)'))

# Append Q10 data to model simulations, interpolating via loess on log-log coordinates where needed
f0dls$vs = 0
for (model in levels(vs$model)) {
	f0dls[f0dls$model == model, ]$vs = 10^approx(log10(vs[vs$model == model, ]$freq), log10(vs[vs$model == model, ]$vs), log10(f0dls[f0dls$model == model, ]$F0*8))[[2]]
}

# Load in Q10 data
freqs = npyLoad('nofigure/tuning_curves/cfs.npy')
freqs = as.data.frame(freqs)
colnames(freqs) = c('freq')
q10 = data.frame()
for (model in c('Heinz2001', 'Zilany2014', 'Verhulst2018')) {
	temp = npyLoad(file.path('nofigure/tuning_curves', paste0(model, '_q10s.npy')))
	temp = as.data.frame(temp)
	colnames(temp) = c('q')
	temp = cbind(freqs, temp)
	temp$model = model
	q10 = rbind(q10, temp)
}
q10$model = factor(as.factor(q10$model), levels=c('Heinz2001', 'Zilany2014', 'Verhulst2018'), c('Heinz et al. (2001)', 'Zilany et al. (2014)', 'Verhulst et al. (2018)'))

# Append Q10 data to model simulations, interpolating via loess on log-linear coordinates where needed
f0dls$q10 = 0
for (model in levels(q10$model)) {
	f0dls[f0dls$model == model, ]$q10 = approx(log10(q10[q10$model == model, ]$freq), q10[q10$model == model, ]$q, log10(f0dls[f0dls$model == model, ]$F0*8))[[2]]
}

# At this point, f0dls should contain, for each test frequency, predicted thresholds for each model and decoding type, Q10 values, and vector strenght values

plot_vs <- function(models=c("Heinz et al. (2001)", "Zilany et al. (2014)", "Verhulst et al. (2018)"), 
			        filename="fig7d_vector_strength_new", width=7, height=3) {
#' Plots a correlation between all-information thresholds and vector strength at the corresponding frequency in each model
#' @param models A vector of model names to plot
#' @param filename The filename (without a file extension) of the plot to be saved in the plots folder
#' @param width Passed to ggsave
#' @param height passed to ggsave

# Filter the data to only include no-roving simulations at 30 dB SPL 
filtered_data = f0dls %>%
	filter(roving_type == 'None') %>%
	filter(nominal_level == 30) %>%
	filter(decoding_type == 'AI') %>%
	filter(model %in% models) 

# Extract the correlation between log-transformed thresholds and log-transformed reciprocal of vector strength
corrs = filtered_data %>% group_by(model) %>% summarize(corr=cor(log10(threshold/(F0)*100), log10(1/(vs*2000))))

# Fit LM between log-transformed thresholds and log-transformed reciprocal of vector strength
filtered_data$rvs = 1/filtered_data$vs
vs_vs_pred_model = lm(log10(threshold/(F0)*100) ~ log10(rvs), data=filtered_data)
beta_0 = coef(vs_vs_pred_model)[1]
beta_1 = coef(vs_vs_pred_model)[2]

# Plot the data
filtered_data %>%
	# Aesthetics calls
	ggplot(aes(x=F0, y=threshold/(F0)*100, shape=decoding_type)) +
	# Geoms
	geom_hline(yintercept=0.004, linetype='dotted', color='gray') +
	geom_vline(xintercept=5000, linetype='dotted', color='gray') + 
	geom_point(aes(y=10^(beta_0 + beta_1*log10(rvs))), color='blue') +
	#geom_smooth(se=FALSE, size=size_smooth) +
	geom_point(size=size_point*2, shape=1) +
	# Axes
	scale_y_log10(breaks=breaks, labels=labels, sec.axis=sec_axis(~ 1/10^((log10(.)-beta_0)/(beta_1)), name="Vector strength")) +
	scale_x_log10(breaks=breaks, labels=labels) +
	# Faceting
	facet_grid(. ~ model) +
	# Theme
	theme_bw() +
	theme(axis.text.y=element_text(size=1*font_scale),   # axis tick label font size
	  axis.text.x=element_text(size=1*font_scale),
	  axis.title.y=element_text(size=1.2*font_scale),    # axis label font size
	  axis.title.x=element_text(size=1.2*font_scale),
	  legend.text=element_text(size=1*font_scale),     # legend text font size
	  legend.title=element_text(size=1.2*font_scale),  # legend title font size
	  strip.text.x=element_text(size=1.2*font_scale),    # facet label font size
	  strip.text.y=element_text(size=1*font_scale),    # facet label font size
	  plot.title=element_text(size=1.5*font_scale),      # figure title font size
	  panel.grid.major=element_blank(),
	  panel.grid.minor = element_blank(),
	  axis.ticks.x=element_line(size=ticksizes),
	  legend.spacing.y=unit(0.05, 'cm'),
	  legend.margin=unit(0, 'cm')) +
	# Labels and legends
	xlab("F0 (Hz)") +
	ylab("F0DL (%)") +
	guides(color=FALSE,
	       shape=FALSE,
	       linetype=guide_legend(title="Decoding Type")) +
	scale_color_manual(values=c('#8dd3c7', '#eded51', '#bebada'))

# Save plot to disk
ggsave(paste0("plots/", filename, ".png"), width=width, height=height)
}


# Define function to plot dual y-axis q10
plot_q10 <- function(models=c("Heinz et al. (2001)", "Zilany et al. (2014)", "Verhulst et al. (2018)"), 
			        filename="fig7d_tuning_new", width=7, height=3) {
#' Plots a correlation between rate-place thresholds and Q10 at the corresponding frequency in each model
#' @param models A vector of model names to plot
#' @param filename The filename (without a file extension) of the plot to be saved in the plots folder
#' @param width Passed to ggsave
#' @param height passed to ggsave

# Start by filtering data and calculating some correlations
filtered_data = f0dls %>%
	filter(roving_type == 'None') %>%
	filter(nominal_level == 30) %>%
	filter(decoding_type == 'RP') %>%
	filter(model %in% models) %>%
	filter(model != 'Verhulst et al. (2018)')  # TODO: remove once I having tuning curves again!

# Extract the correlation between log-transformed thresholds and log-transformed reciprocal of q10
corrs = filtered_data %>% group_by(model) %>% summarize(corr=cor(log10(threshold/(F0)), log10(1/(q10))))

# Fit LM between log-transformed thresholds and log-transformed reciprocal of Q10
filtered_data$rq10 = 1/filtered_data$q10
q10_vs_pred_model = lm(log10(threshold/(F0)*100) ~ log10(rq10), data=filtered_data)
beta_0 = coef(q10_vs_pred_model)[1]
beta_1 = coef(q10_vs_pred_model)[2]

# Plot the data
filtered_data %>%
	# Aesthetics calls
	ggplot(aes(x=F0, y=threshold/(F0)*100, shape=decoding_type)) +
	# Geoms
	geom_hline(yintercept=0.12, linetype='dotted', color='gray') +
	geom_vline(xintercept=4000, linetype='dotted', color='gray') + 
	geom_point(aes(y=10^(beta_0 + beta_1*log10(rq10))), color='blue') +
	#geom_smooth(se=FALSE, size=size_smooth) +
	geom_point(size=size_point*2, shape=0) +
	# Axes
	scale_y_log10(breaks=breaks, labels=labels, sec.axis=sec_axis(~ 1/10^((log10(.)-beta_0)/(beta_1)), name="Q10")) +
	scale_x_log10(breaks=breaks, labels=labels) +
	# Facets
	facet_grid(. ~ model) +
	# Theme
	theme_bw() +
	theme(axis.text.y=element_text(size=1*font_scale),   # axis tick label font size
	  axis.text.x=element_text(size=1*font_scale),
	  axis.title.y=element_text(size=1.2*font_scale),    # axis label font size
	  axis.title.x=element_text(size=1.2*font_scale),
	  legend.text=element_text(size=1*font_scale),     # legend text font size
	  legend.title=element_text(size=1.2*font_scale),  # legend title font size
	  strip.text.x=element_blank(),    # facet label font size
	  strip.background=element_blank(),
	  strip.text.y=element_text(size=1*font_scale),    # facet label font size
	  plot.title=element_text(size=1.5*font_scale),      # figure title font size
	  panel.grid.major=element_blank(),
	  panel.grid.minor = element_blank(),
	  axis.ticks.x=element_line(size=ticksizes),
	  legend.spacing.y=unit(0.05, 'cm'),
	  legend.margin=unit(0, 'cm')) +
	# Labels and legends
	xlab("F0 (Hz)") +
	ylab("F0DL (%)") +
	guides(color=FALSE,
	       shape=FALSE,
	       linetype=guide_legend(title="Decoding Type")) +
	scale_color_manual(values=c('#8dd3c7', '#eded51', '#bebada'))

# Save to disk
ggsave(paste0("plots/", filename, ".png"), width=width, height=height)
}

# Now, call these functions to generate our plots! 
#plot_vs(c("Zilany et al. (2014)"), "fig7d_vector_strength_zilany_only", width=2.5, height=2.0)
#plot_q10(c("Zilany et al. (2014)"), "fig7d_tuning_zilany_only", width=2.5, height=2.0)
plot_q10()
plot_vs()
#plot_correlations(width=3.5, height=3)