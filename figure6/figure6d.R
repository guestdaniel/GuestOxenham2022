source('config.R')
library(RcppCNPy)

# Load simulations
sims = list.files('figure6', pattern='.csv')
fdls = data.frame()
for (sim in 1:length(sims)) {
	# Import each simulation CSV
	temp = read.csv(file.path('figure6', sims[sim]))
	# If level is numeric, that means it's a phase roving simulation --- change level to str
	if (class(temp$level) == 'numeric') {
		temp$level = as.character(temp$level)
	}
	fdls = bind_rows(fdls, temp)
}
fdls$threshold = fdls$result
fdls$nominal_level = factor(fdls$nominal_level)
fdls$roving_type = factor(fdls$roving_type, levels=c("none", "level", "phase"),
						   labels=c("None", "Level Roved", "Phase Randomized"))
fdls$model = factor(fdls$model, levels=c("Heinz2001", "Zilany2014", "Verhulst2018"),
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
fdls$vs = 0
for (model in levels(vs$model)) {
	fdls[fdls$model == model, ]$vs = 10^approx(log10(vs[vs$model == model, ]$freq), log10(vs[vs$model == model, ]$vs), log10(fdls[fdls$model == model, ]$freq))[[2]]
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
fdls$q10 = 0
for (model in levels(q10$model)) {
	fdls[fdls$model == model, ]$q10 = approx(log10(q10[q10$model == model, ]$freq), q10[q10$model == model, ]$q, log10(fdls[fdls$model == model, ]$freq))[[2]]
}

# Define function to plot dual y-axis vector strength
plot_vs <- function(models=c("Heinz et al. (2001)", "Zilany et al. (2014)", "Verhulst et al. (2018)"), 
			        filename="fig6d_vector_strength_new", width=7, height=3) {
filtered_data = fdls %>%
	filter(roving_type == 'None') %>%
	filter(nominal_level == 30) %>%
	filter(decoding_type == 'AI') %>%
	filter(model %in% models) 
corrs = filtered_data %>% group_by(model) %>% summarize(corr=cor(log10(threshold/(freq)*100), log10(1/(vs*2000))))

filtered_data %>%
	ggplot(aes(x=freq, y=threshold/(freq)*100, color=model, shape=decoding_type)) +
	geom_point() +
	geom_point(aes(y=1/(vs*2000))) +
	# Geoms
	geom_smooth(se=FALSE, size=size_smooth) +
	geom_point(size=size_point*2) +
	geom_label(data=corrs, aes(label=paste("r = ", round(corr, 2)), x=3000, y=0.07, shape=NULL), color='black') +
	# Annotations
	# Axes
	scale_y_log10(breaks=breaks, labels=labels, sec.axis=sec_axis(~ 1/(2000*.), name="Vector strength")) +
	scale_x_log10(breaks=breaks, labels=labels) +
	facet_grid(. ~ model) +
	# Theme
	theme_bw() +
	theme(axis.text.y=element_text(size=1*font_scale),   # axis tick label font size
	  axis.text.x=element_text(size=1*font_scale),
	  axis.title.y=element_text(size=1.2*font_scale),    # axis label font size
	  axis.title.x=element_text(size=1.2*font_scale),
	  legend.text=element_text(size=1*font_scale),     # legend text font size
	  legend.title=element_text(size=1.2*font_scale),  # legend title font size
	  strip.text.x=element_text(size=1*font_scale),    # facet label font size
	  strip.text.y=element_text(size=1*font_scale),    # facet label font size
	  plot.title=element_text(size=1.5*font_scale),      # figure title font size
	  panel.grid.major=element_blank(),
	  panel.grid.minor = element_blank(),
	  axis.ticks.x=element_line(size=ticksizes),
	  legend.spacing.y=unit(0.05, 'cm'),
	  legend.margin=unit(0, 'cm')) +
	# Labels
	xlab("Vector Strength") +
	ylab("FDL (%)") +
	guides(color=FALSE,
	       shape=FALSE,
	       linetype=guide_legend(title="Decoding Type")) +
	scale_color_manual(values=c('#8dd3c7', '#eded51', '#bebada'))
ggsave(paste0("plots/", filename, ".png"), width=width, height=height)
}

plot_vs()
plot_vs(c("Zilany et al. (2014)"), "fig6d_vector_strength_zilany_only", width=3, height=2)


# Define function to plot dual y-axis q10
plot_q10 <- function(models=c("Heinz et al. (2001)", "Zilany et al. (2014)", "Verhulst et al. (2018)"), 
			        filename="fig6d_tuning_new", width=7, height=3) {
# Start by filtering data and calculating some correlations
filtered_data = fdls %>%
	filter(roving_type == 'None') %>%
	filter(nominal_level == 30) %>%
	filter(decoding_type == 'RP') %>%
	filter(model %in% models)
corrs = filtered_data %>% group_by(model) %>% summarize(corr=cor(log10(threshold/(freq)*100), log10(1/(q10*1))))

filtered_data %>%
	ggplot(aes(x=freq, y=threshold/(freq)*100, color=model, shape=decoding_type)) +
	geom_point() +
	geom_point(aes(y=1/(q10*1))) +
	# Geoms
	geom_smooth(se=FALSE, size=size_smooth) +
	geom_point(size=size_point*2) +
	geom_label(data=corrs, aes(label=paste("r = ", round(corr, 2)), x=3000, y=0.07, shape=NULL), color='black') +
	# Annotations
	# Axes
	scale_y_log10(breaks=breaks, labels=labels, sec.axis=sec_axis(~ 1/(1*.), name="Q10")) +
	scale_x_log10(breaks=breaks, labels=labels) +
	#facet_grid(. ~ model) +
	# Theme
	theme_bw() +
	theme(axis.text.y=element_text(size=1*font_scale),   # axis tick label font size
	  axis.text.x=element_text(size=1*font_scale),
	  axis.title.y=element_text(size=1.2*font_scale),    # axis label font size
	  axis.title.x=element_text(size=1.2*font_scale),
	  legend.text=element_text(size=1*font_scale),     # legend text font size
	  legend.title=element_text(size=1.2*font_scale),  # legend title font size
	  strip.text.x=element_text(size=1*font_scale),    # facet label font size
	  strip.text.y=element_text(size=1*font_scale),    # facet label font size
	  plot.title=element_text(size=1.5*font_scale),      # figure title font size
	  panel.grid.major=element_blank(),
	  panel.grid.minor = element_blank(),
	  axis.ticks.x=element_line(size=ticksizes),
	  legend.spacing.y=unit(0.05, 'cm'),
	  legend.margin=unit(0, 'cm')) +
	# Labels
	xlab("Vector Strength") +
	ylab("FDL (%)") +
	guides(color=FALSE,
	       shape=FALSE,
	       linetype=guide_legend(title="Decoding Type")) +
	scale_color_manual(values=c('#8dd3c7', '#eded51', '#bebada'))
ggsave(paste0("plots/", filename, ".png"), width=width, height=height)
}
plot_q10()
plot_q10(c("Zilany et al. (2014)"), "fig6d_tuning_zilany_only", width=3, height=2)

# Plot summmary figure
plot_correlations <- function(models=c("Heinz et al. (2001)", "Zilany et al. (2014)", "Verhulst et al. (2018)"), 
			        filename="fig6d_correlations", width=4, height=4) {
filtered_data = fdls %>%
	filter(roving_type == 'None') %>%
	filter(nominal_level == 30) %>%
	filter(model %in% models) 

corrs_vs = filtered_data %>% group_by(model, decoding_type) %>% summarize(corr=cor(log10(threshold/(freq)*100), log10(1/(vs*2000))))
corrs_vs$source = "Vector strength"
corrs_q10 = filtered_data %>% group_by(model, decoding_type) %>% summarize(corr=cor(log10(threshold/(freq)*100), log10(1/(q10*1))))
corrs_q10$source = "Q10"
corrs = rbind(corrs_vs, corrs_q10)

corrs %>% ggplot(aes(x=model, y=corr, shape=decoding_type)) + 
	geom_point(size=size_point*2) +
	geom_hline(yintercept=0, linetype="dashed", color="gray", size=1) + 
	# Theme
	theme_bw() +
	theme(axis.text.y=element_text(size=1*font_scale),   # axis tick label font size
	  axis.text.x=element_text(size=1*font_scale, angle=45, hjust=1),
	  axis.title.y=element_text(size=1.2*font_scale),    # axis label font size
	  axis.title.x=element_text(size=1.2*font_scale),
	  legend.text=element_text(size=1*font_scale),     # legend text font size
	  legend.title=element_text(size=1.2*font_scale),  # legend title font size
	  strip.text.x=element_text(size=1*font_scale),    # facet label font size
	  strip.text.y=element_text(size=1*font_scale),    # facet label font size
	  plot.title=element_text(size=1.5*font_scale),      # figure title font size
	  panel.grid.major=element_blank(),
	  panel.grid.minor = element_blank(),
	  axis.ticks.x=element_line(size=ticksizes),
	  legend.spacing.y=unit(0.05, 'cm'),
	  legend.margin=unit(0, 'cm')) +
	# Labels
	xlab("Model") +
	ylab("Correlation (r)") +
	guides(color=FALSE,
	       shape=FALSE,
	       linetype=guide_legend(title="Decoding Type")) +
	facet_grid(decoding_type ~ source)
ggsave(paste0("plots/", filename, ".png"), width=width, height=height)

}

plot_correlations()