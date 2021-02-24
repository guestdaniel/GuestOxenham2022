source('config.R')
library(RcppCNPy)

# Load simulations
sims = list.files(file.path(root_directory, '/figure5/'), pattern='.csv')
f0dls = data.frame()
for (sim in 1:length(sims)) {
	# Import each simulation CSV
	temp = read.csv(file.path(root_directory, '/figure5/', sims[sim]))
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

# Plot F0DLs vs vector strength
f0dls %>%
	filter(roving_type == 'None') %>%
	filter(nominal_level == 30) %>%
	filter(decoding_type == 'AI') %>%
	ggplot(aes(x=vs, y=threshold/(F0)*100, color=model, shape=decoding_type)) +
	# Geoms
	geom_smooth(se=FALSE, size=size_smooth) +
	geom_point(size=size_point*2) +
	geom_vline(xintercept=0.1, linetype='dashed', color='gray') +
	# Axes
	scale_y_log10(breaks=breaks, labels=labels) +
	scale_x_log10(breaks=breaks, labels=labels) +
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
	ylab("F0DL (%)") +
	guides(color=guide_legend(title="Model"),
	       shape=guide_legend(title="Decoding Type"),
	       linetype=guide_legend(title="Decoding Type")) +
	scale_color_manual(values=c('#8dd3c7', '#eded51', '#bebada'))
ggsave('plots/fig5d_phase_locking.png', width=4.25*0.85, height=2.25*0.85)

# Plot F0DLs vs Q10
f0dls %>%
	filter(roving_type == 'None') %>%
	filter(nominal_level == 30) %>%
	filter(decoding_type == 'RP') %>%
	ggplot(aes(x=q10, y=threshold/(F0)*100, color=model, shape=decoding_type)) +
	# Geoms
	geom_smooth(se=FALSE, size=size_smooth) +
	geom_point(size=size_point*2) +
	# Axes
	scale_y_log10(breaks=c(0.07, 0.08, 0.09, 0.1, 0.2, 0.3), labels=c('0.07', '0.08', '0.09', '0.1', '0.2', '0.3'), limits=c(0.07, 0.3)) +
	#scale_x_log10(breaks=breaks, labels=labels) +
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
	xlab("Q10") +
	ylab("FDL (%)") +
	guides(color=guide_legend(title="Model"),
	       shape=guide_legend(title="Decoding Type"),
	       linetype=guide_legend(title="Decoding Type")) +
	scale_color_manual(values=c('#8dd3c7', '#eded51', '#bebada'))
ggsave('plots/fig5d_tuning.png', width=4.25*0.85, height=2.25*0.85)
