source('config.R')

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

# Construct plot
baseline = fdls[fdls$roving_type == 'None', ]
baseline1 = baseline
baseline2 = baseline
baseline3 = baseline
baseline2$roving_type = 'Level Roved'
baseline3$roving_type = 'Phase Randomized'
fdls_temp1 = rbind(baseline1, baseline2, baseline3)
fdls_temp2 = fdls
fdls_temp1$comparison = 'No Roving'
fdls_temp2$comparison = 'Roved'
fdls_temp = rbind(fdls_temp1, fdls_temp2)
fdls_temp %>%
	# Filter out simulations and roving types to only get what we want
	filter(model == 'Zilany et al. (2014)') %>%
	filter(roving_type != 'None') %>%
	filter(nominal_level == 30) %>%
	# Aesthetics
	ggplot(aes(x=freq, y=threshold/(freq)*100, shape=decoding_type, linetype=comparison)) +
	# Geoms
	geom_vline(xintercept=c(280*8, 1400*8), linetype="dashed", color="gray") +
	geom_smooth(se=FALSE, size=size_smooth) +
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
	xlab("Freq (Hz)") +
	ylab("FDL (%)") +
	guides(color=guide_legend(title="Level per\ncomponent\n(dB re:\nthreshold)"),
	       shape=guide_legend(title="Decoding Type"),
	       linetype=guide_legend(title="Roving Type")) +
	# Line types
	scale_linetype_manual(values=c('dotted', 'solid')) +
	# Facets
	facet_grid(. ~ roving_type)
ggsave('plots/fig6b.png', width=6, height=2.5)
