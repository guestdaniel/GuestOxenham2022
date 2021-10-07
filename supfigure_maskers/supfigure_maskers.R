source('config.R')

# Load simulations
sims = list.files(file.path(root_directory, '/supfigure_maskers/'), pattern='.csv')
f0dls = data.frame()
for (sim in 1:length(sims)) {
	# Import each simulation CSV
	temp = read.csv(file.path(root_directory, '/supfigure_maskers/', sims[sim]))
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
f0dls$stimulus = factor(f0dls$stimulus, levels=c("iso", "geom"), labels=c("ISO", "GEOM"))
f0dls$nominal_F0 = factor(f0dls$nominal_F0, levels=c(280, 1400), labels=c(280, 1400))

# Construct plot
f0dls %>% 
	filter(nominal_level %in% c(30)) %>%
	filter(roving_type == 'None') %>%
	ggplot(aes(x=stimulus, y=threshold/(as.numeric(as.character(nominal_F0)))*100, color=nominal_F0, shape=nominal_level)) + 
	# Geoms
	#geom_vline(xintercept=c(280, 1400), linetype="dashed", color="gray") + 
	geom_smooth(se=FALSE, size=size_smooth) + 
	geom_point(size=size_point) + 
	geom_line(aes(group=interaction(nominal_level, nominal_F0))) + 
	# Axes
	scale_y_log10(breaks=breaks, labels=labels) + 
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
	xlab("F0 (Hz)") + 
	ylab("F0DL (%)") + 
	guides(color=guide_legend(title="Level per\ncomponent\n(dB re:\nthreshold)"),
	       shape=guide_legend(title="Decoding Type"),
	       linetype=guide_legend(title="Decoding Type")) +
	scale_color_manual(values=c("#7fc97f", "#a28ac1"), guide=guide_legend(title="nominal_F0 (Hz)")) +
	# Facets
	facet_grid(decoding_type ~ model)
ggsave(file.path('plots', 'supfigure_maskers_a.png'), width=8, height=3)
