source('config.R')
library(tidyr)

# Load simulations
sims = list.files(file.path(root_directory, '/supfigure_maskers/masker_interval_simulations/'), pattern='.csv')
f0dls = data.frame()
for (sim in 1:length(sims)) {
	# Import each simulation CSV
	temp = read.csv(file.path(root_directory, '/supfigure_maskers/masker_interval_simulations/', sims[sim]))
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
f0dls$masker_interval = factor(f0dls$masker_interval)

# Set plotting parameters
n_breaks = 10
breaks = c(seq(2, 10, 10/n_breaks) %o% 10^(-8:3))
labels = as.character(breaks)
labels[!(log10(breaks)%%1==0) & breaks != 2000] = ''
ticksizes = rep(.25, length(breaks))
ticksizes[log10(breaks)%%1==0] = 1

# Construct plot #1 (just single example)
f0dls %>% 
	filter(nominal_level %in% c(30)) %>%
	filter(delta %in% c(1e-6, 1e-4, 1e-2)) %>%
	filter(roving_type == 'None') %>%
	pivot_wider(values_from=threshold, names_from=stimulus, id_cols=c(masker_interval, decoding_type, nominal_F0, delta)) %>%
	mutate(ratio=GEOM/ISO) %>%
	ggplot(aes(x=masker_interval, y=ratio, shape=decoding_type, color=nominal_F0)) + 
	# Geoms
	#geom_vline(xintercept=c(280, 1400), linetype="dashed", color="gray") + 
	geom_smooth(se=FALSE, size=size_smooth) + 
	geom_point(size=size_point) + 
	geom_line(aes(group=interaction(nominal_F0, decoding_type))) + 
	# Axes
	#scale_y_log10(breaks=breaks, labels=labels) + 
	#scale_x_log10(breaks=breaks, labels=labels) + 
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
	xlab("Masker interval (ST)") + 
	ylab("F0DL (%)") + 
	guides(color=guide_legend(title="F0 (Hz)"),
	       shape=guide_legend(title="Decoding type")) +
	#scale_color_manual(values=c("#7fc97f", "#a28ac1"), guide=guide_legend(title="nominal_F0 (Hz)")) +
	facet_grid(. ~ delta)
	# Facets
ggsave(file.path('plots', 'supfigure_maskers_2b.png'), width=8, height=3)
