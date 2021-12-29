source('config.R')

# Load simulations
sims = list.files(file.path(root_directory, 'nofigure/parity/'), pattern='.csv')
f0dls = data.frame()
for (sim in 1:length(sims)) {
	# Import each simulation CSV
	temp = read.csv(file.path(root_directory, 'nofigure/parity/', sims[sim]))
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

# Define function to take two thresholds and add them optimally
comb <- function(t1, t2) {
	1/sqrt(1/t1^2 + 1/t2^2)
}

# Create combined information 
temp = f0dls[f0dls$parity == "both", ] 
temp$threshold = comb(f0dls[f0dls$parity == "odd", "threshold"], f0dls[f0dls$parity == "even", "threshold"])
temp$parity = "Combined"
f0dls = rbind(f0dls, temp)

# Modify both thresholds by appropriate factor of 2
f0dls[f0dls$parity == "both", "threshold"] = f0dls[f0dls$parity == "both", "threshold"] * 1/sqrt(2)

# Set parity factor
f0dls$parity = factor(f0dls$parity, levels=c('both', 'odd', 'even', 'Combined'), labels=c('Normal', 'Odd Only', 'Even Only', 'Odd + Even'))

# Construct plot
f0dls %>% 
	filter(roving_type == 'None') %>%
	filter(decoding_type == "AI") %>%
	ggplot(aes(x=F0, y=threshold/(F0)*100, color=parity)) + 
	# Geoms
	geom_vline(xintercept=c(280, 1400), linetype="dashed", color="gray") + 
	geom_smooth(se=FALSE, size=size_smooth) + 
	geom_point(size=size_point) + 
	# Axes
	scale_y_log10(breaks=breaks, labels=labels, limits=c(1e-4, 1e-1)) + 
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
	# Scales
	scale_color_manual(values=c('black', 'red', 'blue', 'darkgray')) +
	# Labels
	xlab("F0 (Hz)") + 
	ylab("F0DL (%)") + 
	guides(color=guide_legend(title="Parity"),
	       shape=guide_legend(title="Decoding Type"),
	       linetype=guide_legend(title="Decoding Type")) +
	# Facets
	facet_grid(. ~ nominal_level)
ggsave(file.path('plots', 'parity.png'), width=8, height=3)
