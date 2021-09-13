source('config.R')

# Define a function to implement predictions from Micheyl, Xiao, and Oxenham (2012)
predict_micheyl = function(freq, dur, level, gamma_f=0.82, gamma_d=-0.42, gamma_s=-1.09, beta_f=0.38, beta_d=0.42,
						   beta_s=0.37, alpha=-0.38) {
	# Args:
	#	freq (numeric): frequency in Hz
	#   dur (numeriuc): duration in ms
	#   level (numeric): sensation level in dB
	beta_f * (freq / 1000) ^ gamma_f + beta_d * (dur / 100) ^ gamma_d + beta_s * (level / 10) ^ gamma_s + alpha
}

# Predict data from Micheyl at al. model
f = 8*10^seq(log10(280) - 0.2, log10(1400) + 0.1, length.out=100)
model_data = data.frame(freq=f, threshold=10^(predict_micheyl(f, 200, 25))/f*100)

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
fdls %>% 
	filter(nominal_level %in% c(20, 30, 40)) %>%
	filter(roving_type == 'None') %>%
	ggplot(aes(x=freq, y=threshold/(freq)*100)) +
	# Geoms
	geom_vline(xintercept=c(280*8, 1400*8), linetype="dashed", color="gray") +
	geom_smooth(aes(x=freq, y=threshold*2e-3), se=FALSE, size=size_smooth, data=model_data, color='black', shape=1) +
	geom_smooth(aes(color=nominal_level, shape=decoding_type), se=FALSE, size=size_smooth) +
	geom_point(aes(color=nominal_level, shape=decoding_type), size=size_point) +
	# Axes
	scale_y_log10(breaks=breaks, labels=labels, sec.axis=sec_axis(~ . * 5e2, breaks=breaks, labels=labels)) +
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
	       linetype=guide_legend(title="Decoding Type")) +
	# Facets
	facet_grid(. ~ model)
# Save plot to disk
ggsave('plots/fig6a.png', width=8, height=3.5)
