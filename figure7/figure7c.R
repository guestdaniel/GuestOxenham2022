source('config.R')

# Load simulations
sims = list.files(file.path(root_directory, '/figure7/'), pattern='.csv')
f0dls = data.frame()
for (sim in 1:length(sims)) {
	# Import each simulation CSV
	temp = read.csv(file.path(root_directory, '/figure7/', sims[sim]))
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

# Load prior discrimination data
load('data/F0DL_data.RData')
pd_f0dl = data

# Compute ratio between 8.5 kHz and 2 kHz for the model data
temp_comp = f0dls %>% mutate(threshold=threshold/F0*100) %>%
	group_by(roving_type, decoding_type, nominal_level, model) %>%
	summarize(ratio=10^predict(loess(log10(threshold) ~ log10(F0)), log10(1400))/
                    10^predict(loess(log10(threshold) ~ log10(F0)), log10(280))) %>%
	filter(roving_type=='None') %>%
	filter(nominal_level==30)

# Compute ratio between 8.5 kHz and 2 kHz for the behavioral data
temp_behavior = pd_f0dl %>%
		group_by(src) %>%
		summarize(ratio=10^approx(log10(f), log10(t), log10(1400))$y/
                        10^approx(log10(f), log10(t), log10(280))$y)
# Add additional info to the behavioral dataframe
temp_behavior$model = temp_behavior$src
temp_behavior$decoding_type = 'Behavior'

# Bind together the model data and behavioral data and re-label factors
temp = bind_rows(temp_comp, temp_behavior)
temp$decoding_type = factor(temp$decoding_type, levels=c('Behavior', 'AI', 'RP'), labels=c('Behavior', 'All-information', 'Rate-place'))
temp$model = factor(temp$model,
					levels=c('Heinz et al. (2001)', 'Zilany et al. (2014)', 'Verhulst et al. (2018)', 'Lau2017', 'Gockel2018', 'Gockel2020', 'Guest2020'),
					labels=c('Heinz et al. (2001)', 'Zilany et al. (2014)', 'Verhulst et al. (2018)', 'Lau et al. (2017)', 'Gockel and Carlyon (2020)', 'Gockel et al. (2020)', 'Present data'))

# Play some tricks to space out the various data points along the x-axis
temp$mod_num = as.numeric(temp$model)
temp[temp$mod_num > 3, ]$mod_num = temp[temp$mod_num > 3, ]$mod_num + 1  # move behavior numbers right by one

# Construct plot
temp %>% ggplot(aes(x=mod_num, y=ratio, shape=decoding_type)) +
		geom_hline(yintercept=1) +
		geom_point(size=size_point*1.5) +
		# Annotate
		annotate('rect', xmin=0.5, xmax=3.5, ymin=0.5, ymax=15, alpha=0, linetype='dashed', size=1, color='#b3cde3') +
		annotate('label', x=1.95, y=15, size=3, label.size=1, label='AN', color='#b3cde3') +
		annotate('rect', xmin=4.5, xmax=8.5, ymin=2.5, ymax=9, alpha=0, linetype='dashed', size=1, color='#ccebc5') +
		annotate('label', x=6.5, y=9, size=3, label.size=1, label='Behavior', color='#ccebc5') +
		# Axes
		scale_y_log10(breaks=breaks, labels=labels, limits=c(0.4, 20)) +
		scale_x_continuous(breaks=c(1, 2, 3, 5, 6, 7, 8),
						   labels=c('Heinz et al. (2001)', 'Zilany et al. (2014)', 'Verhulst et al. (2018)', 'Lau et al. (2017)', 'Gockel and Carlyon (2020)', 'Gockel et al. (2020)', 'Present data')) +
		# Theme
		theme_bw() +
		theme(axis.text.y=element_text(size=1*font_scale),   # axis tick label font size
		  axis.text.x=element_text(size=0.7*font_scale, angle=45, hjust=1),
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
		xlab("Source (model or paper)") +
		ylab("Ratio (1.4 kHz / 0.28 kHz)") +
		guides(shape=guide_legend(title="Type"))
ggsave('plots/fig7c.png', width=4, height=2.75)