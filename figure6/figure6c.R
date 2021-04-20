source('config.R')

# Load simulations
sims = list.files('', pattern='.csv')
fdls = data.frame()
for (sim in 1:length(sims)) {
	# Import each simulation CSV
	temp = read.csv(file.path('', sims[sim]))
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

# Load prior discrimination data
load('data/FDL_data.RData')
pd_fdl = data

# Compute ratio between 8.5 kHz and 2 kHz for the model data
temp_comp = fdls %>% mutate(threshold=threshold/freq*100) %>%
	group_by(roving_type, decoding_type, nominal_level, model) %>%
	summarize(high=10^predict(loess(log10(threshold) ~ log10(freq)), log10(8500)),
			  low=10^predict(loess(log10(threshold) ~ log10(freq)), log10(2000)),
			  ratio=10^predict(loess(log10(threshold) ~ log10(freq)), log10(8500))/10^predict(loess(log10(threshold) ~ log10(freq)), log10(2000))) %>%
	filter(roving_type=='None') %>%
	filter(nominal_level==30)

# Compute ratio between 8.5 kHz and 2 kHz for the behavioral data
temp_behavior = pd_fdl %>%
		group_by(src) %>%
		summarize(high=10^approx(log10(f), log10(t), log10(8500))$y,
				  low=10^approx(log10(f), log10(t), log10(2000))$y,
				  ratio=10^approx(log10(f), log10(t), log10(8500))$y/10^approx(log10(f), log10(t), log10(2000))$y)
# Add additional info to the behavioral dataframe
temp_behavior$model = temp_behavior$src
temp_behavior$decoding_type = 'Behavior'

# Bind together the model data and behavioral data and re-label factors
temp = bind_rows(temp_comp, temp_behavior)
temp$decoding_type = factor(temp$decoding_type, levels=c('Behavior', 'AI', 'RP'), labels=c('Behavior', 'All-information', 'Rate-place'))
temp$model = factor(temp$model,
					levels=c('Heinz et al. (2001)', 'Zilany et al. (2014)', 'Verhulst et al. (2018)', 'Moore1973', 'Moore2012', 'Lau2017', 'Gockel2020'),
					labels=c('Heinz et al. (2001)', 'Zilany et al. (2014)', 'Verhulst et al. (2018)', 'Moore (1973)', 'Moore and Ernst (2012)', 'Lau et al. (2017)', 'Gockel et al. (2020)'))

# Play some tricks to space out the various data points along the x-axis
temp$mod_num = as.numeric(temp$model)
temp[temp$mod_num > 3, ]$mod_num = temp[temp$mod_num > 3, ]$mod_num + 1  # move behavior numbers right by one

# Construct plot
temp %>% ggplot(aes(x=mod_num, y=ratio, shape=decoding_type)) +
		# Geoms
		geom_hline(yintercept=1) +
		geom_point(size=size_point*1.5) +
		# Annotate
		annotate('rect', xmin=0.5, xmax=3.5, ymin=0.5, ymax=90, alpha=0, linetype='dashed', size=1, color='#b3cde3') +
		annotate('label', x=1.95, y=90, size=3, label.size=1, label='AN', color='#b3cde3') +
		annotate('rect', xmin=4.5, xmax=8.5, ymin=1.8, ymax=40, alpha=0, linetype='dashed', size=1, color='#ccebc5') +
		annotate('label', x=6.5, y=40, size=3, label.size=1, label='Behavior', color='#ccebc5') +
		# Axes
		scale_y_log10(breaks=breaks, labels=labels, limits=c(0.4, 100)) +
		scale_x_continuous(breaks=c(1, 2, 3, 5, 6, 7, 8),
						   labels=c('Heinz et al. (2001)', 'Zilany et al. (2014)', 'Verhulst et al. (2018)', 'Moore (1973)', 'Moore and Ernst (2012)', 'Lau et al. (2017)', 'Gockel et al. (2020)')) +
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
		ylab("Ratio (8.5 kHz / 2.0 kHz)") +
		guides(shape=guide_legend(title="Type"))
ggsave('plots/fig6c.png', width=4, height=2.75)