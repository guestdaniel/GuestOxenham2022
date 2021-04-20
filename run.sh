conda activate GuestOxenham2021_code

# Code check:   4/14/2021
# Output check: 4/14/2021
echo 'Calculating absolute thresholds'
python3 nofigure/absolute_thresholds/absolute_thresholds.py

# Code check:   4/14/2021
# Output check: 4/14/2021
echo 'Calculating tuning_curves'
python3 nofigure/tuning_curves/estimate_freq_level_functions.py
python3 nofigure/tuning_curves/extract_tuning_curves.py
python3 nofigure/tuning_curves/estimate_q10.py
python3 nofigure/tuning_curves/estimate_q10_bm_clicks.py

# Code check:   4/14/2021
# Output check: 4/14/2021
echo 'Calculating vector strength'
python3 nofigure/vector_strength_curves/vector_strength_curves.py

# Code check:   4/14/2021
# Output check: 4/14/2021
echo 'Generating Figure 0'
python3 figure0/figure0.py

# Code check:   4/19/2021
# Output check: 4/19/2021
echo 'Generating Figure 3'
python3 figure3/figure3.py

# Code check:   4/19/2021
# Output check: 4/19/2021
echo 'Generating Figure 3.5'
python3 figure3pt5/figure3pt5.py

# Code check:   4/19/2021
# Output check: 4/19/2021
echo 'Generating Figure 6'
python3 figure4/figure4a.py
Rscript figure4/figure4a.R
python3 figure4/figure4b.py
Rscript figure4/figure4b.R
Rscript figure4/figure4c.R
Rscript figure4/figure4d.R

# Code check:   4/19/2021
# Output check: 4/19/2021
echo 'Generating Figure 6'
python3 figure7/figure7a.py
Rscript figure7/figure7a.R
python3 figure7/figure7b.py
Rscript figure7/figure7b.R
Rscript figure7/figure7c.R
Rscript figure7/figure7d.R

# Code check:   4/19/2021
# Output check: 4/19/2021
echo 'Generating Figure 7 (supplemental figure 1)'
python3 supplemental_figure_1/figure7a.py
python3 supplemental_figure_1/figure7b.py
python3 supplemental_figure_1/figure7c.py
python3 supplemental_figure_1/figure7d/generate_ISIs.py
python3 supplemental_figure_1/figure7d/figure7d.py
python3 supplemental_figure_1/figure7e/generate_ISI_histograms.py
python3 supplemental_figure_1/figure7e/figure7e.py