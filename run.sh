conda activate GuestOxenham2021_code

echo 'Calculating absolute thresholds'
python3 nofigure/absolute_thresholds/absolute_thresholds.py

echo 'Calculating tuning_curves'
python3 nofigure/tuning_curves/estimate_freq_level_functions.py
python3 nofigure/tuning_curves/extract_tuning_curves.py
python3 nofigure/tuning_curves/estimate_q10.py
python3 nofigure/tuning_curves/estimate_q10_bm_clicks.py

echo 'Calculating vector strength'
python3 nofigure/vector_strength_curves/vector_strength_curves.py

echo 'Generating Figure 3'
python3 figure3/figure3.py

echo 'Generating Figure 4'
python3 figure4/figure4a.py
Rscript figure4/figure4a.R
python3 figure4/figure4b.py
Rscript figure4/figure4b.R
Rscript figure4/figure4c.R
Rscript figure4/figure4d.R


echo 'Generating Figure 5'
python3 figure5/figure5a.py
Rscript figure5/figure5a.R
python3 figure5/figure5b.py
Rscript figure5/figure5b.R
Rscript figure5/figure5c.R
Rscript figure5/figure5d.R


echo 'Generating Figure 6'
python3 figure6/figure6.py

echo 'Generating Figure 7 (supplemental figure 1)'
python3 supplemental_figure_1/figure7a.py
python3 supplemental_figure_1/figure7b.py
python3 supplemental_figure_1/figure7c.py
python3 supplemental_figure_1/figure7d/generate_ISIs.py
python3 supplemental_figure_1/figure7d/figure7d.py
python3 supplemental_figure_1/figure7e/generate_ISI_histograms.py
python3 supplemental_figure_1/figure7e/figure7e.py