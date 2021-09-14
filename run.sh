# Replicates every figure in Guest and Oxenham (2021), including reproducing underlying simulations

# Get the behavioral data before we start running scripts
mkdir data
wget https://zenodo.org/record/4750384/files/data_archive.zip?download=1 -O data/data_archive.zip
unzip data/data_archive.zip
cp data_archive/* data
rm -r data_archive
rm data/data_archive.zip

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
echo 'Generating Figure 1'
python3 figure0/figure1.py

# Code check:   4/20/2021
# Output check: 4/20/2021
echo 'Generating Figure 2'
Rscript figure2/figure2.R

# Code check:   4/19/2021
# Output check: 4/19/2021
echo 'Generating Figure 3'
python3 figure3/figure3.py

# Code check:   4/20/2021
# Output check: 4/20/2021
echo 'Generating Figure 4'
Rscript figure4/figure4.R

# Code check:   5/11/2021
# Output check: 5/11/2021
echo 'Generating Figure 5'
python3 figure5/figure5a.py
python3 figure5/figure5b.py

# Code check:   4/19/2021
# Output check: 4/19/2021
echo 'Generating Figure 6'
python3 figure6/figure6a.py
Rscript figure6/figure6a.R
python3 figure6/figure6b.py
Rscript figure6/figure6b.R
Rscript figure6/figure6c.R
Rscript figure6/figure6d.R

# Code check:   4/19/2021
# Output check: 4/19/2021
echo 'Generating Figure 7'
python3 figure7/figure7a.py
Rscript figure7/figure7a.R
python3 figure7/figure7b.py
Rscript figure7/figure7b.R
Rscript figure7/figure7c.R
Rscript figure7/figure7d.R

# Code check:   4/19/2021
# Output check: 4/19/2021
echo 'Generating Supplemental Figure 1'
python3 supfigure1/supfigure1a.py
python3 supfigure1/supfigure1b.py
python3 supfigure1/supfigure1c.py
python3 supfigure1/supfigure1d/generate_ISIs.py
python3 supfigure1/supfigure1d/supfigure1d.py
python3 supfigure1/supfigure1e/generate_ISI_histograms.py
python3 supfigure1/supfigure1e/supfigure1e.py