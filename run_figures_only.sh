# Replicates every figure in Guest and Oxenham (2021), does not re-run simulations

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

# Code check:   4/19/2021
# Output check: 4/19/2021
echo 'Generating Figure 5'
python3 figure5/figure5.py

# Code check:   4/19/2021
# Output check: 4/19/2021
echo 'Generating Figure 6'
Rscript figure6/figure6a.R
Rscript figure6/figure6b.R
Rscript figure6/figure6c.R
Rscript figure6/figure6d.R

# Code check:   4/19/2021
# Output check: 4/19/2021
echo 'Generating Figure 7'
Rscript figure7/figure7a.R
Rscript figure7/figure7b.R
Rscript figure7/figure7c.R
Rscript figure7/figure7d.R

# Code check:   4/19/2021
# Output check: 4/19/2021
echo 'Generating Supplemental Figure 1'
python3 supfigure1/supfigure1a.py
python3 supfigure1/supfigure1b.py
python3 supfigure1/supfigure1c.py
python3 supfigure1/supfigure1d/supfigure1d.py
python3 supfigure1/supfigure1e/supfigure1e.py