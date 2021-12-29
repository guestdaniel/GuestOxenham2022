cd plots
for file in fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8
do
    convert $file.png $file.tiff
done