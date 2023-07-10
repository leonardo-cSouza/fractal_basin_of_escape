# fractal_basin_of_escape
Code repository accompanying the publication "Fractal and Wada escape basins in the chaotic particle drift motion in tokamaks with electrostatic fluctuations".

This project contains the code to generate and plot the data from all figures.

## Figure 1

To obtain the data and plot the Figure 1, use ``` python perfis.py ```.

## Figure 2

To obtain the data and plot the Figure 2, use ``` python plot_phase_space.py ```. Vary the parameter phi in the code to obtain all the figures (a), (b), (c) and (d).

## Figure 3

To obtain the data used in Figure 3, use ``` gfortran dig_horton.f90 -o exe.x ``` and execute ```exe.x```, and then ``` python plot_fig3.py ``` to plot the data.

## Figure 4

To obtain the data used in Figure 4, use ``` gfortran basins_horton.f90 -o exe.x ``` and execute ```exe.x```, and then ``` python plot_fig4.py ``` to plot the data.

## Figure 5

To obtain the data used in Figure 5, use the same procedure to Figure 4, but limiting the region of interest in the ```basins_horton.f90 ``` use ``` python plot_fig4.py ``` to plot the data, with the new values of xmax, ymax, xmin and ymin.

## Figure 6

To obtain the data used in Figure 6, use the same data used in Figure 4, but now the plotting is different use ``` python plot_fig6.py ```.



