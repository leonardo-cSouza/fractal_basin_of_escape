# fractal_basin_of_escape
Code repository accompanying the publication "Fractal and Wada escape basins in the chaotic particle drift motion in tokamaks with electrostatic fluctuations".

This project contains the code to generate and plot the data from all figures.

## Figure 1

To obtain the data and plot the Figure 1, use ``` python perfis.py ```.

## Figure 2

To obtain the data and plot the Figure 2, use ``` python plot_phase_space.py ```. Vary the parameter phi in the code to obtain all the figures (a), (b), (c) and (d).

## Figure 3

To obtain the data used in Figure 3, use ``` gfortran dig_horton.f90 -o exe.x ``` and execute ```exe.x```, and then```python plot_fig3.py ``` to plot the data.

## Figure 4

To obtain the data used in Figure 4, use ``` gfortran basins_horton.f90 -o exe.x ``` and execute ```exe.x```, and then combine the file 1 and 2 for exemple in unix ```type exit-1.dat exit-2.dat > file.dat ``` with file.dat being call in the plot_fig4.py use ``` python plot_fig4.py ``` to plot the data.

## Figure 5

To obtain the data used in Figure 5, use the same procedure to Figure 4, but limiting the region of interest in the ```basins_horton.f90 ``` then combine the two exits and use ``` python plot_fig4.py ``` to plot the data, with the new values of xmax, ymax, xmin and ymin.

## Figure 6

To obtain the data used in Figure 6, use the same data used in Figure 4, namelly file.dat, but now the plotting is different use ``` python plot_fig6.py ```.

## Figure 7

To obtain the data used in Figure 7, use ``` gfortran manifolds.f90 -o exe.x ``` and execute ```exe.x```, and then```python plot_fig7.py ``` to plot the data.

## Figure 8

To obtain the data used in Figure 8....

## Figure 9

To obtain the data used in Figure 9, use ``` gfortran basins_wada.f90 -o exe.x ``` and execute ```exe.x```, and then```python plot_fig9.py ``` to plot the data.
The code ```basins_wada.f90``` is a slight modification of ```basins_horton.f90``` with the addition of a third exit. Moreover, ```plot_fig9.py``` is the same code of ```plot_fig4.py``` plotting also the unstable manifold calulated with ```manifolds.f90```.

## Figure 10

To obtain the data used in Figure 10, use ``` gfortran grid_wada.f90 -o exe.x ``` and execute ```exe.x```, and then```python plot_fig10.py ``` to plot the data.

## Figure 11

To obtain the data used in Figure 11, adding the col term in the manilfods.f90 use ``` gfortran manilfods.f90 -o exe.x ``` the and execute ```exe.x```, and then```python plot_fig11.py ``` to plot the data.




