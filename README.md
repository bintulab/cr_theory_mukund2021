# Temporal signaling, population control, and information processing through chromatin-mediated gene regulation

Code + data associated with Mukund, Bintu 2021: Temporal signaling, population control, and
information processing through chromatin-mediated gene regulation.

## Using this repository

Code to produce the plots for each figure are in the directories corresponding to each figure. Data
used for multiple figures -- namely, the Gillespie algorithm simulations for Fig. 2 and S1; the
incoherent feedforward loop simulations for Fig. 3 and S2; and the information theory flow cytometry
data for Fig. 4 and S3 -- are all in their respective directories. 

For some of these plots, there's some element of randomness involved (i.e. I call some function like
``rand()`` or the like) -- and so if the plots vary very slightly from bit to bit that's likely
why. It's not enough to affect the conclusions we draw, but there's no need to panic. The plots for
Panel S1B can take a while to re-run enough times to get all 9 subpanels looking nice, though. 

If you have any questions please email me -- amukund [at] stanford [dot] edu. If you use, adapt,
modify, etc. any of this code in a paper, please cite our work.

## (Re)-producing the plots
### Figure 1
Just run the Julia script to produce plots for Figure C-F.

### Figure 2
Run the Julia script to produce plots B-E, and the R script to produce plots F-G.

### Figure 3
Run the Julia script to produce plots A-E, and the R script to produce plots F-G.

### Figure 4
Run the Julia script to produce plots C-D, G-H.

### Figure S1
Run the R script to produce plots A-D.

### Figure S2
Panel A: do this in the folder ``fig_s2a``. To produce the dataframe from the plot, run the
``simulate_crs.jl`` script. To use the existing dataframe, merge the ``split_randomized_df`` files
(e.g. using pandas) and then use the resulting merged csv file in the Jupyter notebook
``2020.06.10_cr_paramsweep_plots.ipynb``. This will produce plot A.

Panel B: do this in the folder ``fig_s2b``. Run the Julia script to produce the graded response csv
file, and then run the R script to produce plot B.

Run the Julia script in the main Figure S2 folder to produce plots C-E, and the R script to produce
plots F-H.

### Figure S3
Run the Julia script to produce plots A-C.

## Processing/producing the datasets

## Gillespie simulations
Run the Jupyter notebook (pretty sure you can just run it straight through) to reproduce the
datasets. 

## Incoherent feedforward loop simulations
These plots are produced by scripts in the ``fig_2.jl`` -- the important functions are
``get_hill_dbase()``, ``get_response_dbase_3f()``, and ``get_response_dbase()``. Otherwise, not much
to see here. 

## Information-theoretic analysis of flow cytometry data
Just run the ``full_channel_capacity_pipeline.jl`` file and it should produce the channel capacity
outputs, which are directly used in the code for the ``fig_4.jl`` script. 
