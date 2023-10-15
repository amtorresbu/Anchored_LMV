# Anchored_LMV

This code in this repository is used to simulate triangulation of one line correspondence and a number of point correspondences across a fixed number of images in the julia language using homotopy continuation methods. 

-- Start-Systems.jl initializes start systems for homotopy continuation for each of the reconstruction approaches.

-- Reconstruction-Fcns.jl defines all julia functions that simulate triangulation for each approach.

-- Figures-Tables.jl creates all figures and values for all tables displayed in the paper.

-- EDD-Fcns.jl contain the functions that numerically calculate Euclidean distance degrees for the anchored multiview varieties.
