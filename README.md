# Generating-synthetic-ground-motions-reaching-target-spectrum-with-the-optimization-approach
MATLAB code for Generating synthetic ground motions reaching target spectrum with the optimization approach

Author: Pouya Tavakoli
PhD student in structural engineering at McGill University
Pouya.tavakoli@mail.mcgill.ca

To enhance readers' understanding of the concept of generating artificial earthquakes, this section presents the MATLAB code used in the article titled "Generating Synthetic Ground Motions Reaching Target Spectrum with Optimization Approach."
The MATLAB files comprise three distinct components: CMAES (the optimization algorithm), ObjFF (the objective function representing the difference between the acceleration response spectrum of the generated earthquake and the target spectrum), and Plott (utilized for plotting information related to the generated earthquake).
Initially, the CMAES and ObjFF files are executed to generate the earthquake (the CMAES file has to be run first). Subsequently, in the second stage, the output from the preceding stage is imported into the Plott file to extract detailed information about the generated earthquake.
