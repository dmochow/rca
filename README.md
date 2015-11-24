# rca
This toolbox provides MATLAB source code to perform “Reliable Components Analysis”, a technique to reduce the dimensionality and increase the interpretability of your EEG/MEG/fMRI data by projecting multivariate data into a few “components”.  Each of these components is a linear combination of the individual sensors, and is characterized by both a spatial topography and associated time course.  While your data set may span hundreds of electrodes, sensors, or voxels, the _reliable_ portion of your data is often captured by just a few of these components.  This makes it easier to evaluate effects and test hypotheses.

Importantly, the criterion used to form the components is physiologically plausible: the components are formed to maximize the trial-to-trial covariance of the data.  Anything that is not consistent from trial-to-trial is generally (but not always) “noise”.  

There are two basic ways of using the code: 

(1) run the data on a samples-by-channels-by-subjects data “cube” where the number of samples is generally much larger than the number of channels or subjects.  This data often comes from recordings of brain activity during long, naturalistic stimuli such as movies.

(2) run the data on a samples-by-channels-by-trials data “cube” where the number of trials is often >100.  This scenario corresponds to the conventional evoked-response design where short trials are repeated many times for a given subject/condition.  

The file rcaDemo.m illustrates a sample usage.

Please cite at least the following paper when using the toolbox:

Dmochowski, J. P., Sajda, P., Dias, J., & Parra, L. C. (2012). Correlated components of ongoing EEG point to emotionally laden attention–a possible marker of engagement?. Frontiers in human neuroscience, 6.

Other papers that have employed the technique:

Dmochowski, J. P., Bezdek, M. A., Abelson, B. P., Johnson, J. S., Schumacher, E. H., & Parra, L. C. (2014). Audience preferences are predicted by temporal reliability of neural processing. Nature communications, 5.

Dmochowski, J. P., Greaves, A. S., & Norcia, A. M. (2015). Maximally reliable spatial filtering of steady state visual evoked potentials. NeuroImage, 109, 63-72.

Dmochowski J. P., Norcia A. M. (2015) Cortical Components of Reaction-Time during Perceptual Decisions in Humans. PLoS ONE 10(11): e0143339. 

(c) Jacek P. Dmochowski, 2015
jdmochowski@ccny.cuny.edu
