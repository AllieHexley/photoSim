% README
% Created by ACH 21/10/2020

To run the photosim toolbox, you need to first simulate your reference database, which you can do by running simulateReferenceDatabase.m first.
You then need to run the metrics to compare the displays to the real-world spectra, which you do by running runMetrics.m
You can visualize the results by running the plotting scripts. The component functions can be found in the functions folder.

Note you must have Psychtoolbox installed and on your path in order to run this toolbox.

Simulations:
- simulateReferenceDatabase: (first script to run) Generate real-world spectra, display spectra, and calculate photoreceptor excitations and MacLeod-Boynton coorindates from real-world spectra
- runMetrics: (second script to run) Calculate "photosim" metrics for given displays

Data:
-99Reflectances.csv: Data on 99 reflectance samples from IES TM-30-2015 used in this analysis
-99Reflectances_sampleTypes.csv: Type of surface for each reflectance in the 99Reflectances database
-401Illuminants.csv: Data on 401 illuminants from the Houser et al., 2013 databse used in this analysis
-401Illuminants_sampleTypes.csv: Type of illuminant for each illuminant in the 401Illuminants database
-lin2012xyz10e_1_7sf.csv: CIE 2015 XYZ functions
-CRT: Radiance spectrum of primaries from NEC CRT Color Monitor MultiSync FP2141SB (NEC Display Solutions, Tokyo, Japan)
-LCD: Radiance spectrum of primaries from Dell U2715H LCD Monitor (Dell, Austin, TX, USA)
-Display++: Radiance spectrum of primaries from CRS Display++ LCD Monitor (Cambridge Research Systems, Rochester, UK)

Functions:
-colorcet: Function to make perceptually unifrom colormaps (Peter Kovesi. Good Colour Maps: How to Design Them.
arXiv:1509.03700 [cs.GR] 2015)
-getSimulatedSpectra: Generate radiance spectra of the 401 illuminants from the 99 surfaces
-getSpectralLocusSpectra: Generate spectral locus over specified wavelength range
-normIllSpd: Normalize illuminant spectra to habe unit area
-XYZToxyY: Convert CIE 2015 XYZ coordinates to CIE xyY chromaticity coordinates (http://psychtoolbox.org/docs/XYZToxyY)
-GetCIES026: Return CIE026 spectral sensitivity functions (https://github.com/spitschan/SilentSubstitutionToolbox/blob/master/PhotoreceptorSensitivity/GetCIES026.m)
-getDistortions: Calculate distortions introduced to each photoreceptor signal when you constrain reproduction on a 3-primary display
-getColourGamut: Return % of real-world chromaticities that can be displayed in gamut for each display
-getPSRM: Return photoreceptor reproduction metric, PSRM, i.e. % of real-world (undistorted) LMSRI quintuplets that can be reproduced by a display (to a specified degree of accuracy)
-getPSDM: Return photoreceptor distortion metric, PSDM, for each display
-getPCDM: Return correlation distortion metric, PCDM, for each display

Plotting scripts:
-Plots figures from corresponding manuscript: https://doi.org/10.1101/2021.02.27.433203; see manuscript for figure descriptions
