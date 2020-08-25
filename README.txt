% README

Simulations:
- simulatePhotoreceptorSignals: Simulate photoreceptor signals from the displays and the real world (should be the first script ran) TO DO: make this 3 primary displays only
- simulatePhotoreceptorDistortions: Calculate distortions introduced to photoreceptor signals reproduced on displays
 after you have transformed your five primary signal space to a 3 primary representation on displays (should be the second script ran) TO DO for publishing: make this LMS only

Data:
-99Reflectances.csv: Data on 99 reflectance samples from IES ... used in this analysis
-99Reflectances_sampleTypes.csv: Type of surface for each reflectance in the 99Reflectances database
-401Illuminants.csv: Data on 401 illuminants from the Houser et al., databse used in this analysis
-401Illuminants_sampleTypes.csv: Type of illuminant for each illuminant in the 401Illuminants database
-daylightSpectra.csv: daylight spectra from ?K to ?K used in this analysis
-lin2012xyz10e_1_7sf.csv: CIE 2015 XYZ functions
-CRT: Radiance spectra of primaries from ??? CRT display
-LCD: Radiance spectra of primaries from ??? LCD display
-Display++: Radiance spectra of primaries from ??? Display++ display

Functions:
-colorcet: Function to make perceptually unifrom colormaps (CITE)
-getChromaticityReproductionMetric: Return % of real-world chromaticities that can be display in gamut for each display
-GetCIES026: Return CIE026 spectral sensitivity functions (CITE)
-getDaylightSpectra: Loads in daylight spectra
-getDistortionReproductionMetric: Return % of distorted LMSRI quintuplets that can be reproduced by a display
-getPhotoreceptorCorrelationDistortionMetric: Return correlation distortion metric for each display
-getPhotoreceptorSignalDistortionMetric: Return photoreceptor distortion metric for each display
-getPhotoreceptorSignalDistortions: Calculate distortions introduced to each photoreceptor signal when you constrain reproduction on a 3-primary display
-getPSDAcrossChromaticities: Find mean photoreceptor distortion metric for all of the real-world spectra that fall within each gridpoint in chromaticity coordinates
-getRealWorldReproductionMetric: Return % of real-world (undistorted) LMSRI quintuplets that can be reproduced by a display (to a specified degree of accuracy)
-getSimulatedSpectra: Generate radiance spectra of the 401 illuminants from the 99 surfaces
-getSpectralLocusSpectra: Generate spectral locus over specified wavelength range
-normIllSpd: Normalize illuminant spectra to habe unit area
-XYZToxyY: Convert CIE 2015 XYZ coordinates to CIE xyY chromaticity coordinates (CITE)

Plotting scripts:
- plotChromaticityReproduction: plot real-world chromaticities on CIE xy diagram along with CIE xy gamut of displays
- visualizeChromaticityReproduction: visualize the MacBeth colour checker chart chromaticities after distorting photoreceptor signals to visualize chromatic distortions that may have been introduced - still TO DO
- plotPhotoreceptorSignalDistortions: plot reproduced photoreceptor signal vs real-world photoreceptor signals 
- plotPSDMHistograms: plot hisotgrams showing the spread of the photoreceptor signal distortion metric
- plotPSDMAcrossChromaticity: plot photoreceptor signal distortion metric across chromaticity diagram
- plotPhotoreceptorSignalDistortionsMB: plot real-world and distorted photoreceptor signals in MacLeod-Boynton space
- plotMetrics: plot bar graphs showing metrics for each display
- plotPhotoreceptorCorrelationDistortions: plot real-world photoreceptor correlations and photoreceptor correlation distortion metric - still TO DO
- plotFundamentals: plot spectral sensitivities, primaries for each display, and XYZ - TO DO
TO ADD: something that looks at illuminant and reflectance type
