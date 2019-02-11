# QCPA-Matlab (Copyright by Cedric van den Berg)

Written and compiled by Cedric van den Berg. This version started development in Feb. 2016 after a series of adaptations and modifications of Matlab code written by Prof. John Endler following initial coding initiated in Feb. 2015 by Cedric van den Berg during his master thesis at UQ. QCPA is designed to use the output from the MICA toolbox to run Colour Adjacency Analysis & Visual Contrast Analysis (Endler 2012, Endler & Mielke 2005) via a dedicated graphical user interface. It also includes a simple approach to modelling spatial acuity via a median filter kernel as well as spatiochromatic plots in colour space. It has been adapted into JAVA in October 2018 as part of the MICA toolbox (Troscianko & Stevens 2015) and substantially expanded. This version of the QCPA only works with di- and tri-chromatic viewers but can easily be adaptapted to tetra-chromatic vision. Chromaticity and hue are calculated using circular metrics in a maxwell triangle where hue corresponds to the angle and chromaticity (Saturation) to the distance from the achromatic point.

Refs:

Endler, J.A. & Mielke, P.P.W. 2005. Comparing entire colour patterns as birds see them. Biol. J. Linn. Soc. 86: 405–431

Endler, J.A. 2012. A framework for analysing colour pattern geometry: adjacent colours. Biol. J. Linn. Soc. 107: 233–253

Troscianko, J. & Stevens, M. 2015. Image calibration and analysis toolbox – a free software suite for objectively measuring reflectance, colour and pattern. Methods Ecol. Evol. 1320–1331

# Initial Parameters

The user can specify the key parameters for the analysis in the top section of the code:

'target':       This parameter is used for spatial acuity modelling as well as setting the pixels/mm ratio to a common standard. This is                 necessary to achieve the same spatial sampling across all images in a dataset. 

'pixsp':        This parameter describes the sampling density (The space between each horizontal and vertical transect).

'acd':          This parameter gives the spatial acuity in cycles/degree.

'dist':         Specify the viewing distance in mm

'blr':          Set to '1' if spatial acuity should be modelled or to '0' if not. If '1' then QCPA will provide both outputs for the                     'blurred image' as well as the 'unblurred image'

'pure':         This parameter lets the user define the number of times a 3x3 pixel median filter should be run across an image to 
                remove 'salt & pepper' noise.

'cutoff':       This parameter lets the user define a minimum cluster size. Below this threshold the code interpolates members of                       clusters based on their neighbouring pixels.

# Functions of the QCPA

Data import:            The QCPA code is designed to read the input from the multispectral image calibration and analysis toolbox.                               Namely: 1. A clustered .txt image (Zone map) 2. The corresponding cluster ID file 3. A reconstructed RGB file                           (.jpg or .png)

Image orientation:      As pattern analyses are very sensitive to the orientation of a pattern the code lets the user specify the top an                         bottom of an organism along which the image will be orientated.

Image re-scaling:       The parameters of the pattern analyses are very sensitive to the sampling density and as such all images of a                           dataset should be rescaled to a common pixel/mm ratio. The code achieves this by using a user specified size                             standard in the image and a user defined 'target' resolution.

Image segmentation:     Through a series of GUI interfaces the user is able to select animal & background and refine the selections                             through dedicated tools, such as a paintbrush and a magnetic lasso. The code automatically moves the outline of                         the animal 50% of a visual angle inside the animal and 50% of an angle away from the animal. This creates a                             visual acuity and viewing distance specified 'Animal-Background border' or 'Donut'. 

Image interpolation:    The QCPA has the ability to interpolate the 'Donut' based on the neighbouring pixels. This is useful in case of                         shadows cast by the organism. However, this process is very computing intensive and should be avoided if                                 possible.

Pattern analysis:       The QCPA has the ability to calculate adjacency parameters as well as visual contrast parameters of either an                           animal, the animal-background transition zone or a background or an animal vs. a background. Output parameters                           are explained in the output parameter chapter for each analysis.

Data visualisation:     The QCPA provides the user with a plot showing the location of each colour pattern in a maxwell triangle with                           the relative size of the pattern elements indicated by the symbol size.

Output files:           The QCPA provides the user with the following output files:  1. Greyscale zone map of animal, background & donut
                                                                                     2. RGB image of animal, background & donut
                                                                                     3. Maxwell triangle for animal, background & donut
                                                                                     4. Pattern statistics for animal, background, donut                                                                                         and animal vs. background
                                                                                     5. Transition matrices for animal, background &                                                                                             donut
                                                                                     6. Patch specific parameters for animal, background                                                                                         & donut
                                                                                    
# Pattern Statistics (See Endler 2012 and Endler & Mielke 2005 for more details and equations)

ClrCmplx:             'Colour Complexity'. Max=2 When norm_TDiv and norm_CDiv are 1, referring to the maximum possible eveness in a                            colour pattern.

ClrCmplxDiff:         'Colour Complexity Difference'. Absolute difference between animal and background.

rel_CLrCmplxDiff:     'Relative Colour Complexity Difference'. Relative difference between animal and background in %.

Cmplx:                'Complexity'. The proportion of non-synonymous transitions in a transition matrix. Max=1 when every pixel has a                          different neighbour.

CmplxDiff:            'Complexity Difference'. Absolute difference between animal and background.

rel_CmplxDiff:        'Relative Complexity Difference'. Relative difference between animal and background in %.

Asptr:                'Aspect Ratio'. pxUD divided by pxLR. A value greater than 1 is vertically elongated, below 1 is horizontally                            elongated.

AspDiff:              'Aspect Ratio Difference'. The absolute difference between animal and background. 

rel_AspDiff:          'Relative Aspect Ratio Difference'. The realtive difference between animal and background in %.

pxUD:                 'Average Vertical Patch Size'. Average vertical (UD=Up/Down) colour pattern element length (patch size) in pixels.

pxUDDiff:             'Average Vertical Patch Size Difference'. Absolute difference between animal and background. 

rel_PxUDDiff:         'Relative Average Vertical Patch Size Difference'. Relative difference between animal and background in %.

pxLR:                 'Average Horizontal Patch Size'. Average horizontal (LR=Left/Right) colour pattern element length (patch size) in                        pixels.

pxLRDiff:             'Average Horizontal Patch Size Difference'. Absolute difference between animal and background. 

rel_PxLRDiff:         'Relative Average Horizontal Patch Size Difference'. Relative difference between animal and background in %.

AvgPatchDiff_mm2:     'Average Patch Size Difference in mm^2'. Absolute Patch size difference between animal and background.

rel_avgPatchDiff_mm:  'Relative Average Patch Size Difference'. Relative Patch size DIfference in %.

TDiv:                 'Transition Diversity'. Iverse Simpson Diversity Index of the off-diagonals in the transition matrix.

norm_TDiv:            'Normalised Transition Diversity'. Max=1

TdivDiff:             'Transition Diversity Difference' between animal and background.

normTdivDiff:         'Normalised Transition Diversity Difference' between animal and background.

rel_normTdivDiff:     'Relative Normalised Transition Diversity Difference' between animal and background in %.

CDiv:                 'Colour Diversity' Inverse Simpson Diversity Index of the diagonal in the transition matrix.

norm_CDiv:            'Normalised Colour Diversity'. Max=1.

CDivDiff:             'Colour Diversity Difference' between animal and background.

normCDivDiff:         'Normalised Colour Diversity Difference' between animal and background.

rel_normCDivDiff:     'Relative Normalised Colour Diversity Difference' between animal and background in %.

CVLum:                'Coeffient of Luminance Variance' 

CVLumDiff:            'Coefficient of Luminance Variance Difference' between animal and background.

rel_CVLumDiff:        'Relative Coefficient of Luminance Variance Difference' between animal and background in %.

CVChr:                'Coefficient of Chromaticity Variance'

CVChrDiff:            'Coefficient of Chromaticity Variance Difference' between animal and background.

rel_CVChrDiff:        'Relative Coefficient of Chromaticity Variance Difference' between animal and background in %.

CVHue:                'Coefficient of Hue Variance'

CVHueDiff:            'Coefficient of Hue Variance Difference' between animal and background.

rel_CVHueDiff:        'Relative Coefficient of Hue Variance Difference' between animal and background in %.

MnLum:                'Mean Luminance' of a colour pattern using the relative abundance weighted pattern average.

MnLumDiff:            'Mean Luminance Difference' between animal and background.

rel_MnLumDiff:        'Relative Mean Luminance Difference' between animal and background in %.

MnChr:                'Mean Chromaticity' of a colour pattern using the relative abundance weighted pattern average.

MnChrDiff:            'Mean Chromaticity Difference' between animal and background.

rel_MnChrDiff:        'Relative Mean Chromaticity Difference' between animal and background in %.

MnHue:                'Mean Hue' of a colour pattern using the relative abundance weighted pattern average.

MnHueDiff:            'Mean Hue Difference' between an animal and background.

rel_MnHueDiff:        'Relative Mean Hue Difference' between an animal and background in %.

# Example Data

I have attached three seperate examples of input and output data for the user to play around with. All are nudibranch molluscs photographed against their natural backgrounds. The cone catch quanta have been modelled according to a natural 5m depth light spectrum (slightly greenish) and tri-chromatic triggerfish vision with the double cone as the luminance channel using the MICA toolbox. by Jolyon Troscianko. 

# Using QCPA on a 4K screen

QCPA runs into scaling issues on 4k screens. One solution is to change the screen resolution. However, I have been able to fix these issues and a 4K version is available on request.




