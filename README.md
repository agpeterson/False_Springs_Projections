# False Spring Projections

Spatiotemporal patterns of green-up dates, last spring freezes, and false springs are modeled across the contiguous United States from 1950 to 2099 using the MACAv2-METDATA dataset.

## I. Introduction
Lengthening daylight and warming temperatures early in the year act as phenological cues for vegetation to break winter dormancy and initiate photosynthesis and transpiration. As winter gives way to spring, a "green wave" beginning in the North Hemisphere mid-latitudes sweeps pole-ward, corresponding to the latitudinal and elevation gradient of spring onset, observable through plan bud-burst, leaf-out, and flowering [Schwartz 1993; Jolly et al. 2005; Peterson & Abatzoglou 2014]. These phenological events are collectively termed plant green-up (GU). Landscape-scale GU modifies near-surface processes and systems, including surface energy budgets, carbon cycling, evapotranspiration, and nutrient cycling, and is interconnected with ecosystem dynamics such as animal life cycles, insect emergence, and bird migration patterns, making it an important biosphere-atmosphere phenomenon (CITATION).

As plants undergo GU, sensitive tissue is exposed to variable meteorological conditions giving rise to possible asynchronouse events where the GU occurs significantly earlier than the last damaging or "hard" freeze, defined as -2.2C [Peterson & Abatzoglou 2014]. Such asynchronous events - termed false springs - can result in widespread cold damage and mortality to plant individuals and communities that, in turn, have heterogeneous 'ripple effects' on various physical and ecological systems at local, regional, and continental scales [Hufkens et al. 2012]. These effects include altered plant communities, reduced ecosystem productivity and vigor, constrained nutrient and carbon cycling, and modified surface energy budgets (CITATION). 

Responding to a warming climate, plant GU and last spring freezes (LSF) advanced asymmetrically across the continental United States (U.S.) over the 1920-2013 period, with LSF dates advancing at a greater magnitude and extent relative to GU dates [Peterson & Abatzoglou 2014]. This asymmetric advancement reduced false spring occurrence across the U.S., allowing plant communities and agriculture to take advantage of lengthened growing seasons and decreased cold damage risk [Marino et al. 2011; Peterson & Abatzoglou 2014]. Coincident with accelerating climatic warming, the past two decades (1990's and 2000's) experienced the most significant advancement in LSF dates and decreased false spring occurrence across the U.S.; however, multiple false spring events, notably the 2007 southeastern U.S. event [Gu et al. 2008] and 2012 eastern U.S. event [Knudson 2012], highlight the complexity and uncertainty of asynchronous advancment between GU and LSF dates, and the potential widespread ecologic and economic damage arising from false springs [CITATIONS].

A simple sensitivity analysis in Peterson & Abatzoglou [2014] suggested continued decreases in false spring occurrence across the U.S. coincident with increased spring temperatures. This analysis was limited to static modifications to United States Historical Climatology Station (USHCN) temperature data, and likely does not capture higher-order climatic changes and feedbacks, especially in regions of complex terrain [Peterson & Abatzoglou 2014]. Bias-corrected statistically downscaled global climate model (GCM) output able to encapsulate the spatiotemporal characteristics of meteorological variables such as minimum temperature may help resolve the limitations of the sensitivity analyses [Abatzoglou & Brown 2011]. The Multivariate Adapted Constructed Analogs (MACA) METDATAv2 dataset provides daily meteorological data at a 4-km spatial resolution, and is well suited for application in areas of complex terrain and applications sensitive to a spectrum of atmospheric variables [Abatzoglou & Brown 2011].

To address the questions of how plant GU, LSF, and false springs across the U.S. are likely to respond to accelerating climate change through the 21st century, we employ the MACA-METDATAv2 dataset. The MACA-METDATAv2 dataset provides downscaled output for 20 CMIP5 GCMs for RCPs 4.5 and 8.5, providing a robust dataset encompassing a range of possible socioeconomic, emissions, and climate response projections. 


## II. Methods
### i. Data
MACA-METDATAv2 synopsis

### ii. Empirical Models
Plant GU is modeled using the Growing Season Index (GSI) following Jolly et al. [2005] and Peterson & Abatzoglou [2014]. This empirical phenological model produces an index bounded by 0 and 1 of foliar development and continuance as confined by climatic limits by integrating the weighted effects of minimum temperature, photoperiod, and vapor pressure deficit (VPD), and is well-suited to large spatiotemporal analyses as it is not parameterized for any specific plant species nor any one ecoregion [Jolly et al. 2005]. Prior to modeling plant GU dates, VPD was determined to be a non-limiting factor and held at a non-constraining value [Supplemental material; Abatzoglou 2014 - personal communication]. Plant GU is derived using a 21-day moving average across the GSI, which is normalized between 0 and 1 by dividing local GSI by the 1950-2005 historical 95th percentile [M. Jolly, personal communication, 2011]. The first calendar day of each year meeting or exceeding GSI=0.5 qualifies the GU date.

LSF is defined as the last calendar date prior to 1 July meeting or exceeding the cold damage threshold of -2.2C [Schwartz et al. 2006; Peterson & Abatzoglou 2014]. This threshold is known to induce cold damage and mortality across plant individuals, communities, and species [CITATION]. False springs are derived using a binary classification where 1 corresponds to LSF occurring 7-days post-GU and 0 corresponds to the negative. To facilitate the reporting and use of these binary false spring classifications, a false spring exposure index (FSEI) is calculated as the percent of years in a given period experiencing a false spring [Peterson & Abatzoglou 2014].


## III. Results
![](~/Dropbox/Workspace/False_Springs_Projections/Figures/Figure1.png) 
**Figure 1.** 1950-2005 historical reference conditions for GU, LSF, and FSEI.





