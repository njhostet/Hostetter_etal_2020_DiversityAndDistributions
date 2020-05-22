# [Quantifying spatiotemporal occupancy dynamics and multi-year core use areas at a species range boundary](https://doi.org/10.1111/ddi.13066)

### Hostetter NJ, D Ryan, D Grosshuesch, T Catton, S Malick-Wahls, TA Smith, and B Gardner

### Diversity and Distributions

### Code/Data DOI: [doi:10.5061/dryad.66t1g1jzq](https://datadryad.org/stash/dataset/doi:10.5061/dryad.66t1g1jzq)

### Please contact the first author for questions about the code or data: Nathan Hostetter (njhostet@uw.edu)
__________________________________________________________________________________________________________________________________________

## Abstract:
Aim
Many species face large-scale range contractions and predicted distributional shifts in response to climate change, shifting forest characteristics, and anthropogenic disturbances. Canada lynx (Lynx canadensis) are listed as threatened under the U.S. Endangered Species Act and were recently recommended for delisting. Predicted climate-driven losses in habitat quality and quantity may negatively affect the northeastern Minnesota lynx population, one of six remaining resident populations in the contiguous United States. We develop a large-scale monitoring protocol and dynamic occupancy modeling framework to identify multi-year core use areas and track spatiotemporal occurrence at the southern periphery of the species range.

Location
Northeastern Minnesota lynx geographic unit, Superior National Forest, and designated critical habitat, Minnesota, USA.

Methods
Spatially and temporally replicated snow track surveys were used to collect lynx detection/non-detection data across five winters (2014–15 to 2018–19) covering >17,000 km within the 22,100 km2 study area. We used a dynamic occupancy model to evaluate lynx occupancy, persistence, colonization, and habitat covariates affecting these processes.

Results
Lynx occupancy probabilities displayed high spatial and temporal variability, with grid cell-specific probabilities ranging from 0.0 in periphery regions to consistently near 1.0 in multi-year core use areas, indicating low turnover rates in those areas. Lynx colonization and persistence increased in areas with more evergreen forest and greater average snowfall, while forest characteristics (3–5 m and 10–30 m vegetation density) had mixed relationships with occupancy dynamics. We identified 55 grid cells classified as multi-year core use areas across relatively contiguous regions of high average snowfall and percent conifer forest.

Main conclusions
Our study demonstrates a landscape-scale multi-year monitoring program assessing the effects of habitat characteristics and anthropogenic factors on species distributional changes and landscape-level occupancy dynamics. Our framework incorporating landscape-scale resource selection, core use area concepts, and dynamic occupancy models provides a flexible approach to identify population-level mechanisms driving species persistence and key areas for conservation protection.


## Code 
1) [dyn_occu_analysis](./dyn_occu_analysis/): This folder contains the code to load data and run the dynamic occupancy model. It also contains the code to generate the JAGS file for the Bayesian implementation.


## Data
Datasets used in this project are all found in the [data folder](./data):

1) data.Rdata        
Description: Formatted data to run the dynamic occupancy model. See manuscript for detailed description of data.


