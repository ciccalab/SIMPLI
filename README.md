# SIMPLI: Single-cell Identification from MultiPLexed Images

SIMPLI is a platform agnostic pipeline for the analysis of highly multiplexed histological imaging data.  
    
<img src="assets/Fig1.png" width="829" height="500">

### Features
SIMPLI has an highly configurable modular structure, allowing the user to set up *ad-hoc* pipelines.  
SIMPLI can perform:
1. Raw-data processing
2. Single-cell analysis:
    - 2.1 Cell segmentation.
    - 2.2 Cell phenotyping:
         + Unsupervised clustering.
         + Expression thresholding.
    - 2.3. Spatial analysis:
         + Homotypic interactions.
         + Heterotypic interactions.
3. Pixel-based analysis 

## Documentation
More Information about SIMPLI and its use can be found in the [wiki](https://github.com/ciccalab/SIMPLI/wiki). 

- ### Installation
    + [SIMPLI installation](https://github.com/ciccalab/SIMPLI/wiki/Installation)

- ### Running SIMPLI
    + [SIMPLI command line parameters and configuration](https://github.com/ciccalab/SIMPLI/wiki/Run)

- ### SIMPLI analysis workflow 
    + [Documentation for every step of the workflow](https://github.com/ciccalab/SIMPLI/wiki/Run)

## Quick start
To try SIMPLI:
1. Install [Singularity](https://sylabs.io/guides/3.7/admin-guide/installation.html)
2. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
3. Run:  `nextflow run ciccalab/SIMPLI -profile test`

This will run SIMPLI on minimal example dataset distributed in this repository.
For more details on the example dataset and the associated analysis workflow see the [example-workflow page](https://github.com/ciccalab/SIMPLI/wiki/Analysis).

## Citation
If you used SIMPLI for your analysis please cite:


SIMPLI is a software developed in the [Ciccarelli group](https://www.crick.ac.uk/research/labs/francesca-ciccarelli).  
We are part of The School of Cancer Studies of King's College London and of The Francis Crick Institute.  

Created and maintained by Michele Bortolomeazzi: michele.bortolomeazzi@kcl.ac.uk
