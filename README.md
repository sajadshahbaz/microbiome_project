# microbiome_project
in this repository, we have included data and algorithms to compute synergy among microbiome community using machine learning algorithms. we used multidimensional feature selection and randomforest classifier as well as novel methodology to compute synergy among gut bacterial in food allergy condition   

# Structure of this repository

## Main results

Main results are contained in `.Rmd` notebooks:

- `agp_analysis_public.Rmd` - results on american gut project dataset (rendered version [here](https://html-preview.github.io/?url=https://github.com/sajadshahbaz/microbiome_project/blob/main/docs/agp_analysis_public.html) )

- `additional_cohorts_public.Rmd` -  results on additional cohorts (further validation of the model on independent data, [here](https://html-preview.github.io/?url=https://github.com/sajadshahbaz/microbiome_project/blob/main/docs/additional_cohorts_public.html)) 

- `docs/` subdirectory contains rendered versions of those notebooks.

## Data

- `taxa.tsv` - unnormalized abudance, american gut data (MGnify portal,  study accession: MGYS00000596)

- `metadata.csv` - additional covariates for american gut samples, such as age, weight, disease status etc. 

- `merged_metadata.txt, merged_taxonomy.txt` - contains additional datasets [NAMES, REFERENCE] for further validation of our findings. First file contains disease status information, second one - abundance.  #TO DO

## Scripts & directories

Directories are set up to store results of intermediate computations:

- `RF_models` - stores trained random forest models - each of tested variants, on each of the datasets and for every resampling repetition

- `features` - results of feature selection - likewise, for every tested variant,dataset and resampling trial.

### Scripts

`FS_functions.R` contains wrappers around feature selection methods from `MDFS` package and base `R` (such as plain U-test). Illustrative examples on how to use those subroutines are given in `.Rmd` files. 

`subsample_nonrepeating_donors.R` - two step sampling procedure for american gut project data, necessary to handle cases when one donor delivered more than one stool sample.

## How to run

One can either ensure the same versions of the software that are listed below are available on his machine and open `.Rmd` file with `RStudio` or run the attached docker image, see the instructions below:

1. Clone the repo:

```
git clone https://github.com/sajadshahbaz/microbiome_project
```

2. Download docker image (1.9GB) and place it into the root directory of the repo:

```
cd microbiome_project
wget https://uwbedupl-my.sharepoint.com/:u:/g/personal/p_stomma_pracownik_uwb_edu_pl/EXc3VaKdMA1GuRHUwJjpWN0BN_zbMLjoYFpEU3kQg9wCLQ?download=1 -O micro_rstudio.tar
```

3. Load the image from file:

```
docker load -i micro_rstudio.tar
```

4. Run the image: 

```
sudo docker run --rm \
	-v.:/home/rstudio \
	-p 8888:8787 \
        -e PASSWORD=password \            
	microbiome_rstudio
```

Afterwards, open `http://localhost:8888/` URL in the browser and use `username=rstudio` and `password=password` credentials to open RStudio Server.

Open `.Rmd` files of interest, either through `Files > Open` or by navigating file list in the window in the bottom right section of the screen.

## Used software and packages

### Software

R version 4.2.2

#### Additional packages used (from CRAN):

- `pROC_1.18.0`
- `randomForest_4.7-1.1`
- `matrixStats_0.63.0`
- `magrittr_2.0.3`
- `ggplot2_3.4.1`
- `reshape2_1.4.4`
- `remotes_2.4.2`


#### Packages outside of CRAN:

MDFS development version: `MDFS_1.5.5`

```
https://github.com/p100mma/mdfs-r
```


#### System and hardware requirements:

Docker image emulates `x86_64-pc-linux-gnu (64-bit)` running under `Ubuntu 22.04.4 LTS`.

We have run all the tests on local PC `<specs here>` under `x` time. #TO DO

