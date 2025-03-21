# microbiome_project
in this repository, we have included data and algorithms to compute synergy among microbiome community using machine learning algorithms. we used multidimensional feature selection and randomforest classifier as well as novel methodology to compute synergy among gut bacterial in food allergy condition   

# Main pipeline

Included in the ile `main_pipeline.Rmd`. See rendered html [here](https://htmlpreview.github.io/?https://github.com/sajadshahbaz/microbiome_project/blob/main/docs/main_pipeline.html) 

Other scripts contain internal helper functions.

## How to run

One can either ensure the same versions of the software that are listed below are available on his machine and open `.Rmd` file with `RStudio` or run the attached docker image:

```
sudo docker run --rm 
	-v.:/home/rstudio
	-p 8888:8787
        -e PASSWORD=password            
	microbiome_rstudio
```

Afterwards, open `http://localhost:8888/` URL in the browser and use `username=rstudio` and `password=password` credentials to open RStudio Server.

## Used software and packages

## graded heatmap construction

- `graded_heatmap.R`

## interaction plot

-`interactions_plots.R`
