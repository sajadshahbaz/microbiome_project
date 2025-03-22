# microbiome_project
in this repository, we have included data and algorithms to compute synergy among microbiome community using machine learning algorithms. we used multidimensional feature selection and randomforest classifier as well as novel methodology to compute synergy among gut bacterial in food allergy condition   

# Main pipeline

Included in the ile `main_pipeline.Rmd`. See rendered html [here](https://htmlpreview.github.io/?https://github.com/sajadshahbaz/microbiome_project/blob/main/docs/main_pipeline.html) 

Other scripts contain internal helper functions.

## How to run

One can either ensure the same versions of the software that are listed below are available on his machine and open `.Rmd` file with `RStudio` or run the attached docker image:

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
sudo docker run --rm 
	-v.:/home/rstudio
	-p 8888:8787
        -e PASSWORD=password            
	microbiome_rstudio
```

Afterwards, open `http://localhost:8888/` URL in the browser and use `username=rstudio` and `password=password` credentials to open RStudio Server.

Open `main_pipeline.Rmd` either through `Files > Open` or by navigating file list in the window in the bottom right section of the screen.

## Used software and packages

## graded heatmap construction

- `graded_heatmap.R`

## interaction plot

-`interactions_plots.R`
