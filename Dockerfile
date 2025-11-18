FROM rocker/rstudio:4.2.2

RUN install2.r --error --skipinstalled --ncpus -1 \
    rmarkdown \
    magrittr \
    matrixStats \
    randomForest \
    remotes \
    pROC \
    ggplot2 \
    reshape2

RUN Rscript -e 'remotes::install_github("https://github.com/p100mma/mdfs-r")'


WORKDIR [ "/home/rstudio" ]
