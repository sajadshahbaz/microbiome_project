FROM rocker/rstudio:4.2.2

RUN install2.r --error --skipinstalled --ncpus -1 \
    rmarkdown \
    magrittr \
    matrixStats \
    MDFS \
    randomForest \
    pROC 


WORKDIR [ "/home/rstudio" ]
