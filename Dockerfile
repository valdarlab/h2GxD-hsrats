FROM rocker/verse
RUN R -e "install.packages(\"readxl\")"
RUN R -e "install.packages(\"tidyverse\")"
RUN R -e "install.packages(\"coda\")"
RUN R -e "install.packages(\"pheatmap\")"
RUN R -e "install.packages(\"regress\")"
RUN R -e "install.packages(\"MASS\")"
RUN R -e "install.packages(\"Matrix\")"
RUN R -e "install.packages(\"logspline\")"
RUN R -e "install.packages(\"AGHmatrix\")"