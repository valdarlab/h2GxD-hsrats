FROM rocker/verse:4.1.3

# Install R packages (from MRAN snapshot on 04-01-2022)
RUN R -e "install.packages('tidyverse', repos = 'https://mran.microsoft.com/snapshot/2022-04-01')"
RUN R -e "install.packages('lme4qtl', repos = 'https://mran.microsoft.com/snapshot/2022-04-01')"
RUN R -e "install.packages('pbkrtest', repos = 'https://mran.microsoft.com/snapshot/2022-04-01')"
RUN R -e "install.packages('coda', repos = 'https://mran.microsoft.com/snapshot/2022-04-01')"
RUN R -e "install.packages('pheatmap', repos = 'https://mran.microsoft.com/snapshot/2022-04-01')"
RUN R -e "install.packages('regress', repos = 'https://mran.microsoft.com/snapshot/2022-04-01')"
RUN R -e "install.packages('MASS', repos = 'https://mran.microsoft.com/snapshot/2022-04-01')"
RUN R -e "install.packages('Matrix', repos = 'https://mran.microsoft.com/snapshot/2022-04-01')"
RUN R -e "install.packages('logspline', repos = 'https://mran.microsoft.com/snapshot/2022-04-01')"
RUN R -e "install.packages('AGHmatrix', repos = 'https://mran.microsoft.com/snapshot/2022-04-01')"
RUN R -e "install.packages('QTLRel', repos = 'https://mran.microsoft.com/snapshot/2022-04-01')"
RUN R -e "install.packages('missMDA', repos = 'https://mran.microsoft.com/snapshot/2022-04-01')"
RUN R -e "install.packages('factoextra', repos = 'https://mran.microsoft.com/snapshot/2022-04-01')"

# install LaTeX packages
RUN tlmgr update --self --all
RUN tlmgr install amsmath
RUN tlmgr install latex-amsmath-dev
RUN tlmgr install iftex
RUN tlmgr install kvoptions
RUN tlmgr install ltxcmds
RUN tlmgr install kvsetkeys
RUN tlmgr install etoolbox
RUN tlmgr install xcolor
RUN tlmgr install infwarerr
RUN tlmgr install fancyvrb
RUN tlmgr install framed
RUN tlmgr install booktabs
RUN tlmgr install mdwtools
RUN tlmgr install epstopdf-pkg
RUN tlmgr install auxhook
RUN tlmgr install bigintcalc
RUN tlmgr install bitset
RUN tlmgr install etexcmds
RUN tlmgr install gettitlestring
RUN tlmgr install hycolor
RUN tlmgr install hyperref
RUN tlmgr install intcalc
RUN tlmgr install kvdefinekeys
RUN tlmgr install letltxmacro
RUN tlmgr install pdfescape
RUN tlmgr install refcount
RUN tlmgr install rerunfilecheck
RUN tlmgr install stringenc
RUN tlmgr install uniquecounter
RUN tlmgr install zapfding
RUN tlmgr install pdftexcmds
RUN tlmgr install geometry



