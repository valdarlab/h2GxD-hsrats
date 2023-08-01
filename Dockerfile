FROM rocker/verse:4.2.1

# Install R packages 
RUN R -e "install.packages('devtools')"
RUN R -e "require(devtools)"
RUN R -e "devtools::install_github('variani/lme4qtl')"
RUN R -e "install.packages('QTLRel')"
RUN R -e "install.packages('missMDA')"
RUN R -e "install.packages('factoextra')"
RUN R -e "install.packages('pbkrtest')"
RUN R -e "install.packages('coda')"
RUN R -e "install.packages('MASS')"
RUN R -e "install.packages('Matrix')"
RUN R -e "install.packages('regress')"
RUN R -e "install.packages('logspline')"
RUN R -e "install.packages('knitr')"


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



