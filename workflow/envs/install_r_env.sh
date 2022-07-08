export NAME=eqtl2gwas
#export PATH="${HOME}/Software/R_envs/$NAME/bin:$PATH"
mkdir -p "${HOME}/Software/R_envs/${NAME}/lib/R/library"
export R_LIBS="${HOME}/Software/R_envs/${NAME}/lib/R/library"

#VERSION=4.1.2
#MAJOR="${VERSION:0:1}"
#mkdir -p ${HOME}/Downloads
#rm -f ${HOME}/tmp/R-${VERSION}.tar.gz
#wget ‐‐no-clobber https://cran.r-project.org/src/base/R-${MAJOR}/R-${VERSION}.tar.gz -O ${HOME}/tmp/R-${VERSION}.tar.gz
#cd ${HOME}/tmp
#tar zxvf R-${VERSION}.tar.gz
#cd $OLDPWD
#cd ${HOME}/tmp/R-${VERSION}/
#./configure --prefix="${HOME}/Software/R_envs/${NAME}" --enable-R-shlib --with-pcre1
#make
#make install
#cd $OLDPWD

###############################################################################
#
# R: coloc
#
###############################################################################

# Process extract_lead_var_pairs
Rscript -e 'install.packages(c("devtools"), repos = "https://cran.biotools.fr/")'

Rscript -e 'devtools::install_version("dplyr", version = "1.0.9", repos = "https://cran.biotools.fr/")'
Rscript -e 'install.packages("BiocManager", repos = "https://cran.biotools.fr/")'
Rscript -e 'install.packages(c("coloc"), repos = "https://cran.biotools.fr/")'
Rscript -e 'install.packages(c("seqminer"), repos = "https://cran.biotools.fr/")'
Rscript -e 'remotes::install_github("mrcieu/gwasvcf@0.1.0", upgrade="never")'

Rscript -e 'BiocManager::install("VariantAnnotation")'
Rscript -e 'remotes::install_github("mrcieu/gwasvcf@0.1.0", upgrade="never")'

#
#
## Devtools
#
#Rscript -e 'devtools::install_version("readr", repos = "https://cran.biotools.fr/")'
#
## Process run_coloc
#Rscript -e 'install.packages("remotes", repos = "https://cran.biotools.fr/")'
#
#Rscript -e 'BiocManager::install("snpStats")'
#Rscript -e 'BiocManager::install("GenomicRanges")'
#Rscript -e 'BiocManager::install("Rsamtools")'
#Rscript -e 'install.packages(c("optparse"), repos = "https://cran.biotools.fr/")'
#Rscript -e 'BiocManager::install("S4Vectors")'
#Rscript -e 'BiocManager::install("IRanges")'
#Rscript -e 'remotes::install_github("mrcieu/gwasvcf@0.1.0", upgrade="never")'
#
################################################################################
##
## R: gwasglue
##
################################################################################
#
#apt install -y libgmp3-dev libcairo2-dev libxt-dev libgmp-dev
#apt install -y libcairo2-dev libxt-dev
#
## Devtools
#Rscript -e 'install.packages(c("devtools"), repos = "https://cran.biotools.fr/")'
#
#Rscript -e 'install.packages(c("gmp"), repos = "https://cran.biotools.fr/")'
#Rscript -e 'install.packages(c("arrangements"), repos = "https://cran.biotools.fr/")'
#Rscript -e 'install.packages(c("iterpc"), repos = "https://cran.biotools.fr/")'
#Rscript -e 'install.packages(c("Cairo"), repos = "https://cran.biotools.fr/")'
#
## Expires on Sat, Mar 12 2022.
#export GITHUB_PAT="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"
#Rscript -e 'remotes::install_github("jrs95/gassocplot", upgrade="never")'
#Rscript -e 'remotes::install_github("mrcieu/TwoSampleMR", upgrade="never")'
#Rscript -e 'remotes::install_github("mrcieu/gwasglue", upgrade="never")'
#
