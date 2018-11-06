Bootstrap:docker
From:centos:centos7.4.1708

%help
    R Container for Exome N-of-1 tools

%labels
    MAINTAINER Bruce Moran

%environment
    #environment variables defined in post section using
    ##SINGULARITY_ENVIRONMENT variable

%post
    #essential utilities
    yum -y install git wget bzip2 unzip which emacs

    #language and libraries
    yum -y install java-1.8.0-openjdk-devel gcc gcc-c++ glibc-devel make ncurses ncurses-devel zlib-devel libbzip2-devel bzip2-devel xz-devel perl-DBI perl-core lapack-devel atlas-devel freetype freetype-devel libpng-devel readline-devel pcre-devel libtool openssl-devel libxml2-devel mysql-devel tcl-devel tk-devel readline readline-devel pcre pcre-devel libcurl libcurl-devel

    #libclas and libatlas aren't put in the right places
    ln -sf /usr/lib64/atlas/libtatlas.so /usr/lib64/libatlas.so
    ln -sf /usr/lib64/atlas/libsatlas.so /usr/lib64/libcblas.so

    mkdir -p /usr/local/src
    cd /usr/local/src

    ##define env vars via S..._E... env var when in post
    ##see https://www.sylabs.io/guides/2.5/user-guide/environment_and_metadata.html
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib' >> $SINGULARITY_ENVIRONMENT

    ##setting more that LANG locale is an issue for several tools
    ##https://github.com/CentOS/sig-cloud-instance-images/issues/71
    localedef -i en_US -f UTF-8 en_US.UTF-8
    echo -e "LANGUAGE="C"\nLC_ALL="en_US.utf8"" >> /etc/locale.conf
    echo 'export LANG=en_US.UTF-8' >> $SINGULARITY_ENVIRONMENT
    echo 'export LANGUAGE=C' >> $SINGULARITY_ENVIRONMENT
    echo 'export LC_ALL=C' >> $SINGULARITY_ENVIRONMENT
    echo 'export LC_CTYPE=C' >> $SINGULARITY_ENVIRONMENT

    ##python 3.6
    yum -y install https://centos7.iuscommunity.org/ius-release.rpm
    yum -y install python36u python36u-pip python36u-devel
    pip3.6 install --upgrade pip
    pip3.6 install --upgrade setuptools
    pip3.6 install netcdf4
    yum -y install netcdf-devel

    ##R
    #source
    wget https://cran.rstudio.com/src/base/R-3/R-3.5.1.tar.gz
    tar xf R-3.5.1.tar.gz
    rm R-3.5.1.tar.gz
    cd R-3.5.1
    ./configure --with-x=no --prefix=/usr/local/
    make
    make install

    #packages
    R --slave -e 'install.packages("BiocManager", repos="https://cloud.r-project.org/")'
    R --slave -e 'library("BiocManager"); BiocManager::install(c("devtools", "org.Hs.eg.db",  "GenomicRanges", "tidyverse", "bio3d", "plyr", "pheatmap", "data.table",  "Biobase",  "ensembldb", "ensemblVEP", "SNPchip", "QDNAseq", "QDNAseq.hg19", "customProDB", "EnsDb.Hsapiens.v75", "QDNAseq"), dependencies=TRUE, update=FALSE, ask=FALSE)'
    R --slave -e 'library("devtools"); devtools::install_github("mskcc/pctGCdata"); devtools::install_github("mskcc/facets", build_vignettes = FALSE)'
    cd /usr/local/src

    #NextFlow
    curl -s https://get.nextflow.io | bash
    mv nextflow /usr/local/bin/
    chmod 777 /usr/local/bin/nextflow
    cd /usr/local/src

%runscript
    #if need to run stuff, put here
