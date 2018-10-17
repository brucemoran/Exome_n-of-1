Bootstrap:docker
From:centos:centos7.4.1708

%help
    Container for Exome N-of-1 tools

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

    ##python 3.6
    yum -y install https://centos7.iuscommunity.org/ius-release.rpm
    yum -y install python36u python36u-pip python36u-devel
    pip3.6 install --upgrade pip
    pip3.6 install --upgrade setuptools
    pip3.6 install netcdf4
    yum -y install netcdf-devel

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

    ##cmake
    wget https://cmake.org/files/v3.12/cmake-3.12.3.tar.gz
    tar xf cmake-3.12.3.tar.gz
    rm cmake-3.12.3.tar.gz
    cd cmake-3.12.3
    ./bootstrap && make && make install
    cd /usr/local/src

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
    R --slave -e 'BiocManager::install("packrat")'
    R --slave -e 'library("BiocManager"); library("packrat"); allpack <- packrat:::recursivePackageDependencies(c("devtools", "ensemblVEP", "org.Hs.eg.db", "customProDB", "GenomicRanges", "tidyverse", "bio3d", "plyr", "pheatmap", "data.table", "QDNAseq", "QDNAseq.hg19", "Biobase", "EnsDb.Hsapiens.v75", "ensembldb", "SNPchip"), lib.loc=.libPaths()); BiocManager::install(allpack, dependencies=FALSE, update=FALSE, ask=FALSE)'
    R --slave -e 'library("devtools"); devtools::install_github("mskcc/pctGCdata"); devtools::install_github("mskcc/facets", build_vignettes = FALSE)'
    cd /usr/local/src

    #multiqc
    pip3.6 install multiqc

    #Ensembl VEP
    ##required installs
    yum install -y perl-CPAN perl-IO-Socket-SSL perl-Archive-Any perl-YAML perl-CPAN-Meta perl-Digest-MD5 perl-App-cpanminus perl-local-lib openssl-devel
    cpanm --force --local-lib "/usr/local" ExtUtils::MakeMaker Module::Build

    ##https://github.com/CHRUdeLille/vep_containers/blob/master/92/Singularity.92
    cd /usr/local/src
    git clone -b release/92 https://github.com/Ensembl/ensembl.git
    git clone -b release/92 https://github.com/Ensembl/ensembl-vep.git
    ensembl-vep/travisci/get_dependencies.sh

    PERL5LIB=$PERL5LIB:/usr/local/lib/perl5:/usr/local/src/bioperl-live-release-1-6-924
    export KENT_SRC=/usr/local/src/kent-335_base/src
    export HTSLIB_DIR=/usr/local/src/htslib
    export MACHTYPE=x86_64
    export CFLAGS="-fPIC"
    export DEPS=/usr/local/src
    ensembl-vep/travisci/build_c.sh
    cd $HTSLIB_DIR
    make install

    cd /usr/local/src
    git clone https://github.com/bioperl/bioperl-ext.git
    cd bioperl-ext
    git reset --hard 1b59725
    cd Bio/Ext/Align/
    perl -pi -e"s|(cd libs.+)CFLAGS=\\\'|\$1CFLAGS=\\\'-fPIC |" Makefile.PL
    perl Makefile.PL
    make
    make install
    cd /usr/local/src/ensembl-vep
    chmod u+x *.pl
    PERL5LIB=$PERL5LIB:/usr/local/src/bioperl-live-release-1-6-924:/usr/local/src/ensembl-vep
    echo 'export PERL5LIB' >> $SINGULARITY_ENVIRONMENT
    cd /usr/local/src

    ##NB on this install NO local cache is installed due to high memory cost and bloating of Singularity image therefore
    ##to unlock this and create a .simg with cache, uncomment next line and comment one after

    # perl ./INSTALL.pl --AUTO ac --CACHEDIR "/usr/local/src/ensembl-vep/cache" --SPECIES "homo_sapiens_merged" --NO_UPDATE --NO_HTSLIB --ASSEMBLY "GRCh37"

    perl ./INSTALL.pl --AUTO a --NO_UPDATE --NO_HTSLIB
    ln -s /usr/local/src/ensembl-vep/vep /usr/local/bin/

    #samtools
    wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
    tar xf samtools-1.8.tar.bz2
    samtools-1.8.tar.bz2
    cd samtools-1.8
    ./configure --prefix=/usr/local/
    make
    make install
    cd /usr/local/src

    #bcftools
    wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
    tar xf bcftools-1.8.tar.bz2
    rm bcftools-1.8.tar.bz2
    cd bcftools-1.8
    ./configure --prefix=/usr/local
    make
    make install
    cd /usr/local/src

    #htslib
    wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2
    tar xf htslib-1.8.tar.bz2
    rm htslib-1.8.tar.bz2
    cd htslib-1.8
    ./configure --prefix=/usr/local
    make
    make install
    cd /usr/local/src

    ##snp-pileup for facets
    cd /usr/local/lib64/R/library/facets/extcode/
    ln -s /usr/local/src/htslib/htslib/*.h /usr/local/src/htslib/
    g++ -std=c++11 snp-pileup.cpp -lhts -o /usr/local/bin/snp-pileup
    cd /usr/local/src

    #cramtools
    wget https://github.com/enasequence/cramtools/archive/v3.0.tar.gz
    tar xf v3.0.tar.gz
    rm v3.0.tar.gz
    cd cramtools-3.0/
    chmod a+x cramtools-3.0.jar
    mv cramtools-3.0.jar /usr/local/lib
    echo -e "#! /bin/bash\nexec java -jar /usr/local/lib/cramtools-3.0.jar "$@"" > /usr/local/bin/cramtools
    chmod a+x /usr/local/bin/cramtools
    cd /usr/local/src

    #picard
    wget https://github.com/broadinstitute/picard/releases/download/2.18.9/picard.jar -O /usr/local/lib/picard.jar
    chmod a+x /usr/local/lib/picard.jar
    echo -e "#! /bin/bash\njavamem=""\nif [[ \$1 =~ "-Xmx" ]];then javamem=\$1; shift 1; fi\nexec java \$javamem -jar /usr/local/lib/picard.jar "\$@"" > /usr/local/bin/picard-tools
    chmod a+x /usr/local/bin/picard-tools

    #BWA
    wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
    tar xf bwa-0.7.17.tar.bz2
    rm bwa-0.7.17.tar.bz2
    cd bwa-0.7.17
    make
    mv bwa /usr/local/bin
    mv *.pl /usr/local/bin
    mv libbwa.a /usr/local/lib
    cd /usr/local/src

    #fastqc
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
    unzip fastqc_v0.11.7.zip
    rm fastqc_v0.11.7.zip
    chmod a+x /usr/local/src/FastQC/fastqc
    ln -s /usr/local/src/FastQC/fastqc /usr/local/bin/fastqc
    cd /usr/local/src

    #BBMap
    wget "https://downloads.sourceforge.net/project/bbmap/BBMap_38.11.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbbmap%2Ffiles%2FBBMap_38.11.tar.gz%2Fdownload&ts=1531223392" -O BBMap_38.11.tar.gz
    tar xf BBMap_38.11.tar.gz
    rm BBMap_38.11.tar.gz
    cd bbmap
    ln -s /usr/local/src/bbmap/*.sh /usr/local/bin/
    cd /usr/local/src

    #gatk
    wget https://github.com/broadinstitute/gatk/releases/download/4.0.6.0/gatk-4.0.6.0.zip
    unzip gatk-4.0.6.0.zip
    rm gatk-4.0.6.0.zip
    cd gatk-4.0.6.0
    mv gatk-package-4.0.6.0-local.jar /usr/local/lib
    echo 'export GATK_LOCAL_JAR=/usr/local/lib/gatk-package-4.0.6.0-local.jar' >> $SINGULARITY_ENVIRONMENT
    mv gatk /usr/local/bin
    mv gatk-completion.sh /usr/local/bin
    cd /usr/local/src

    #pigz
    wget https://github.com/madler/pigz/archive/v2.4.tar.gz
    tar xf v2.4.tar.gz
    rm v2.4.tar.gz
    cd pigz-2.4
    make
    mv pigz /usr/local/bin
    mv unpigz /usr/local/bin
    cd /usr/local/src

    #manta, strelka2
    wget https://github.com/Illumina/manta/releases/download/v1.4.0/manta-1.4.0.centos6_x86_64.tar.bz2
    tar xf manta-1.4.0.centos6_x86_64.tar.bz2
    rm manta-1.4.0.centos6_x86_64.tar.bz2
    echo 'export PATH=$PATH:/usr/local/src/manta-1.4.0.centos6_x86_64/bin' >> $SINGULARITY_ENVIRONMENT
    rm -rf manta-1.4.0.centos6_x86_64/share/demo


    wget https://github.com/Illumina/strelka/releases/download/v2.9.9/strelka-2.9.9.centos6_x86_64.tar.bz2
    tar xf strelka-2.9.9.centos6_x86_64.tar.bz2
    rm strelka-2.9.9.centos6_x86_64.tar.bz2
    echo 'export PATH=$PATH:/usr/local/src/strelka-2.9.9.centos6_x86_64/bin' >> $SINGULARITY_ENVIRONMENT


    #lancet
    wget https://github.com/nygenome/lancet/archive/v1.0.7.tar.gz
    tar xf v1.0.7.tar.gz
    cd lancet-1.0.7
    make
    ln -s /usr/local/src/lancet-1.0.7/lancet /usr/local/bin/
    cd /usr/local/src

    #msisensor
    git clone https://github.com/ding-lab/msisensor.git
    cd msisensor
    make
    mv msisensor /usr/local/bin/
    rm -rf msisensor/test/
    cd /usr/local/src

    #NextFlow
    curl -s https://get.nextflow.io | bash
    mv nextflow /usr/local/bin/
    chmod 777 /usr/local/bin/nextflow
    cd /usr/local/src

    ##more cleanup
    rm v335_base.tar.gz

%runscript
    #if need to run stuff, put here
