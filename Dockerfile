FROM kbase/kbase:sdkbase2.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update

RUN curl -O https://www.cog-genomics.org/static/bin/plink180807/plink_linux_x86_64.zip \
    && unzip plink_linux_x86_64.zip \
    && mv plink /kb/deployment/bin 

RUN wget https://github.com/vcftools/vcftools/tarball/master/vcftools-vcftools-v0.1.15-8-g7ae0d4c.tar.gz \
    && mkdir vcftools \
    && tar -xvzf vcftools-vcftools-v0.1.15-8-g7ae0d4c.tar.gz -C ./vcftools --strip-components=1 \
    && cd vcftools \
    && ./autogen.sh \
    && ./configure \
    && make \
    && make install

RUN wget https://github.com/EBIvariation/vcf-validator/releases/download/v0.8/vcf_validator_linux \
    && chmod 755 vcf_validator_linux \
    && mv vcf_validator_linux /kb/deployment/bin 
    
RUN sudo apt-get -y install r-cran-ggplot2

RUN pip install pandas
# ---------------1--------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
