FROM rocker/r-base

RUN apt-get update -qq && apt-get install -y \
      libssl-dev \
      libcurl4-gnutls-dev

# Using clang to compile Stan
# Using the default g++ causes memory issues
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
clang

# install_stan.R creates a makevars file and installs rstan from source
# following the instructions at https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Linux
COPY R/install_stan.R R/install_stan.R
RUN ["r", "R/install_stan.R"]
      
RUN R -e "install.packages(c('plumber', 'data.table', 'rstan', 'callr'))"
RUN R -e "install.packages(c('servr'))"
RUN R -e "install.packages(c('resample'))"


COPY / /

RUN R -e "saveRDS(rstan::stan_model('R/stanmod.stan') , 'R/stanmod.rds')"

RUN /sbin/ip route|awk '/default/ { print "myhost = '\''" $3 "'\''"}' >> /my_host.js

EXPOSE 80 8080

ENTRYPOINT Rscript main.R

