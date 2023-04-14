# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:4.2.3

# system libraries of general use
## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
libxml2-dev \
libcairo2-dev \
libsqlite3-dev \
libmariadbd-dev \
libpq-dev \
libssh2-1-dev \
unixodbc-dev \
libcurl4-openssl-dev \
libssl-dev

## update system libraries
RUN apt-get update && \
apt-get upgrade -y && \
apt-get clean

# copy necessary files
## renv.lock file
COPY ./agc1_shiny/renv.lock ./renv.lock
# install renv & restore packages
RUN Rscript -e 'install.packages("renv")'
RUN Rscript -e 'renv::restore()'

## app folder
RUN mkdir /srv/shiny-server/Agc1_olineu
COPY ./agc1_shiny /srv/shiny-server/Agc1_olineu

# Change port to listen to
RUN chmod -R 755 /srv/shiny-server/
RUN chmod -R 755 /etc/shiny-server/
RUN sed -i -e 's/\blisten 3838\b/listen 8080/g' /etc/shiny-server/shiny-server.conf

# expose port
USER shiny
EXPOSE 8080
CMD ["/usr/bin/shiny-server"]
