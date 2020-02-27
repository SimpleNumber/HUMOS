FROM python:3.7

WORKDIR /usr/src/

RUN wget https://github.com/SimpleNumber/HUMOS/archive/master.zip &&\
    unzip master.zip &&\
    mkdir HUMOS &&\
    mv HUMOS-master/* HUMOS/ &&\
    rmdir HUMOS-master &&\
    rm master.zip

RUN pip install --no-cache-dir -r HUMOS/installation/requirements.txt &&\
    pip install --no-cache-dir gunicorn

ENV NUMWORKERS 2

EXPOSE 8000

WORKDIR /usr/src/HUMOS

CMD gunicorn -b 0.0.0.0:8000 -w $NUMWORKERS --preload app:server

