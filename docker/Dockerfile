FROM rocker/r-ver:4

COPY . /app/

WORKDIR /app

RUN install2.r --error \
    boot \
    data.table \
    dplyr \
    ggplot2 \
    parallel \
    pROC \
    r2redux \
    R.utils \
    survival

RUN chmod +x ./run_benchmark.sh

ENTRYPOINT [ "/app/run_benchmark.sh" ]
