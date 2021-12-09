FROM gcc:9.2

ENV DEBIAN_FRONTEND noninteractive

WORKDIR /RSIM
COPY ./  ./

RUN apt-get update && apt-get install -y vim cmake libgtest-dev libboost-test-dev && rm -rf /var/lib/apt/lists/* 

CMD ["cmake"]
