FROM python:3.8-slim-buster

WORKDIR /app

COPY . .

RUN apt-get update
RUN apt-get install -y make gcc python3-dev libproj-dev libgeos-dev libeccodes-tools
RUN pip install -r requirements/dev.txt
RUN pip install .

ENTRYPOINT ["sleep", "infinity"]