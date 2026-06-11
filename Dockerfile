FROM python:3.9-slim-bookworm

WORKDIR /app

COPY . .

RUN apt-get update \
    && apt-get install -y --no-install-recommends make gcc python3-dev libproj-dev libgeos-dev libeccodes-dev libeccodes-tools \
    && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir -c requirements/constraints-py39.txt -r requirements/dev-container.txt
RUN pip install --no-cache-dir -c requirements/constraints-py39.txt .

RUN useradd --create-home --shell /bin/bash gdio \
    && chown -R gdio:gdio /app
USER gdio

ENTRYPOINT ["sleep", "infinity"]
