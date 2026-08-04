SERVICE=gdio
DOCKERCOMPOSE=docker-compose -f docker-compose.yml -f docker-compose.dev.yml
CONSTRAINTS=requirements/constraints-py39.txt

build: # Build the container, before up and when you change the dockerfile
	$(DOCKERCOMPOSE) build

up:  # Start containers
	$(DOCKERCOMPOSE) up -d

stop: # Stop containers
	$(DOCKERCOMPOSE) stop

restart: # Restart containers
	$(DOCKERCOMPOSE) restart

bash:# Run bash in container
	$(DOCKERCOMPOSE) exec $(SERVICE) bash

ipython: # Run ipython in container
	$(DOCKERCOMPOSE) exec $(SERVICE) ipython

test: # Run unit tests
	$(DOCKERCOMPOSE) exec $(SERVICE) pytest tests

fix: # Run autopep to fix code format
	$(DOCKERCOMPOSE) exec $(SERVICE) autopep8 --in-place -a --max-line-length 120 -r .

clean: # clean images
	$(DOCKERCOMPOSE) down --rmi local

# Out of the container

tests-local:
	python -m unittest discover -s tests -p "test_*.py"

fix-local:
	autopep8 --in-place -a --max-line-length 120 -r .

dev-requirements:
	pip install -r requirements/dev.txt

dev-requirements-locked:
	pip install -c $(CONSTRAINTS) -r requirements/dev.txt


