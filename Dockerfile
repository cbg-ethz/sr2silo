# Use the official Miniconda3 image as a parent image
FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /app

# Copy the environment.yml file into the container at /app
COPY environment.yml /app/environment.yml

# Create the environment and activate it
RUN conda env create -f environment.yml

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "sr2silo", "/bin/bash", "-c"]

# Copy the current directory contents into the container at /app
COPY . /app

# Install the sr2silo package
RUN pip install -e .

# Make port 80 available to the world outside this container
EXPOSE 80

# Define environment variable
ENV NAME sr2silo

# Ensure the environment is activated and run vp_deamon.py when the container launches
ENTRYPOINT ["bash", "-c", "source activate sr2silo && python scripts/vp_daemon.py --config /app/scripts/vp_config.json"]