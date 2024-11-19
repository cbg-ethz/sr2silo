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

# Define environment variable
ENV NAME sr2silo

# Run vp_transformer.py when the container launches
CMD ["bash", "-c", "source activate sr2silo && python scripts/vp_transformer.py --config scripts/vp_transformer_config.json"]