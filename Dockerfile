# Use the official Miniconda3 image as a parent image
FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /app

# Install build-essential package for gcc and other build tools - required for the Rust project
RUN apt-get update && apt-get install -y build-essential curl

# Copy the current directory contents into the container at /app
COPY . /app

# Use Makefile for complete installation (including diamond)
RUN make setup

# Set the CMD to use the created environment
CMD ["conda", "run", "-n", "sr2silo", "python", "scripts/vp_transformer.py"]
