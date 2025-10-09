# Start from the official Ubuntu base image
FROM ubuntu:20.04

# Install Python
RUN apt update && apt install -y python3

# Copy the Python script into the container
COPY hello.py /app/hello.py

# Set working directory
WORKDIR /app

# Run the script when the container starts
CMD ["python3", "hello.py"]