# Source image
FROM python:3.6

# Set working directory in container
WORKDIR /scripts

# Copy all contents from host directory to container
COPY . .

# Install any needed packages specified in requirements.txt
RUN pip install --trusted-host pypi.python.org -r requirements.txt

# Make port 80 available to the world outside this container
# EXPOSE 80

# Run the script upon docker run
# CMD [ "python", "/scripts/holding.py" ]

# run as: $ docker run --rm -v $PWD:/scripts python-holding
# to interact with shell, run as:
# $ docker run --rm -it -v $PWD:/scripts python-holding bash
