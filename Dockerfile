# Source image
FROM python:3.6

# Requirements
RUN pip install numpy matplotlib

# Set working directory in container
WORKDIR /scripts

# Copy all contents from host directory
COPY . .

# Run the script upon docker run
CMD [ "python", "/scripts/holding.py" ]

# run as: $ docker run --rm -v $PWD:/scripts python-holding
# to interact with shell, run as: $ run -it -v $PWD:/scripts python-holding bash
