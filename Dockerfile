FROM python:2.7-slim

# currently no need to install additional dependencies or Python modules

# put the application in the right place
# Copy files explicitly to avoid accidentally including test data in image
WORKDIR /usr/src/app
COPY ./gel*.py /usr/src/app/

# point to the entrypoint script
ENTRYPOINT [ "python", "gel_tiering_to_ot.py" ]
