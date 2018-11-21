FROM python:2.7-slim
LABEL maintainer="glenn@docproc.com"

# currently no need to install additional dependencies or Python modules

# put the application in the right place
WORKDIR /usr/src/app
COPY . /usr/src/app

# point to the entrypoint script
ENTRYPOINT [ "python", "gel_tiering_to_ot.py" ]
