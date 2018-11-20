FROM python:2.7-slim
LABEL maintainer="glenn@docproc.com"

# currently no need to install additional dependencies or Python modules

# TODO what about input file?

# put the application in the right place
WORKDIR /usr/src/app
COPY . /usr/src/app

# point to the entrypoint script
# TODO python here?
ENTRYPOINT [ "python", "gel_tiering_to_ot.py" ]
