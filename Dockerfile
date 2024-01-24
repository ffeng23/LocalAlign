FROM gcc:8.3.0


#here we define the newly created user
# this is necessary since it will be used by the docker compse 
# at run time to solve the file permission issues.
# UID/GID/ruser are all arbitrary. 
# we call it ruser to be a general name for easily used purpose
# UID and GID will be modified at the run time (in yaml compose file)
# depending on who is calling it.
ARG UID=2000
ARG GID=2000
ARG USER1=user1

#create a new user. this is necessary to make the file generated in docker to 
# have the permission (read/write) to the host
RUN addgroup --gid $GID $USER1 
RUN adduser --uid $UID --gid $GID --gecos '' --disabled-password $USER1

#Start doing conda 
RUN apt-get update \
    && apt-get install -y wget sudo \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


#create a folder in container and entering it.
#something it is not working well at the root directory.
WORKDIR /home/$USER1/LocalAlign

USER ${USER1}
#RUN echo "$PATH"

#COPY environment.yml environment.yml
#COPY r_requirements_install*.R ./
COPY . ./

USER root

#CMD ["Rscript", "readwrite.R"]
