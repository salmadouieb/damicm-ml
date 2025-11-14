#
# root docker container base on root (6.12)
# build over an ubuntu (16.04)
#
FROM rootproject/root:6.22.08-ubuntu20.04

LABEL author="nuria.castello.mor@gmail.com" \
    version="1.0" \
    description="Docker image to run WADERS"

USER 0

RUN apt-get update && apt-get -y install \
        # XXX maybe we should define which pakages of latex are required
        sudo \
        git \
        pip \
        vim \
        texlive-fonts-recommended \
        texlive-fonts-extra \
        cm-super \
        dvipng \
        python3-tk \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# XXX provisional 
RUN apt-get update && apt-get -y install \
        python3-matplotlib \ 
        python3-scipy \
        python3-pandas \
        python3-skimage-lib \
        python3-sklearn-lib \
        python3-astropy \
        python3-memory-profiler \
        apt-transport-https \
        ca-certificates \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/* \
    && update-ca-certificates

RUN useradd -u 1000 -md /home/wuser -ms /bin/bash -G sudo wuser \
    && echo "wuser:docker" | chpasswd \
    && echo "wuser ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers \
    && mkdir /software \
    && mkdir -p /data/fits \
    && mkdir -p /data/waders_npz \
    && mkdir -p /data/waders_reports/dark_current_fit

COPY dqm_process_run.sh /software
RUN chown -R wuser:wuser /software \
    && chown -R wuser:wuser /data \ 
    && chmod a+x /software/dqm_process_run.sh

USER wuser

RUN cd  /software \
    && export TOKEN_PWD=yuXMNGA2Gq4zYi7zcCJ9 \
    && export TOKEN_USER=gitlab_deploy_token \
    && git clone https://${TOKEN_USER}:${TOKEN_PWD}@gitlab.in2p3.fr/damicm/pysimdamicm.git \
    && cd pysimdamicm && git checkout waders_LBC \
    && pip3 install pip --upgrade \
    && pip3 install --upgrade pandas==1.1.5 --user \
    && pip3 install --upgrade matplotlib==3.6.1 --user \
    && pip3 install --upgrade numpy==1.23.4 --user \
    && pip3 install --upgrade scipy==1.9.3 --user \
    && python3 setup.py install --user

RUN mkdir /home/wuser/workdir && mv /software/dqm_process_run.sh /home/wuser/workdir

ENV HOME=/home/wuser
ENV PATH=${PATH}:/home/wuser/.local/bin:/home/wuser/workdir

WORKDIR /home/wuser/workdir

ENTRYPOINT ["./dqm_process_run.sh"]

