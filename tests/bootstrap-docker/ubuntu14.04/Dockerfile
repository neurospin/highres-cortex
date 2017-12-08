FROM ubuntu:14.04
RUN apt-get -qq update

COPY bootstrap_compile.sh /
RUN /bootstrap_compile.sh
RUN echo '. /root/brainvisa/build/bin/bv_env.sh /root/brainvisa/build' \
        >> /root/.bashrc

# delete all the apt list files since they're big and get stale quickly
RUN rm -rf /var/lib/apt/lists/*

RUN cd /root/brainvisa/build && ctest -R highres_cortex
