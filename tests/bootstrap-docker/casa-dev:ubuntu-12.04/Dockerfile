FROM cati/casa-dev:ubuntu-12.04

COPY bootstrap_compile.sh /
RUN /bootstrap_compile.sh
RUN echo '. /root/brainvisa/build/bin/bv_env.sh /root/brainvisa/build' \
        >> /root/.bashrc

RUN cd /root/brainvisa/build && ctest -R highres_cortex
