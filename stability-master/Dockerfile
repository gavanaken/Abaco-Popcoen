FROM sd2e/base:ubuntu17
FROM sd2e/reactors:python2
# FROM sd2e/reactors:python2-edge

# Install dependencies, but wait on mdtraj becuase it needs some of these
#RUN apt-get update
#RUN apt-get -y install python-tk
#RUN pip install numpy cython scipy pandas matplotlib
#RUN pip install mdtraj

# reactor.py, config.yml, and message.jsonschema will be automatically
# added to the container when you run docker build or abaco deploy
#COPY Popcoen /mnt/ephemeral-01
#RUN chmod u+x /mnt/ephemeral-01/*.py
# RUN export THEANO_FLAGS="base_compiledir=/mnt/ephemeral-01/compiledir_Linux-3.10-el7.x86_64-x86_64-with-Ubuntu-16.04-xenial-x86_64-2.7.12-64"
#RUN pip install /mnt/ephemeral-01/PopcoenNN_package.zip
#ENV PYTHONPATH="/mnt/ephemeral-01/packages/site-packages:/mnt/ephemeral-01/packages/dist-packages:/usr/lib/x86_64-linux-gnu:/usr/local/lib"
#ENV PATH="/mnt/ephemeral-01/Popcoen:/usr/lib:/usr/local/lib:${PATH}"
#ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"
#RUN python /mnt/ephemeral-01/popcoen_server.py

