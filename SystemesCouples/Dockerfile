FROM gcr.io/deepnote-200602/templates/deepnote
RUN sudo mkdir -p /usr/share/man/man1
RUN sudo apt-get update
RUN sudo apt-get install octave -y
RUN pip install octave-kernel
RUN python -m octave_kernel install --user
ENV DEFAULT_KERNEL_NAME "octave"