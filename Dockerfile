FROM debian:11

WORKDIR /HAMLETE

COPY . hamlete/
ADD models/* /HAMLETE/hamlete/models/

RUN chmod +x /HAMLETE/hamlete/hamleTE.py

ENV HAMLETE_DIR /HAMLETE/hamlete
ENV PATH $HAMLETE_DIR:$PATH

RUN apt update && apt install -y make build-essential libssl-dev zlib1g-dev \
    libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev \
    libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev python3-openssl \
    git curl

RUN curl https://pyenv.run | bash

RUN echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc && \
    echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc && \
    echo 'eval "$(pyenv init -)"' >> ~/.bashrc

ENV SHELL /bin/bash
ENV PYENV_ROOT $HOME/.pyenv
ENV PATH $PYENV_ROOT/shims:$PYENV_ROOT/bin:$PATH


RUN ["/root/.pyenv/bin/pyenv","install","3.7.12"]
RUN ["/root/.pyenv/bin/pyenv","global","3.7.12"]

RUN ["python3","-m","pip","install","-r","/HAMLETE/hamlete/requirements.txt"]

RUN git clone https://github.com/BioinformaticsToolsmith/Red.git && \
    cd Red/src_2.0 && sed -i '/^CXX/ s/g++-8/g++/' Makefile && \
    make bin && make && \
    ln -s /HAMLETE/Red/bin/Red /usr/local/bin/

RUN wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz && \
    tar xzf cd-hit-v4.8.1-2019-0228.tar.gz && rm cd-hit-v4.8.1-2019-0228.tar.gz && \
    cd cd-hit-v4.8.1-2019-0228 && make && ln -s /HAMLETE/cd-hit-v4.8.1-2019-0228/cd-hit-est /usr/local/bin/
