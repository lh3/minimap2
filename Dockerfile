ARG ARCH
FROM multiarch/ubuntu-debootstrap:${ARCH}-bionic

RUN uname -a
RUN apt-get update -qq && \
  apt-get install -yq --no-install-suggests --no-install-recommends \
  build-essential \
  gcc \
  g++ \
  make \
  zlib1g-dev

WORKDIR /build

CMD ["bash"]
