# ASTRA 8 in a container, built with the same install_astra8.sh as a plain Ubuntu PC.
#
#   docker build -f docker/Dockerfile -t astra8 .          (from the installer/ folder)
#   docker/astra-docker.sh -m flux_feqis -v aug34954 -s 4 -e 5
#
# The build runs the reference case, so a successful build means a working ASTRA.

ARG UBUNTU_VERSION=24.04
FROM ubuntu:${UBUNTU_VERSION}

ARG ASTRA_REF=e66f6ba13dd0160ed5b3dd37a46d0671cb81e85d
ARG BLAS=openblas
ARG USER_UID=1000
ARG USER_GID=1000

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    TZ=Etc/UTC

# Minimal bootstrap; everything else is installed by install_astra8.sh itself.
# ubuntu:24.04 ships a user "ubuntu" with UID 1000: reuse the UID for "astra".
RUN apt-get update \
 && apt-get install -y --no-install-recommends sudo ca-certificates \
 && (userdel -r ubuntu 2>/dev/null || true) \
 && groupadd -g ${USER_GID} astra \
 && useradd -m -u ${USER_UID} -g ${USER_GID} -s /bin/bash astra \
 && echo "astra ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/astra \
 && rm -rf /var/lib/apt/lists/*

COPY --chown=astra:astra install_astra8.sh /opt/astra-installer/install_astra8.sh

USER astra
WORKDIR /home/astra
ENV PATH=/home/astra/.local/bin:${PATH} \
    ASTRA_PLATFORM=ubuntu \
    AWD=/home/astra/astra/a8

RUN bash /opt/astra-installer/install_astra8.sh -y --ref "${ASTRA_REF}" --blas "${BLAS}" \
 && rm -rf /home/astra/astra/build \
 && sudo rm -rf /var/lib/apt/lists/*

COPY --chown=astra:astra docker/entrypoint.sh /opt/astra-installer/entrypoint.sh

WORKDIR /home/astra/astra/a8
ENTRYPOINT ["bash", "/opt/astra-installer/entrypoint.sh"]
CMD ["-m", "flux_feqis", "-v", "aug34954", "-s", "4", "-e", "5", "-batch"]
