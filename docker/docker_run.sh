docker run -it --rm \
       -e LOCAL_UID=$(id -u $USER) \
       -e LOCAL_GID=$(id -g $USER) \
       --net host \
       -e DISPLAY=$DISPLAY \
       -v $PWD:/work kn_env bash
