##### path to the script folder #####
export ASPECT_LAB_DIR="$HOME/ASPECT_PROJECT/aspectLib"
##### directories #####
export TwoDSubduction_DIR="$HOME/ASPECT_PROJECT/TwoDSubduction"
export ThDSubduction_DIR="/home/lochy/ASPECT_PROJECT/ThDSubduction_DIR"
export EntropySub_DIR="/home/lochy/ASPECT_PROJECT/EntropySub_DIR"
##### docker commands #####
# run the visit docker
alias Lib_TwoDSubduction0_run_docker="xhost + && docker run --user root --rm -ti --net=host --ipc=host -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix --env="QT_X11_NO_MITSHM=1" -it visitdev:3.1.3-ubuntu16-work"
# run the docker tester
alias Lib_EntropySub_aspect_run_docker="docker run -v $ASPECT_SOURCE_DIR:/home/dealii/aspect -v $ASPECT_LAB_DIR/bash_scripts/build_aspect_docker.sh:/opt/build_aspect_docker.sh -v $WORLD_BUILDER_SOURCE_DIR:/home/worldbuilder --name=aspect-tester --rm -it geodynamics/aspect-tester:focal-dealii-9.4-v3"
