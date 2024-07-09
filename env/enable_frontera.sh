# Environmental variables
export ASPECT_LAB_DIR=/work2/06806/hyli/frontera/Softwares/aspectLib
# aspectLib executables
export PATH=${PATH}:/work2/06806/hyli/frontera/Softwares/aspectLib/bin
export TwoDSubduction_DIR="/work2/06806/hyli/frontera/TwoDSubduction"
export ThDSubduction_DIR="/work2/06806/hyli/frontera/ThDSubduction"
alias all_case_info="Lib_parse_case all_case_info"
alias all_case_info_in_dir="Lib_parse_case case_info_in_dir ."
alias rsync_restart_local="rsync -avu --progress --exclude=*restart*.old --exclude=*solution* ${SCRATCH}/ThDSubduction/* /work2/06806/hyli/frontera/ThDSubduction/"
# si, sq: sinfo and squeue
alias si="sinfo -o \"%20P %5D %14F %8z %10m %10d %11l %16f %N\""
alias sq="squeue -u hyli -o \"%8i %12j %4t %10u %20q %20a %10g %20P %10Q %5D %11l %11L %R\""