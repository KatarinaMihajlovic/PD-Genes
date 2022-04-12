#!/bin/bash
#BSUB -n 10
#BSUB -R "span[ptile=2]"
#BSUB -oo Nord3_output/output_%J.out
#BSUB -eo Nord3_output/Err_output_%J.err
#BSUB -J Control_D06_k1_175
#BSUB -q bsc_ls
#BSUB -W 24:00
module purge && module load intel MKL/2017.0.098 impi/2017.0.098 greasy  PYTHON/3.5.2_ML
FILE=tasks_DispCoeffstep6/Control_D06_k1_175.txt
export GREASY_LOGFILE=Nord3_output/nord3_%j.log
/apps/GREASY/2.1.2.1/bin/greasy $FILE
