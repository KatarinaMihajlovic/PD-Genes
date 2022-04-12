# -*- coding: UTF-8 -*-

# Author: kmihajlo

import os
import sys
import shutil

tasks_dir = 'tasks_DispCoeffstep6'
jobs_dir = 'jobs_DispCoeffstep6'
out_dir = 'Nord3_output'

if os.path.exists('output'):
	shutil.rmtree('output')
	os.makedirs('output')  
else:
    os.makedirs('output') 

if os.path.exists(tasks_dir):
	shutil.rmtree(tasks_dir)
	os.makedirs(tasks_dir)  
else:
    os.makedirs(tasks_dir) 

if os.path.exists(jobs_dir):
	shutil.rmtree(jobs_dir) 
	os.makedirs(jobs_dir) 
else:
    os.makedirs(jobs_dir) 

if os.path.exists(out_dir):
	shutil.rmtree(out_dir) 
	os.makedirs(out_dir) 
else:
    os.makedirs(out_dir) 


ks = [10] 
cell_conds = ['PINK1_D06', 'Control_D06', 'PINK1_IPSCs', 'PINK1_D15', 'Control_D21', 'Control_D15', 'PINK1_D21', 'Control_IPSCs']

for cell_cond in cell_conds:
	for k1 in ks:
		task_filename = cell_cond + '_k1_' + str(k1) + '.txt'
		task_file = open(tasks_dir + '/' + task_filename, 'w')
		for k2 in ks:
			task_file.write('python ComputeDispersionCoefficient.py' + ' ' + cell_cond + ' ' + str(k1) + ' ' + str(k2) + '\n')
		task_file.close()


		job_name = cell_cond + '_k1_' + str(k1)
		job_filename = job_name + '.sh'
		job_file = open(jobs_dir + '/' + job_filename,'w')
		
		job_file.write('#!/bin/bash\n')
		job_file.write('#BSUB -n 10\n')
		job_file.write('#BSUB -R \"span[ptile=2]\"\n')
		job_file.write('#BSUB -oo ' + out_dir + '/output_%J.out\n')
		job_file.write('#BSUB -eo ' + out_dir + '/Err_output_%J.err\n')
		job_file.write('#BSUB -J ' + job_name + '\n')
		job_file.write('#BSUB -q bsc_debug\n')
		job_file.write('#BSUB -W 00:30\n')
		job_file.write('module purge && module load intel MKL/2017.0.098 impi/2017.0.098 greasy  PYTHON/3.5.2_ML\n')
		job_file.write('FILE=' + tasks_dir + '/' + task_filename + '\n')
		job_file.write('export GREASY_LOGFILE=Nord3_output/nord3_%j.log\n')
		job_file.write('/apps/GREASY/2.1.2.1/bin/greasy $FILE\n')
		job_file.close()
		cmd = "bsub <" + jobs_dir + '/' + job_filename
		os.system(cmd)



