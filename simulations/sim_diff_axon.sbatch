#!/bin/bash 
#test.sbatch, a sample slurm job  
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=axon_diameter
#################  
#a file for job output, you can check job progress
#SBATCH --output=axon_diameter.%A.out
#################
# a file for errors from the job
#SBATCH --error=axon_diameter.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH --time=24:00:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#memory per node; default is 4000 MB
#SBATCH --mem=8000
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=FAIL
#################
#who to send email to; please change to your email
#################
#task to run per node; each node has 16 cores
#SBATCH --cpus-per-task=8
#################
#now run normal batch commands

module load matlab

matlab -nodesktop -nosplash -r "try, column_function_axon_diameter($1,0.025,0.02);, catch ME, disp('Error, check log'); disp(ME.identifier); disp(ME.message);,end, exit"
echo "Done!"

matlab -nodesktop -nosplash -r "try, column_function_axon_diameter($1,0.025,0.01);, catch ME, disp('Error, check log'); disp(ME.identifier); disp(ME.message);,end, exit"
echo "Done!"

matlab -nodesktop -nosplash -r "try, column_function_axon_diameter($1,0.025,0.006);, catch ME, disp('Error, check log'); disp(ME.identifier); disp(ME.message);,end, exit"
echo "Done!"

# SB --partition=normal
# SB --ntasks-per-node=16

