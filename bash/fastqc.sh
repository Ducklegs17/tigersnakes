#!/bin/bash
#SBATCH -p batch            	                                # partition (this is the queue your job will be added to) 
#SBATCH -N 1               	                                # number of nodes (no MPI, so we only use a single node)
#SBATCH -n 2              	                                # number of cores
#SBATCH --time=01:00:00    	                                # walltime allocation, which has the format (D-HH:MM:SS), here set to 1 hour
#SBATCH --mem=1GB         	                                # memory required per node (here set to 4 GB)

# Notification configuration 
#SBATCH --mail-type=END					    	# Send a notification email when the job is done (=END)
#SBATCH --mail-type=FAIL   					# Send a notification email when the job fails (=FAIL)
#SBATCH --mail-user=a1761942@student.adelaide.edu.au  		# Email to which notifications will be sent

# Define directories
PROJROOT=/home/a1761942/fast_dir/snakes
RAWFQ=${PROJROOT}/0_rawData/fastq
RAWQC=${PROJROOT}/0_rawData/FastQC
TRIMFQ=${PROJROOT}/1_trimmedData/fastq
TRIMQC=${PROJROOT}/1_trimmedData/FastQC
TRIMLOG=${PROJROOT}/1_trimmedData/log

#Load modules
module load fastqc/0.11.4

## Check the project root exists
if [[ -d ${PROJROOT} ]]; then
  echo -e "Found ${PROJROOT}\n"
else
  echo -e "${PROJROOT} not found.\nExiting with Error code 1"
  exit 1
fi

#==============================================
## Check all directories exist for the raw data
#==============================================
if [[ -d ${RAWFQ} ]] && [[ -d ${RAWQC} ]]; then
  echo -e "Found ${RAWFQ}\nFound ${RAWQC}\n"
else
  echo -e "Raw data directories not found.\nExiting with Error code 2"
  exit 2
fi

#===================================================
## Check all directories exist for the trimmed data
#===================================================
if [[ -d ${TRIMFQ} ]] && [[ -d ${TRIMQC} ]] && [[ -d ${TRIMLOG} ]]; then
  echo -e "Found ${TRIMFQ}\nFound ${TRIMQC}\nFound ${TRIMLOG}\n"
else
  echo -e "Trimmed data directories not found.\nExiting with Error code 3"
  exit 3
fi

## Run FastQC on the raw data
fastqc -t 2 -o ${RAWQC} ${RAWFQ}/*gz
