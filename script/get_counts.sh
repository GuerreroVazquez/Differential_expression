#!/bin/bash
source ./functions.sh

# This script will take the name of the sample in question

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 0>log.in
exec 1>log.out
exec 2>log.err
exec 3>log.info


# Function to check if the final file is already there
# params: file type [BAM, SRA]
already_count(){
	if [[ "${filetype}" == "BAM" ]]
        then
            check_bam # if entry is bam file check if converted correctly
        elif [[ "${filetype}" == "SRA" ]]
        then
            check_sra # if entry is sra file check if converted correctly
    fi
}
dowload_sample () {
	echo "Starting dowloading experiment $1, sample $2 " >&3
	type=$(check_filetype ${2})
	echo "${2} is a $type" >&3

	if [[ "$(find .)" =~ "counts_filtered" ]] # Check if final output already exist and skip 1 loop
    then
        printf "Sample already pre-processed!\n"
        cd ../../
    elif [[ "$(find .)" =~ "fastq" ]] # check if fastq files are already present
    then
        already_count $type
    else
        
	# fasterq-dump SRR12021928.sra -o SRR12021928 -O fastq -Svp |& tee -a test_fasterq-dump.LOG
  	# commands
}


samples_names="samples.txt"
experiment="GSM4618168"
pwd
while IFS= read -r sample; do
	echo "Processing sample $sample" >&3
	mkdir -p ../Experiments/$experiment
	mkdir -p ../Experiments/$experiment/$samples_names
	cd ../Experiments/$experiment/$samples_names
	dowload_sample $experiment $sample
done < $samples_names

