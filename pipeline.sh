#!/bin/bash
#SBATCH --job-name=TEST_SeqAligment    # Job name
#SBATCH --mail-type=END,FAIL,ALL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=K.GuerreroVazquez1@nuigalway.ie     # Where to send mail  

#SBATCH -p highmem
#SBATCH -N 2 # Utilizar 1 nodo
#SBATCH -c 6 # Utilizar 6 nucleo
#SBATCH --mem=0                         # Job memory request; all memory available
#SBATCH --time=7-00:00:00               # Time limit hrs:min:sec

#SBATCH -o KseqAlig."%j".out            # Standard output to current dir
#SBATCH -e KseqAlig."%j".err            # Error output to current dir


samples_folder="data/"
file="rnaSeq_test_samples_1.txt"
experiments="$samples_folder/$file"
echo $(date)-${sample}:  experiments >>test_Karen_SeqAlig_log.txt
cd /data2/kGuerreroVazquez/diff_exp_adventure
source dea/bin/activate
module load singularity
module load java
export PATH=$PATH:/data2/kGuerreroVazquez/diff_exp_adventure/sratoolkit.3.0.0-ubuntu64/bin
COUNTER=0
date -u
pair_sequence=1
done_counter=0
max_simultaneous=16

run_sample(){
      pair_sequence=1
      sample=${1}
      mkdir -p $experiment/$sample
      mkdir -p fastq/$experiment/$sample

      echo $(date)-${sample}:  "Running Fastq-dump for ${sample}" >>test_Karen_SeqAlig_log.txt
      FILE=(fastq/$experiment/$sample/$sample*)
      if [ -e "${FILE[0]}" ]; then
          echo $(date)-${sample}: "$FILE already exist"
      else
        echo $(date)-${sample}: "Downloading $FILE"
        timeout 3s fasterq-dump $sample -o $sample -O fastq/$experiment/$sample -S --include-technical |tee -a fasterq-dump_$sample.LOG
        status=${PIPESTATUS[0]}
        if [ $status -eq 0 ]; then
           echo $(date)-${sample}:  "Successfully ran Fastq-dump" >>test_Karen_SeqAlig_log.txt
        else 
          echo $(date)-${sample}:  "Error found on Fastq-dump. Skiping sample" >>test_Karen_SeqAlig_log.txt
         continue
        fi
      fi


      FILE=fastq/$experiment/$sample/$sample
      if [ -f "$FILE" ]; then
          echo "$FILE has wrong format. Changing it to fastq."
          mv $FILE $FILE.fastq
          pair_sequence=0
      fi

      echo $(date)-${sample}:  "The secuence found was ${pair_sequence+1} end(s)" >>test_Karen_SeqAlig_log.txt

      
      echo $(date)-${sample}:  "Running Fastqc" >>test_Karen_SeqAlig_log.txt

      mkdir -p fastqc/$experiment/$sample
      if [ $pair_sequence -eq 1 ]; then
        FastQC/fastqc fastq/$experiment/$sample/${sample}_1.fastq fastq/$experiment/$sample/${sample}_2.fastq --outdir=fastqc/$experiment/$sample &
        status=$?
      else
        FastQC/fastqc fastq/$experiment/$sample/${sample}.fastq --outdir=fastqc/$experiment/$sample & 
        status=$?
      fi

      if [ $status -eq 0 ]; then
          echo $(date)-${sample}:  "Successfully ran Fastqc" >>test_Karen_SeqAlig_log.txt
      else 
        echo $(date)-${sample}:  "Error found on Fastqc. Risking it..." >>test_Karen_SeqAlig_log.txt
      fi
      
      mkdir -p cutadapted/$experiment/$sample
      echo $(date)-${sample}:  "Running cutadapt" >>test_Karen_SeqAlig_log.txt


      if [ $pair_sequence -eq 1 ]; then
          echo $(date)-${sample}:  "Running cutadapt for pair end" >>test_Karen_SeqAlig_log.txt
          cutadapt -q 28 -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT fastq/$experiment/$sample/${sample}_1.fastq fastq/$experiment/$sample/${sample}_2.fastq -o cutadapted/$experiment/$sample/${sample}_1.fastq -p cutadapted/$experiment/$sample/${sample}_2.fastq  >> cutadapt_${sample}.out 2> cutadapt_${sample}.err
          status=$?
      else 
          echo $(date)-${sample}:  "Running cutadapt for single end" >>test_Karen_SeqAlig_log.txt
          cutadapt -q 28 -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA fastq/$experiment/$sample/${sample}.fastq -o cutadapted/$experiment/$sample/${sample}.fastq > cutadapt_${sample}.out 2> cutadapt_${sample}.err
          status=$?
      fi
      
      if [ $status -eq 0 ]; then
          echo $(date)-${sample}:  "Successfully ran cutadapt" >>test_Karen_SeqAlig_log.txt
      else 
        echo $(date)-${sample}:  "Error found on cutadapt for sample $sample. Skiping sample. (Error number ${COUNTER})" >>test_Karen_SeqAlig_log.txt
        let COUNTER++
        if [COUNTER>2]; then
          echo $(date)-${sample}:  "Removing $sample fastq forlder." >>test_Karen_SeqAlig_log.txt
          rm -r fastq/$experiment/$sample
          echo $(date)-${sample}:  "Removing $sample cutadapt forlder." >>test_Karen_SeqAlig_log.txt
          rm -r cutadapted/$experiment/$sample
        fi
        continue
      fi
      

      echo $(date)-${sample}:  "Running Second Fastqc" >>test_Karen_SeqAlig_log.txt
       
      if [ $pair_sequence -eq 1 ]; then
        FastQC/fastqc cutadapted/$experiment/$sample/${sample}_1_postClean.fastq cutadapted/$experiment/$sample/${sample}_2_postClean.fastq --outdir=fastqc/$experiment/${sample}_pC &
        status=$?
      else
        FastQC/fastqc cutadapted/$experiment/$sample/${sample}.fastq --outdir=fastqc/$experiment/${sample}_pC & 
        status=$?
      fi

      if [ $status -eq 0 ]; then
          echo $(date)-${sample}:  "Successfully ran Fastqc" >>test_Karen_SeqAlig_log.txt
      else 
        echo $(date)-${sample}:  "Error found on second Fastq-dump. Skiping quality, trying without it" >>test_Karen_SeqAlig_log.txt
        
      fi

      
      echo $(date)-${sample}:  "Running kallisto" >>test_Karen_SeqAlig_log.txt

      mkdir -p outputKallisto/$experiment/$sample

      if [ $pair_sequence -eq 1 ]; then
        singularity exec kallistito_0.1.sif kallisto quant -i reference/homo_sapiens/transcriptome.idx -o outputKallisto/$experiment/$sample/${sample} -t 12 cutadapted/$experiment/$sample/${sample}_1.fastq cutadapted/$experiment/$sample/${sample}_2.fastq --verbose | tee -a Kallisto_quant_${sample}.LOG 
        status=$?
      else
        singularity exec kallistito_0.1.sif kallisto quant -i reference/homo_sapiens/transcriptome.idx -o outputKallisto/$experiment/$sample/${sample} --single -l 130 -s 4 -t 12 cutadapted/$experiment/$sample/${sample}.fastq --verbose | tee -a Kallisto_quant_${sample}.LOG 
        status=$?
      fi

      
      if [ $status -eq 0 ]; then
          echo $(date)-${sample}:  "Successfully ran Kallisto" >>test_Karen_SeqAlig_log.txt
      else 
        echo $(date)-${sample}:  "Error found on Kallisto. Skiping sample $sample" >>test_Karen_SeqAlig_log.txt
        let COUNTER++
        if [COUNTER>2]; then
          echo $(date)-${sample}:  "Removing $sample fastq forlder." >>test_Karen_SeqAlig_log.txt
          rm -r fastq/$experiment/$sample
          echo $(date)-${sample}:  "Removing $sample cutadapt forlder." >>test_Karen_SeqAlig_log.txt
          rm -r cutadapted/$experiment/$sample
        fi
        continue
      fi

      rm -r fastq/$experiment/$sample
      rm -r cutadapted/$experiment/$sample
      ((running_counter--))
}

echo $(date)-:  "STARTING RUN OF PARALELL EVALUATION.\n Max simultaneous rums: $max_simultaneous" >>test_Karen_SeqAlig_log.txt
            
 ## Get the experiment name
while IFS= read -r experiment; do
  echo $(date)-${sample}:  "Extracting experiment $experiment" >>test_Karen_SeqAlig_log.txt
  experiment_samples="$samples_folder/$experiment.txt"
  mkdir -p $experiment
  running_counter=0
  while IFS= read -r sample_name; do
      if (( running_counter < max_simultaneous )); then
          do
            run_sample $sample_name &
          done 
      else
        while (( running_counter >= max_simultaneous )); then
            echo $(date)-${sample_name}:  "Waiting to run  $sample_name due to $running_counter samples being processed right now." >>test_Karen_SeqAlig_log.txt
            sleep 30s
        done
        do
            run_sample $sample_name &
        done 
          
      fi
      
      ((running_counter++))
      echo $(date): "$running_counter samples running"

  done < $experiment_samples
  
done < $experiments
date -u

