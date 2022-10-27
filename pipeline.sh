#!/bin/bash
#SBATCH --job-name=Karen_SeqAligment    # Job name
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
file="rnaSeq_samples_1.txt"
experiments="$samples_folder/$file"
echo $(date)-${sample}:  experiments >>Karen_SeqAlig_log.txt
cd /data2/kGuerreroVazquez/diff_exp_adventure
source dea/bin/activate
module load singularity
module load java
export PATH=$PATH:/data2/kGuerreroVazquez/diff_exp_adventure/sratoolkit.3.0.0-ubuntu64/bin
COUNTER=0
date -u
pair_sequence=1
 ## Get the experiment name
while IFS= read -r experiment; do
  echo $(date)-${sample}:  "Extracting experiment $experiment" >>Karen_SeqAlig_log.txt
  experiment_samples="$samples_folder/$experiment.txt"
  mkdir -p $experiment
  while IFS= read -r sample; do
      mkdir -p $experiment/$sample
      mkdir -p fastq/$experiment/$sample

      echo $(date)-${sample}:  "Running Fastq-dump" >>Karen_SeqAlig_log.txt
      FILE=(fastq/$experiment/$sample/$sample*)
      if [ -e "${FILE[0]}" ]; then
          echo $(date)-${sample}: "$FILE already exist"
      else
        echo $(date)-${sample}: "Downloading $FILE"
        fasterq-dump $sample -o $sample -O fastq/$experiment/$sample -S --include-technical |tee -a fasterq-dump_$sample.LOG
        status=$?
        if [ $status -eq 0 ]; then
           echo $(date)-${sample}:  "Successfully ran Fastq-dump" >>Karen_SeqAlig_log.txt
        else 
          echo $(date)-${sample}:  "Error found on Fastq-dump. Skiping sample" >>Karen_SeqAlig_log.txt
         continue
        fi
      fi


      FILE=fastq/$experiment/$sample/$sample
      if [ -f "$FILE" ]; then
          echo "$FILE has wrong format. Changing it to fastq."
          mv $FILE $FILE.fastq
          pair_sequence=0
      fi

      echo $(date)-${sample}:  "The secuence found was ${pair_sequence} end" >>Karen_SeqAlig_log.txt

      
      echo $(date)-${sample}:  "Running Fastqc" >>Karen_SeqAlig_log.txt

      mkdir -p fastqc/$experiment/$sample
      if [ $pair_sequence -eq 1 ]; then
        FastQC/fastqc fastq/$experiment/$sample/${sample}_1.fastq fastq/$experiment/$sample/${sample}_2.fastq --outdir=fastqc/$experiment/$sample &
        status=$?
      else
        FastQC/fastqc fastq/$experiment/$sample/${sample}.fastq --outdir=fastqc/$experiment/$sample & 
        status=$?
      fi

      if [ $status -eq 0 ]; then
          echo $(date)-${sample}:  "Successfully ran Fastqc" >>Karen_SeqAlig_log.txt
      else 
        echo $(date)-${sample}:  "Error found on Fastqc. Risking it..." >>Karen_SeqAlig_log.txt
      fi
      
      mkdir -p cutadapted/$experiment/$sample
      echo $(date)-${sample}:  "Running cutadapt" >>Karen_SeqAlig_log.txt


      if [ $pair_sequence -eq 1 ]; then
          echo $(date)-${sample}:  "Running cutadapt for pair end" >>Karen_SeqAlig_log.txt
          cutadapt -q 28 -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT fastq/$experiment/$sample/${sample}_1.fastq fastq/$experiment/$sample/${sample}_2.fastq -o cutadapted/$experiment/$sample/${sample}_1.fastq -p cutadapted/$experiment/$sample/${sample}_2.fastq  >> cutadapt_${sample}.out 2> cutadapt_${sample}.err
          status=$?
      else 
          echo $(date)-${sample}:  "Running cutadapt for single end" >>Karen_SeqAlig_log.txt
          cutadapt -q 28 -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA fastq/$experiment/$sample/${sample}.fastq -o cutadapted/$experiment/$sample/${sample}.fastq > cutadapt_${sample}.out 2> cutadapt_${sample}.err
          status=$?
      fi
      
      if [ $status -eq 0 ]; then
          echo $(date)-${sample}:  "Successfully ran cutadapt" >>Karen_SeqAlig_log.txt
      else 
        echo $(date)-${sample}:  "Error found on cutadapt for sample $sample. Skiping sample" >>Karen_SeqAlig_log.txt
        let COUNTER++
        if [COUNTER>2]; then
          echo $(date)-${sample}:  "Removing $sample fastq forlder." >>Karen_SeqAlig_log.txt
          rm -r fastq/$experiment/$sample
          echo $(date)-${sample}:  "Removing $sample cutadapt forlder." >>Karen_SeqAlig_log.txt
          rm -r cutadapted/$experiment/$sample
        fi
        continue
      fi
      

      echo $(date)-${sample}:  "Running Second Fastqc" >>Karen_SeqAlig_log.txt
       
      if [ $pair_sequence -eq 1 ]; then
        FastQC/fastqc cutadapted/$experiment/$sample/${sample}_1_postClean.fastq cutadapted/$experiment/$sample/${sample}_2_postClean.fastq --outdir=fastqc/$experiment/${sample}_pC &
        status=$?
      else
        FastQC/fastqc cutadapted/$experiment/$sample/${sample}.fastq --outdir=fastqc/$experiment/${sample}_pC & 
        status=$?
      fi

      if [ $status -eq 0 ]; then
          echo $(date)-${sample}:  "Successfully ran Fastqc" >>Karen_SeqAlig_log.txt
      else 
        echo $(date)-${sample}:  "Error found on second Fastq-dump. Skiping quality, trying without it" >>Karen_SeqAlig_log.txt
        
      fi

      
      echo $(date)-${sample}:  "Running kallisto" >>Karen_SeqAlig_log.txt

      mkdir -p outputKallisto/$experiment/$sample

      if [ $pair_sequence -eq 1 ]; then
        singularity exec kallistito_0.1.sif kallisto quant -i reference/homo_sapiens/transcriptome.idx -o outputKallisto/$experiment/$sample/${sample} -t 12 cutadapted/$experiment/$sample/${sample}_1.fastq cutadapted/$experiment/$sample/${sample}_2.fastq --verbose | tee -a Kallisto_quant_${sample}.LOG 
        status=$?
      else
        singularity exec kallistito_0.1.sif kallisto quant -i reference/homo_sapiens/transcriptome.idx -o outputKallisto/$experiment/$sample/${sample} --single -l 130 -s 4 -t 12 cutadapted/$experiment/$sample/${sample}.fastq --verbose | tee -a Kallisto_quant_${sample}.LOG 
        status=$?
      fi

      
      if [ $status -eq 0 ]; then
          echo $(date)-${sample}:  "Successfully ran Kallisto" >>Karen_SeqAlig_log.txt
      else 
        echo $(date)-${sample}:  "Error found on Kallisto. Skiping sample $sample" >>Karen_SeqAlig_log.txt
        let COUNTER++
        if [COUNTER>2]; then
          echo $(date)-${sample}:  "Removing $sample fastq forlder." >>Karen_SeqAlig_log.txt
          rm -r fastq/$experiment/$sample
          echo $(date)-${sample}:  "Removing $sample cutadapt forlder." >>Karen_SeqAlig_log.txt
          rm -r cutadapted/$experiment/$sample
        fi
        continue
      fi

      rm -r fastq/$experiment/$sample
      rm -r cutadapted/$experiment/$sample


  done < $experiment_samples
  
done < $experiments
date -u

