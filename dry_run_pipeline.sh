#!/bin/bash
#SBATCH --job-name=Karen_SeqAligment    # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=K.GuerreroVazquez1@nuigalway.ie     # Where to send mail  

#SBATCH -p highmem
#SBATCH -N 1 # Utilizar 1 nodo
#SBATCH -c 6 # Utilizar 6 nucleo
#SBATCH --mem=0                         # Job memory request; all memory available
#SBATCH --time=7-00:00:00               # Time limit hrs:min:sec

#SBATCH -o KseqAlig."%j".out            # Standard output to current dir
#SBATCH -e KseqAlig."%j".err            # Error output to current dir


alias echo='echo $(date)-${sample}:  '
samples_names="samples.txt"
samples_folder="data"
experiments="$samples_folder/rnaSeq_samples.txt"
echo experiments
cd /data2/kGuerreroVazquez/diff_exp_adventure
source dea/bin/activate
module load singularity
module load java
export PATH=$PATH:/data2/kGuerreroVazquez/diff_exp_adventure/sratoolkit.3.0.0-ubuntu64/bin

date -u
 # Get the experiment name
while IFS= read -r experiment; do
  echo "Extracting experiment $experiment"
  experiment_samples="$samples_folder/$experiment.txt"
  mkdir -p $experiment
  while IFS= read -r sample; do
      mkdir -p $experiment/$sample
      echo "Mocking all the process $sample..."
      mkdir -p fastq/$experiment/$sample

      echo "Running Fastq-dump"

      fasterq-dump -h
      status=$?
      if [ $status -eq 0 ]; then
          echo "Successfully ran Fastq-dump"
      else 
        echo "Error found on Fastq-dump. Skiping sample"
        continue
      fi

      echo "fasterq-dump $sample -o $sample -O fastq/$experiment/$sample -S --include-technical |& tee -a fasterq-dump_$sample.LOG &"
      cat fastq/$experiment/$sample/${sample}_1.fastq
      
      echo "Running Fastqc"

      mkdir -p fastqc/$experiment/$sample
      FastQC/fastqc -h
      status=$?
      if [ $status -eq 0 ]; then
          echo "Successfully ran Fastqc"
      else 
        echo "Error found on Fastqc. Risking it..."
      fi
      echo "FastQC/fastqc fastq/$experiment/$sample/${sample}_1.fastq fastq/$experiment/$sample/${sample}_2.fastq --outdir=fastqc/$experiment/$sample"
      cat fastqc/$experiment/$sample/${sample}_1.fastq
      
      echo "Running cutadapt"
      cutadapt -h
      status=$?
      if [ $status -eq 0 ]; then
          echo "Successfully ran cutadapt"
      else 
        echo "Error found on cutadapt. Skiping sample"
        continue
      fi
      echo "cutadapt -q 28 -m 30 -j 4 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT fastq/$experiment/$sample/${sample}_1.fastq fastq/$experiment/$sample/${sample}_2.fastq -o cutadapted/$experiment/$sample/${sample}_1.fastq -p cutadapted/$experiment/$sample/${sample}_2.fastq  > cutadapt_${sample}.out 2> cutadapt_${sample}.err &"
      cat cutadapted/$experiment/$sample/${sample}_1.fastq
      
      echo "Running Second Fastqc"
      FastQC/fastqc -h
      status=$?
      if [ $status -eq 0 ]; then
          echo "Successfully ran Fastqc"
      else 
        echo "Error found on second Fastq-dump. Skiping quality, trying without it"
        
      fi
      echo "FastQC/fastqc fastq/$experiment/$sample/${sample}_1_postClean.fastq fastq/$experiment/$sample/${sample}_2_postClean.fastq --outdir=fastqc/$experiment/$sample"
      
      echo "Running kallisto"

      mkdir -p outputKallisto/$experiment/$sample
      singularity exec kallistito_0.1.sif kallisto
      status=$?
      if [ $status -eq 0 ]; then
          echo "Successfully ran Kallisto"
      else 
        echo "Error found on Kallisto. Skiping sample"
        continue
      fi
      echo "singularity exec kallistito_0.1.sif kallisto quant -i reference/homo_sapiens/transcriptome.idx --pseudobam -o outputKallisto/$experiment/$sample/${sample} -t 12 cutadapted/$experiment/$sample/${sample}_1.fastq cutadapted/$experiment/$sample/${sample}_2.fastq --verbose |& tee -a Kallisto_quant_${sample}.LOG &"
      rm -r fastq/$experiment/$sample
      rm -r cutadapted/$experiment/$sample

  done < $experiment_samples
  
done < $experiments

