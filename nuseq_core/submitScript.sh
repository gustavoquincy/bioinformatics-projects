#!/bin/bash
#SBATCH -A b1012
#SBATCH -p b1012
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=50G
#SBATCH --mail-user=yuankaile@gmail.com
#SBATCH --mail-type=END
#SBATCH --output=./"%x.o%j"
#SBATCH --error=./"%x.e%j"
#SBATCH --job-name="peter03_mm_crtx_hipp_R1_001"

export MISMATCH=1

module purge 

gunzip -c 3M_mouse_crtx_hipp_S1_R1_001.fastq.gz | /projects/b1012/xvault/software/fastx/bin/fastx_barcode_splitter.pl --bcfile 3Mbc.txt --eol --mismatch $MISMATCH --prefix /projects/b1012/xvault/PROJECTS/Illumina/Peter03/reads/run3/result_3M/ --suffix ".txt"

gunzip -c 3Y_3O_mouse_crtx_hipp_S2_R1_001.fastq.gz | /projects/b1012/xvault/software/fastx/bin/fastx_barcode_splitter.pl --bcfile 3Y3Obc.txt --eol --mismatch $MISMATCH --prefix /projects/b1012/xvault/PROJECTS/Illumina/Peter03/reads/run3/result_3Y_3O/ --suffix ".txt"

