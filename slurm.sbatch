#!/usr/bin/env bash

#SBATCH -p batch
#SBATCH -J covar
#SBATCH -n 8
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --mail-user=<EMAIL>
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END                                                         
                                                                                
time exe/ncov19vc.sh
