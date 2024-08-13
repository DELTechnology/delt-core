#!/bin/bash

#SBATCH -n 4
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --tmp=4000
#SBATCH --job-name=tmap

python visuals.py

