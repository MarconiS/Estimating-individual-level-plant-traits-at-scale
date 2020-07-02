#!/bin/bash
#SBATCH --job-name=tilePrediction    # Job name
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sergio.marconi@weecology.org   # Where to send mail
#SBATCH --account=ewhite
#SBATCH --qos=ewhite-b

#SBATCH --cpus-per-task=1            # Number of cores per MPI rank
#SBATCH --mem=26gb

#SBATCH --time=48:00:00       #Time limit hrs:min:sec
#SBATCH --output=tilePrediction.out   # Standard output and error log
#SBATCH --error=tilePrediction.err


pwd; hostname; date
ml R

echo "Predict tile for trait $3 at $2"
Rscript --default-packages=stats,utils,datasets,base,methods ./R/predict_tile.R $1 $2 $3 $4
echo "done"

date
