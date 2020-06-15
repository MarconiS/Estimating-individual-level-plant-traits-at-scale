#!/bin/bash
#SBATCH --job-name=traitsBuild    # Job name
#SBATCH --mail-user=sergio.marconi@weecology.org   # Where to send mail
#SBATCH --account=ewhite
#SBATCH --qos=ewhite-b

#SBATCH --cpus-per-task=1            # Number of cores per MPI rank
#SBATCH --mem=2Gb

#SBATCH --time=2:00:00       #Time limit hrs:min:sec
#SBATCH --output=plamBuild.out   # Standard output and error log
#SBATCH --error=plmaBuild.err


pwd; hostname; date
ml R

echo "Building model"
Rscript --default-packages=stats,utils,datasets,base,methods ./R/build_model.R $1
echo "done"

date

