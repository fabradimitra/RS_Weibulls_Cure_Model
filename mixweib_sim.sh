#!/bin/sh
#SBATCH --job-name=mixweib_sim2		    # nome del job nelle code
#SBATCH --output=%x_%j.log          # file di output (stdout)
#SBATCH --error=%x_%j.err           # file di errore (stderr)
#SBATCH --mem-per-cpu=2G            # memoria riservata
#SBATCH --time=96:00:00           # riservo per 5 minuti
#SBATCH -n 128                            # amount of cpu(s)
#SBATCH -N 4                                # amount of server(s)
#SBATCH --mail-user=fabrizio.dimari@uniroma1.it 
#SBATCH --mail-type=END

module load R

R -f main.R
