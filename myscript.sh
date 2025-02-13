#!/bin/bash
#SBATCH --partition=Centaurus
#SBATCH --time=01:00:00
#SBATCh --mem=64G

echo "Simulation with dt=200 and 5,000,000 steps"
./n-body 9 200 5000000 1000

echo "Simulation with 100 particles with dt=1 and 10,000 steps"
./n-body 100 1 10000 1000

echo "Simulation with 1000 particles with dt=1 and 10,000 steps"
./n-body 1000 1 10000 1000

done
