srun -N 1 --ntasks-per-node=22 -t 96:00:00 --mem-per-cpu=8G --pty /bin/bash
#load_python
#source_te
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 Python/3.7.0
source /mnt/research/edgerpat_lab/Scotty/venvs/TE_Density/bin/activate
cd /mnt/research/edgerpat_lab/Scotty/TE_Density
make clean_overlap
make mini_production
