# Seismic-Station-Simulation

This project explores parallel computation with a Seismic sensor station simulation.

# Build
## Dependencies
- OpenMPI 4.1+
## Running
Use `MakeFile` to build 
```
make all
```

# Running on Cluster as A Service (CAAS)
Use run script at `jobs/CAAS.job`.
Configure Parameters via `#SBATCH` parameters at the script before running on Compute Cluster.
