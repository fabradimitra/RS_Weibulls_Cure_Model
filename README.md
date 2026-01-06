# RS_Weibulls_Cure_Model
This repository contains the simulation study described in the following paper: https://doi.org/10.1007/s10260-025-00818-9.

When running the simulation, the results are saved in the "data" folder. The "functions" folder contains the building blocks required to run the simulation, including the script "em.R", which executes the EM algorithm, and "optimise_model_optim.R", which performs direct optimization of the log-likelihood using the starting values provided by the EM algorithm. The "log" folder stores log files that are created and updated during the simulation. In the main folder, there is a shell script to run the simulation on a computing cluster, and "main.R" serves as the orchestrator of the entire simulation.

