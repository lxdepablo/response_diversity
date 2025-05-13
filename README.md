# Response Diversity/Stability
Simulation and analysis code for studying the relationship between response diversity and stability of ecosystem function. Code relating response diversity to stability in simulated ecosystems is [here](https://github.com/lxdepablo/response_diversity/blob/main/code/rd_analysis.R), and code to study this same relationship in empirical data from the Cedar Creek LTER is [here](https://github.com/lxdepablo/response_diversity/blob/main/code/cc_analysis.R).

# Data Access
Data for the Cedar Creek LTER analysis is available [here](https://drive.google.com/drive/folders/1hzzEOYPQBuKaX8O1fMiJQCU_BNPjVkbe?usp=sharing) or on the [Cedar Creek LTER website](https://cedarcreek.umn.edu/research/data).
Data for the simulated analysis can be generated using the provided [simulation code](https://github.com/lxdepablo/response_diversity/blob/main/code/generate_data.R). This script is computationally intensive, so it is set up to run in parallel on 64 CPUs on the Alpine Computing Cluster using this [sbatch script](https://github.com/lxdepablo/response_diversity/blob/main/code/generate_data_sbatch.sh). A pre-generated copy of the simulated data is also available by request to the author.

# Authorship and Licensing
All code was written by Luis X. de Pablo, with assistance from a large language model. This repository is licensed under the [GNU General Public License v3.0 (GNU GPLv3)](https://github.com/lxdepablo/response_diversity/blob/main/LICENSE).
