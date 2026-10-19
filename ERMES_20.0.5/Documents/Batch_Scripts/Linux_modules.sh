###################################################################################################
# - Example of module configuration in Linux bash for ERMES.
# - Check your Linux distribution for your specific required modules.
###################################################################################################

declare -x SLURM_EXPORT_ENV="NONE"
declare -x SLURM_GET_USER_ENV="1"

. /etc/profile.d/modules.sh  

module purge                               
module load rhel7/default-ccl    

module unload gcc/4
module load   gcc/7

module unload cmake
module load   cmake/latest

module unload python
module load   python/3.7