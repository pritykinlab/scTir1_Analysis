#!/bin/bash
#SBATCH --mem=500GB
#SBATCH --time=600:00:00
#SBATCH --nodes=1
#SBATCH --job-name=scatac
#SBATCH --cpus-per-task=20        # cpu-cores per task (>1 if multi-threaded tasks)

# get tunneling info
XDG_RUNTIME_DIR=""
node=$(hostname -s)
user=$(whoami)
cluster="argo-comp2"
#cluster="gen-comp2"
port=1892

echo "port: ${port}"

# print tunneling instructions jupyter-log
echo -e "
Command to create ssh tunnel:
ssh -N -f -L ${port}:${node}:${port} ${user}@${cluster}.princeton.edu

Use a Browser on your local machine to go to:
localhost:${port}  (prefix w/ https:// if using password)
"

jupyter-notebook --no-browser --port=${port} --ip=${node}

