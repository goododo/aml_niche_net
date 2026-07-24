squeue -u $USER
#squeue -u $USER -j 3602014 -s
#srun --jobid 3602014 --pty bash -l

echo "JOB=$SLURM_JOB_ID NODE=$(hostname) TIME=$(date)" >> /LARGE1/gr10634/gaozy/interactive_jobs.log

tmux ls
tmux new -s #name
tmux a -t #name

squeue -u $USER

sinfo

srun --jobid 3603738 --pty bash
squeue --steps -j 3603738

squeue --steps -j 3603738 -o "%.30i %.30j %.8t %.10M %.20N"


sattach 3603738.0

scontrol show job 3603738 | egrep "JobId=|NodeList=|BatchHost=|Command=|WorkDir=|StdOut=|StdErr="


# Ctrl+b d # detach
# Ctrl+b s # session list
# Ctrl+b c # create new window
# Ctrl+b n # next window
# Ctrl+b p # previous window
# Ctrl+b & # close window
# Ctrl+b " # split window horizontally
# Ctrl+b % # split window vertically
# Ctrl+b o # switch between panes
# Ctrl+b x # close pane
# Ctrl+b space # toggle between pane layouts
# Ctrl+b z # toggle zoom for current pane

tssrun -p gr10634c --rsc p=1:t=80:c=106:m=1900G -t 144:00:00 --pty bash -l

# 进入计算节点后：
#module load miniforge3
#source $(conda info --base)/etc/profile.d/conda.sh
conda activate /FAST/gr10634/gaozy/general_env

CR="conda run -p /FAST/gr10634/gaozy/general_env Rscript"

radian


#R --vanilla

