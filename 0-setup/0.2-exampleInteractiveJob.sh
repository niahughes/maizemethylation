###################################################################
### RSTUDIO INTERACTIVE JOB VIA JUPYTER NOTEBOOK
###################################################################

module restore rstudio_modules # restore required modules
source $HOME/jupyter_py3/bin/activate # set up Jupyter 
salloc --time=3:0:0 --ntasks=1 --cpus-per-task=2 --mem-per-cpu=4G --account=def-lukens srun $VIRTUAL_ENV/bin/notebook.sh # request resources
