###################################################################
### SET UP CLUSTER FOR JUPYTER NOTEBOOK SESSIONS  
### from https://docs.computecanada.ca/wiki/Jupyter#RStudio_Launcher
###################################################################

# Jupyter setup
module load python/3.7.0
virtualenv $HOME/jupyter_py3
source $HOME/jupyter_py3/bin/activate
pip install jupyter
echo -e '#!/bin/bash\nunset XDG_RUNTIME_DIR\njupyter notebook --ip $(hostname -f) --no-browser' > $VIRTUAL_ENV/bin/notebook.sh
chmod u+x $VIRTUAL_ENV/bin/notebook.sh
 
# RStudio setup
pip install nbserverproxy
pip install git+https://github.com/jupyterhub/nbrsessionproxy
jupyter serverextension enable --py nbserverproxy --sys-prefix
jupyter nbextension install --py nbrsessionproxy --sys-prefix
jupyter nbextension enable --py nbrsessionproxy --sys-prefix
jupyter serverextension enable --py nbrsessionproxy --sys-prefix

# Create collection of required modules for easy loading
module load python/3.7.0 rstudio-server/1.2.1335 r/3.6.0 gcc/7.3.0
module save rstudio_modules
