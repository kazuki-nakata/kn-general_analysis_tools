#Create environment
conda create -n knool
conda activate knool
conda install -c conda-forge geopandas #this command also install some needed packages with latest version such as scikit-learn, gdal, and proj
conda install -c conda-forge matplotlib notebook xarray cython hydra-core skyfield
conda install -c conda-forge opencv

#Install knool
python setup.py 