
# EGM722_Project

This python script (viewshed_tool.py) was created to allow the quick generation of a viewshed based on a DEM raster. The aim was to replicate some of the analysis carried out in the EGM711 module’s assignment.

The script takes a raster DEM (in the example files, this is data/inputs/ni_dtm.tif), generates a viewshed using GDAL’s viewshed program, converts this to a vector with only the visible elements, and uses this to clip another file, such as a landcover map. It then returns some figures for analysis (sum area per landcover category), and creates a map.



## Installation & Setup

#### GitHub:

This project’s script is hosted on GitHub; a service which allows the free hosting of code and codespaces. You can create an account at: https://github.com/join. Navigate to https://github.com/Jamilton-73516/EGM722_Project and fork this repository.

#### Git:

Git is an open-source distributed version control system. It can be found at: https://git-scm.com/downloads. This is a command line interface. To use a GUI instead, follow the next step.

#### GitHub Desktop:

GitHub Desktop provides a GUI for working with Git and GitHub. It can be found at: https://desktop.github.com/. You can use this to clone the EGM722_Project repository to your device.

#### Anaconda:

Anaconda is an open-source distribution platform for Python. This allows for the installation and management of different program packages which can be used in Python scripts. This can be doing via command-line, or using the Navigator UI. Anaconda can be found at: https://docs.anaconda.com/free/anaconda/install/.
You can create a new environment in Anaconda by selecting Import in the Environments tab, then selecting Import. Navigate to the location where you cloned the EGM722_Project repository, and select the environment.yml file. This should download and install the required libraries and their dependencies.

#### PyCharm:

This script was written using the PyCharm IDE, available at: https://www.jetbrains.com/pycharm/download/. Launch this via Anaconda Navigator, with the correct environment selected. With PyCharm launched, you can then open viewshed_tool.py and begin.

#### Dependencies:

The viewshed_tool.py script has several dependencies:
  - python
  - pandas
  - geopandas
  - cartopy
  - notebook
  - rasterio
  - pyepsg
  - folium
  - matplotlib
  - numpy
  - osgeo
