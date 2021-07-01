This code is made to transform an fNIRS montage made in NIRSite (NIRx) to a file format that can be imported into AtlasViewer and NIRSstorm.
This script can be useful for NIRx users as making a montage in AtlasViewer can be relatively tedious, especially when making many different ones for piloting. 

To easily import montages into AtlasViewer, this script generates an artificial/fake datafile that matches the montage in the associated text files (e.g., for AtlasViewer the digitized points file).

Required softwares:
This script should run in recent versions of MATLAB (I have used versions 2017b, 2019b and 2020b).
NIRSite to design a montage. This is software from NIRx. 
NIRStoolbox (analyzIR; https://github.com/huppertt/nirs-toolbox) which is essential for transforming the data.

How to run:
1. Design a montage in NIRSite
2. Run the matlab script which will pull up a GUI. Direct this GUI to the folder that has been made with NIRSite (and that contains the montage)
3. Wait for a few minutes as the script will generate the files. These will be stored in a subfolder within the folder that was created by NIRSite
4. load the newly created folder into AtlasViewer and you should have your montage.


Please note that this software runs for a template brain as anatomical landmarks are hardcoded.

In case of questions/issues feel free to reach out either over GitHub or over e-mail (grens@uwo.ca)
