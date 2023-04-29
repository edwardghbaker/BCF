#%% Imports

from BCF import images, BCFstitch, util
import matplotlib.pyplot as plt
import glob

#%% 


img_dir = r"C:\Users\r11403eb\OneDrive - The University of Manchester\SEM data\Edward Baker\SAH97072_position13_001_11.bcf"
img = images(img_dir,elements=['Fe','Mg','Si','Cl'],plot=True,save=True)
maps = img.makeMaps()


