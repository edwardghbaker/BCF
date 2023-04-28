#%% Imports

from BCF import images, BCFstitch, util
import glob

#%% 

HAL3_filename = r"C:\Users\r11403eb\OneDrive - The University of Manchester\SEM data\Edward Baker\HAL3\fullsamplemap\HAL3_map.bcf"
HAL3 = images(HAL3_filename,elements=['Fe','Mg','Si','Cl','Na'])
#MAC_fe = MAC.parseAndSlice('Fe')
HAL3.makeMaps()


#%%

stitch = BCFstitch(MAC_directory)

# %%
