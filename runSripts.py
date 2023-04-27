#%% Imports

from BCF import images, BCFstitch, util
import matplotlib.pyplot as plt
import glob

#%% 

HAL3 = r"C:\Users\User\OneDrive - The University of Manchester\SEM data\Edward Baker\HAL3\fullsamplemap\HAL3_map.bcf"

HAL3_image = images(HAL3,elements=['Fe','Mg','Si','Cl','Na'])
#MAC_fe = MAC.parseAndSlice('Fe')
HAL3_image.makeMaps()


#%%
fig, ax = plt.subplots(1,3)
ax[0].imshow(HAL3_image.maps['Fe'])
ax[1].imshow(HAL3_image.maps['Mg'])
ax[2].imshow(HAL3_image.maps['Si'])

# %%
