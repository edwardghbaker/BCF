#%% Imports

from BCF import images, BCFstitch, util
import matplotlib.pyplot as plt
import glob

#%% 


HAL3 = r"C:\Users\r11403eb\OneDrive - The University of Manchester\SEM data\Edward Baker\SAH97072_position13_001_1.bcf"

HAL3_image = images(HAL3,elements=['Fe'])
maps = HAL3_image.makeMaps()



#%%
fig, ax = plt.subplots(1,3)
ax[0].imshow(HAL3_image.maps['Fe'])
ax[1].imshow(HAL3_image.maps['Mg'])
ax[2].imshow(HAL3_image.maps['Si'])

# %%
