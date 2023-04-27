#%% Imports

from BCF import images, BCFstitch, util
import glob

#%% 

MAC_filename = r'C:\Users\r11403eb\OneDrive - The University of Manchester\meteoriteData\Meteorites\MAC 88136\MAC88136_position4_001.bcf'
MAC_filename = r"C:\Users\r11403eb\OneDrive - The University of Manchester\meteoriteData\Meteorites\MAC 88136\MAC88136_position4_001.bcf"
MAC_directory = r"C:\Users\r11403eb\OneDrive - The University of Manchester\meteoriteData\Meteorites\MAC 88136"

MAC = images(MAC_filename)
#MAC_fe = MAC.parseAndSlice('Fe')
MAC.makeMaps()


#%%

stitch = BCFstitch(MAC_directory)

# %%
