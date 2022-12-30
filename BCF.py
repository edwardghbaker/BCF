
#%%
import os as os
import glob as glob
import numpy as np
import hyperspy.api as hs
from hyperspy.io_plugins.bruker import BCF_reader, HyperHeader, bcf_images, dictionarize, EDXSpectrum, bcf_hyperspectra
from hyperspy.misc.elements import elements as elements_db
import xml.etree.cElementTree as et
import matplotlib.pyplot as plt
from tqdm import tqdm
import PIL.Image as Image 




#%% 

class util():
    
    def energyKa(ele):
        return hs.material.elements[ele].Atomic_properties.Xray_lines['Ka']['energy (keV)']


class images():

    def __init__(self,filename, elements=['Fe','Mg','Si']):
        self.data = BCF_reader(filename).parse_hypermap(lazy=True)
        self.data.rechunk({0: -1, 1: -1, 2: 1})
        
        header = BCF_reader(filename).header
        self.header = header
        
        self.x = header.stage_metadata['X']
        self.y = header.stage_metadata['Y']
        self.z = header.stage_metadata['Z']
        
        self.nX = header.dsp_metadata['ImageWidth']
        self.nY = header.dsp_metadata['ImageHeight']
        
        self.x_res = header.x_res
        self.y_res = header.y_res
        
        self.elements = elements
        self.channels = {i: header.get_spectra_metadata().energy_to_channel(util.energyKa(i)) for i in elements}

    def parseAndSlice(self,element):
        return self.data[:,:,self.channels[element]]

    def makeMaps(self):
        elements = self.elements
        maps = [plt.imshow(self.parseAndSlice(i).compute()) for i in elements]
        self.maps = {elements[i]: maps[i] for i in range(len(elements))}
    
    
class BCFstitch():

    def __init__(self,directory,elements=['Fe','Mg','Si']):
        files = glob.glob(os.path.join(directory,'*.bcf'))
        self.files = files
        self.directory = directory
        self.images = [images(i,elements=['Fe','Mg','Si']) for i in files]
        #need to make a list of filenames 
        self.getListOfResolutions(self)
        self.resampleImages(self)

        
    def getListOfResolutions(self):
        resolutions = [images(i).x_res for i in self.files]
        self.resolutions = resolutions
        self.max_res = min(resolutions)
        return resolutions

    def resampleImages(self):
        images = self.images
        max_res = self.max_res
        for i in images:
            i.data = i.data.rebin(scale=(i.x_res/max_res,i.y_res/max_res,1))
        self.images = images

    def makeStitch(self):
        images = self.images
        max_res = self.max_res
        nX = max([i.nX for i in images])
        nY = max([i.nY for i in images])
        x = max([i.x for i in images])
        y = max([i.y for i in images])
        z = max([i.z for i in images])
        nZ = len(images[0].elements)
        data = np.zeros((nX,nY,nZ))
        for i in images:
            data[i.x:i.x+i.nX,i.y:i.y+i.nY,:] = i.data
        return hs.signals.Signal2D(data)

#    def loadImages()
        # call image class to get list of dicts?
        # define highest res and the dimensions of the map
        
#    def resampleImages():
        #take image from list and convert to new image with correct resolution 
        
#    def makeStitch():
        #run functions for each element and save png of image with blank as backspace
    



#%%



MAC_filename = r'C:\Users\User\OneDrive - The University of Manchester\meteoriteData\Meteorites\MAC 88136\MAC88136_position4_001.bcf'

MAC_filename = r"C:\Users\User\OneDrive - The University of Manchester\meteoriteData\Meteorites\MAC 88136\MAC88136_position4_001.bcf"

MAC_directory = r"C:\Users\User\OneDrive - The University of Manchester\meteoriteData\Meteorites\MAC 88136"

MAC = images(MAC_filename)
MAC_fe = MAC.parseAndSlice('Fe')


stitch = BCFstitch(MAC_directory)
stitch.makeStitch().plot()

# %%
