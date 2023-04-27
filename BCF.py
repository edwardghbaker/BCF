
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
from scipy.ndimage import zoom


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
        self.maps = maps
        self.maps_dir = {elements[i]: maps[i] for i in range(len(elements))}
    
class BCFstitch():
    '''
    need to use hyperspy rebinning function. Applies to signal class so must apply before parsing to element maps - this will be memory intensive. 
    '''


    def __init__(self,directory,elements=['Fe','Mg','Si']):
        files = glob.glob(os.path.join(directory,'*.bcf'))
        self.files = files
        self.directory = directory
        self.images = [images(i,elements=['Fe','Mg','Si']) for i in files]
        #need to make a list of filenames 
        self.resolutions = BCFstitch.getListOfResolutions(self)
        BCFstitch.resampleImages(self)
        BCFstitch.makeBlankArea(self)
        BCFstitch.addImagesToBlank(self)

    def getListOfResolutions(self):
        resolutions = [images(i).x_res for i in self.files]
        self.resolutions = resolutions
        self.max_res = min(resolutions)
        return resolutions

    def resampleImages(self):
        images = self.images
        max_res = self.max_res
        for i in tqdm(images):
            i.data = zoom(i.data,max_res/i.x_res,order=1)
        self.images_same_res = images

    def makeBlankArea(self):
        images = self.images_same_res
        x = [i.x for i in images]
        y = [i.y for i in images]
        self.x_min = min(x)
        self.x_max = max(x)
        self.y_min = min(y)
        self.y_max = max(y)
        x_range = max(x) - max(x)
        y_range = max(y) - min(y)
        x_res = images[0].x_res
        y_res = images[0].y_res
        nX = int(x_range/x_res)
        nY = int(y_range/y_res)
        blank = np.zeros((nY,nX))
        self.blank = blank

    def addImagesToBlank(self):
        images = self.images_same_res
        blank = self.blank
        for i in tqdm(images):
            x = int((i.x - self.x_min)/i.x_res)
            y = int((i.y - self.y_min)/i.y_res)
            blank[y:y+i.nY,x:x+i.nX] = i.data[:,:,0]
        self.composit = blank

