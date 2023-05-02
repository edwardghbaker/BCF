
#%%
import os as os
import glob as glob
import numpy as np
import hyperspy.api as hs
from hyperspy.io_plugins.bruker import BCF_reader, HyperHeader, bcf_images, dictionarize, EDXSpectrum, bcf_hyperspectra
from hyperspy.misc.elements import elements as elements_db
import xml.etree.cElementTree as et
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm
import PIL.Image as Image 
from scipy.ndimage import zoom


#%% 

class util():
    
    def energyKa(ele):
        return hs.material.elements[ele].Atomic_properties.Xray_lines['Ka']['energy (keV)']

class images():

    def __init__(self,filename, elements=['Fe','Mg','Si'],plot=False,save=False):
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

        self.plot = plot
        self.save = save
        self.filename = filename

    def parseAndSlice(self,element):
        return self.data[:,:,self.channels[element]]

    def makeMaps(self):
        elements = self.elements
        maps = [self.parseAndSlice(i).compute() for i in elements]  
        self.maps = maps
        self.maps_dir = {elements[i]: maps[i] for i in range(len(elements))}

        if self.plot == True:
            self.plotMaps()
        if self.save == True:
            self.saveMaps()
            
        return maps
    
    def saveMaps(self):
        maps = self.maps
        elements = self.elements
        for i,m in enumerate(maps):
            fig, ax = plt.subplots()
            ax.imshow(m)
            os.makedirs(f"{os.path.dirname(self.filename)}\\{os.path.basename(self.filename).split('.')[0]}", exist_ok=True)
            fig.savefig(f"{os.path.dirname(self.filename)}\\{os.path.basename(self.filename).split('.')[0]}\\{elements[i]}.png")
    
    def plotMaps(self):
        maps = self.maps
        elements = self.elements
        for i,m in enumerate(maps):
            fig, ax = plt.subplots()
            ax.imshow(m)
            plt.show()
    

#%%
class BCFstitch():
    '''
    need to use hyperspy rebinning function. Applies to signal class so must apply before parsing to element maps - this will be memory intensive. 
    '''


    def __init__(self,directory,elements=['Fe','Mg','Si','Cl','Na']):
        files = glob.glob(os.path.join(directory,'*.bcf'))
        self.elements = elements
        self.files = files
        self.directory = directory
        self.maps = [images(i,elements=elements).makeMaps() for i in files]
        print("Images loaded")
        #print(np.shape(self.maps[0][0]))
        self.resolutions = BCFstitch.getListOfResolutions(self)


    def getListOfResolutions(self):
        resolutions = [images(i).x_res for i in self.files]
        self.resolutions = resolutions
        self.max_res = min(resolutions)
        self.scale_factors = [i/self.max_res for i in resolutions]
        return resolutions

    def resampleImages(self):
        maps = self.maps
        max_res = self.max_res
        for f in tqdm(len(files)):
            print(np.shape(i.data))
            i.data = zoom(i.data,i.x_res/max_res,order=1)
            print(np.shape(i.data))
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

#%%

bcf_dir = r'C:\Users\r11403eb\OneDrive - The University of Manchester\SEM data\Edward Baker\2021.11.04'
alh = BCFstitch(bcf_dir)
alh.resampleImages()
# stuff = []
# for i in glob.glob(os.path.join(bcf_dir,'*.bcf')):
#     print(i)
#     stuff.append(images(i,elements=['Fe','Mg','Si'],plot=True,save=True))



# %%
