
#%%
import os as os
import glob as glob
import numpy as np
import hyperspy.api as hs
import pandas as pd
from hyperspy.io_plugins.bruker import BCF_reader, HyperHeader, bcf_images, dictionarize, EDXSpectrum, bcf_hyperspectra
from hyperspy.misc.elements import elements as elements_db
import xml.etree.cElementTree as et
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
    Applies to signal class so must apply before parsing to element maps - this will be memory intensive. 
    '''


    def __init__(self,directory,elements=['Fe','Mg','Si','Cl','Na']):
        files = glob.glob(os.path.join(directory,'*.bcf'))
        self.elements = elements
        self.files = files
        self.directory = directory
        self.images = [images(i,elements=elements) for i in tqdm(files)]
        self.maps = [images(i,elements=elements).makeMaps() for i in tqdm(files)]
        print("Images loaded")
        self.resolutions = BCFstitch.getListOfResolutions(self)

    def getListOfResolutions(self):
        resolutions = [images(i).x_res for i in self.files]
        self.resolutions = resolutions
        self.max_res = min(resolutions)
        self.scale_factors = [i/self.max_res for i in resolutions]
        return resolutions

    def resampleImages(self):
        maps = self.maps
        for f in tqdm(range(len(self.files))):
            for ele in tqdm(range(len(self.elements))):
                maps[f][ele] = zoom(maps[f][ele],self.scale_factors[f],order=1)
        self.images_same_res_df = pd.DataFrame(maps,columns=self.elements)

    def makeBlankArea(self,debug=False):

        print("Making blank area")
        images = self.images
        x = np.array([i.x for i in images])
        y = np.array([i.y for i in images])
        self.x = x
        self.y = y

        nX_img = np.array([i.nX for i in images])
        nY_img = np.array([i.nY for i in images])

        x_res_img = np.array([i.x_res for i in images])
        y_res_img = np.array([i.y_res for i in images])

        x_dist = np.array([i.x_res*i.nX for i in images])
        y_dist = np.array([i.y_res*i.nY for i in images])
        self.x_dist = x_dist
        self.y_dist = y_dist

        bottom_lefts = x,y
        top_lefts = x,y+y_dist
        bottom_rights = x+x_dist,y

        self.left = self.x
        self.right = self.x + self.x_dist
        self.bottom = self.y
        self.top = self.y + self.y_dist

        x_min = min(bottom_lefts[0])
        x_max = max(bottom_rights[0])
        self.x_min = x_min
        self.x_max = x_max

        y_min = min(bottom_lefts[1])
        y_max = max(top_lefts[1])
        self.y_min = y_min
        self.y_max = y_max

        x_range = x_max - x_min
        y_range = y_max - y_min

        x_res = min([i.x_res for i in images])
        y_res = min([i.y_res for i in images])

        nX = int(x_range/x_res)
        nY = int(y_range/y_res)
        blank = np.zeros((nX,nY))

        plt.imshow(blank,cmap=mpl.colormaps['Greys'])
        plt.show()
        self.blank = blank

        if debug == True:
            plt.imshow(blank,cmap=mpl.colormaps['Greys'])
            plt.show()
            plt.scatter(x-min(x),y-min(y),color=['r','g'])
            for left,bottom,width,height,c in zip(x,y,x_dist,y_dist,['r','g']):
                print(width,height)
                rect=mpatches.Rectangle((left-min(x),bottom-min(y)),width,height,fill=False,color=c,linewidth=2)
                plt.gca().add_patch(rect)
            plt.show()

    def addToBlank(self,ele):
        print('Adding to blank')
        blank = self.blank
        print(np.shape(blank))
        for n,i in tqdm(enumerate(self.images)):
            print('---------------Next Step---------------')
            x = int((i.x - self.x_min)/i.x_res)#need to make a new nY and nX for the resampled images 
            y = int((i.y - self.y_min)/i.y_res)
            print('X and Y Coords:',x,y)
            nX_new,nY_new = np.shape(self.images_same_res_df[ele][n])
            print('Pixel numbers for x and y axes:',nX_new,nY_new)
            blank[x:x+nX_new,y:y+nY_new] = alh.images_same_res_df[ele][n]
            plt.imshow(alh.images_same_res_df[ele][n],cmap=mpl.colormaps['Greys'])
            plt.show()
            plt.imshow(blank,cmap=mpl.colormaps['Greys'])
            plt.show()
        self.composit = blank

#%%

bcf_dir = r'ExampleData'
alh = BCFstitch(bcf_dir)
alh.resampleImages()
alh.makeBlankArea(debug=True)
# alh.addToBlank('Si')

# stuff = []
# for i in glob.glob(os.path.join(bcf_dir,'*.bcf')):
#     print(i)
#     stuff.append(images(i,elements=['Fe','Mg','Si'],plot=True,save=True))



# %%


fe_fig, fe_ax = plt.subplots(frameon=False)
fe_ax.imshow(alh.maps[0][0],extent=[alh.left[0],alh.right[0],alh.bottom[0],alh.top[0]],cmap=mpl.colormaps['Greys'])                                        
fe_ax.imshow(alh.maps[1][0],extent=[alh.left[1],alh.right[1],alh.bottom[1],alh.top[1]],cmap=mpl.colormaps['Blues'])                                        
fe_ax.set_xlim([alh.x_min,alh.x_max])
fe_ax.set_ylim([alh.y_min,alh.y_max])
fe_fig.show()
# %%
