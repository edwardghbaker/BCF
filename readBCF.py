# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 17:45:47 2021

@author: r11403eb
"""
import os as os
import glob as glob
import numpy as np
import hyperspy.api as hs
from hyperspy.io_plugins.bruker import BCF_reader, HyperHeader, bcf_images, dictionarize, EDXSpectrum, bcf_hyperspectra
from hyperspy.misc.elements import elements as elements_db
#from hyperspy.io_plugins.EDXSpectrum import EDXSpec
#from lxml import etree, objectify  # not xml.etree; 
import xml.etree.cElementTree as et
import matplotlib.pyplot as plt
from tqdm import tqdm
#import scipy.misc.imresize as imresize
import PIL.Image as Image 
import pyresample as re


# #%% define the functions
# class Bruker:
    
#     def plotWithInfo(directory = r'C:\Users\r11403eb\OneDrive - The University of Manchester\meteoriteData\LAR 12156',
#                      file = 'LAR12156_position9_001.bcf',
#                      MinEle=['Cl','Fe','Si','S','Na']):
        
#         newdir = os.path.join(directory, file.split(".", 1)[0])
#         os.makedirs(newdir,exist_ok=True)
        
#         filename = str(directory+'\\'+file)
#         data = hs.load(filename, select_type='None', lazy=True)
#         signal = data[-1].as_lazy()
#         signal.add_elements(MinEle)
#         elements = signal.metadata.Sample.elements
#         signal.add_lines()
#         print(elements)
        
#         hs.plot.plot_images(data[0].as_lazy())
#         plt.savefig(str(newdir+'\\'+'BSE.PNG'),dpi=1200)
#         plt.close()
#         signal.sum().plot(True)
#         plt.savefig(str(newdir+'\\'+'Spectrum.PNG'),dpi=1200)
#         plt.close()

#         EDSdata = signal.get_lines_intensity()
#         for i in tqdm(range(len(EDSdata))):
#             EDSdata[i].plot()
#             plt.savefig(str(newdir+'\\'+elements[i]+'_info'),dpi=1200)
            
#     def plotNoInfo(directory=r'C:\Users\r11403eb\OneDrive - The University of Manchester\meteoriteData\LAR 12156',
#                    file='LAR12156_position9_001.bcf',
#                    MinEle=['Cl','Fe','Si','S','Na']):
        
#         newdir = os.path.join(directory, file.split(".", 1)[0])
#         os.makedirs(newdir,exist_ok=True)

#         filename = str(directory+'\\'+file)
#         data = hs.load(filename, lazy=True)
#         data[-1].add_elements(MinEle)
#         elements = data[-1].metadata.Sample.elements
#         data[-1].add_lines()
#         print(elements)
#         EDSdata = data[-1].get_lines_intensity()
        
#         for i in range(len(EDSdata)):
#             if elements[i] == 'Cl':
#                 CMAP = 'Greens'
#             elif elements[i] == 'Fe':
#                 CMAP = 'Blues'
#             elif elements[i] == 'Si':
#                 CMAP = 'Reds'
#             else:
#                 CMAP = 'Oranges'
            
#             # hs.plot.plot_images(hs.transpose(EDSdata[i]),cmap=CMAP, colorbar=None,label=None,tight_layout=True,padding={'wspace':0.0, 'hspace':0.0})
#             # plt.axis('off') 
#             # Bruker.plotNoInfo(directory=folder,file=i)
            
#             # elif method == 'Info':
        
#             #     for i in tqdm(listOfFiles):
#             #         Bruker.plotWithInfo(directory=folder,file=i)

#%%   

ALH_filename = r"C:\Users\r11403eb\OneDrive - The University of Manchester\meteoriteData\Meteorites\ALH 84206.9\SEM\ALH.bcf"
#Fe_ka = elements_db['Fe']['Atomic_properties']['Xray_lines']['Ka']['energy (keV)']

MAC_filename = r'C:\Users\r11403eb\OneDrive - The University of Manchester\meteoriteData\Meteorites\MAC 88136\MAC88136_position4_001.bcf'

# MAC_header = BCF_reader(MAC_filename).header
# MAC_data = BCF_reader(MAC_filename).parse_hypermap(lazy=False)
# MAC_metaData = bcf_images(BCF_reader(MAC_filename))
# #MAC = 
# fuckyou = bcf_images(BCF_reader(MAC_filename))[0]['original_metadata']['Stage']
# FU = bcf_images(BCF_reader(MAC_filename))

# MAC_header.get_consistent_min_channels()

# avMAC_EDS = np.average(MAC_data,(0,1))

# metaData = MAC_header.get_spectra_metadata()


def EnergyKa(ele):
    return hs.material.elements[ele].Atomic_properties.Xray_lines['Ka']['energy (keV)']

def makeMap_fromBCF(filename=r'C:\Users\r11403eb\OneDrive - The University of Manchester\meteoriteData\Meteorites\MAC 88136\MAC88136_position4_001.bcf',
                    ele=['Si','Fe']
                    ):
    
    header = BCF_reader(filename).header
    data = BCF_reader(filename).parse_hypermap(lazy=True)
    fig,ax=plt.subplots(len(ele),1,squeeze=True)
    
    for i,e in enumerate(ele):
        print(e)
        channel = header.get_spectra_metadata().energy_to_channel(EnergyKa(e))
        ax[i].imshow(data[:,:,channel],cmap='Pastel1')

    
makeMap_fromBCF(filename=MAC_filename)



# import cv2
 
# img = cv2.imread('/home/img/python.png', cv2.IMREAD_UNCHANGED)
 
# print('Original Dimensions : ',img.shape)
 
# scale_percent = 220 # percent of original size
# width = int(img.shape[1] * scale_percent / 100)
# height = int(img.shape[0] * scale_percent / 100)
# dim = (width, height)
  
# # resize image
# resized = cv2.resize(img, dim, interpolation = cv2.INTER_AREA)
 
# print('Resized Dimensions : ',resized.shape)
 
# cv2.imshow("Resized image", resized)
# cv2.waitKey(0)
# cv2.destroyAllWindows()



#%%


# filename = r"C:\Users\r11403eb\OneDrive - The University of Manchester\meteoriteData\MAC 88136\MAC88136_position_001.bcf"
# MinEle=['Cl','Fe','Si','S','Na']
# data = hs.load(filename, select_type='None', lazy=True)
# signal = data[-1]
# signal.add_elements(MinEle)
# signal.add_lines()
# # x = signal.get_lines_intensity()

# # hs.plot.plot_images(x, axes_decor='off', scalebar='all', lazy=True)
# EDSdata = signal.get_lines_intensity()
# for i in tqdm(range(len(EDSdata))):
#     EDSdata[i].plot()
#     plt.savefig(str('newdir'+'\\'+MinEle[i]+'_info'),dpi=1200)

# makePics(folder=r'C:/Users/r11403eb/OneDrive - The University of Manchester/meteoriteData/SEM/ALH')



#%% Old code to extract header info
# # here we get the single file system handler
# #all file items in this singlefilesystem class instance is held inside
# # dictionary hierarchy, we fetch the header:
# header = file.vfs['EDSDatabase']['HeaderData']  
# #the items in the tree have special method which allows to get selection as BytesIO object:
# bstring = header.get_as_BytesIO_string()
# #rewind:
# bstring.seek(0)
# #for very huge nodes:
# parser = objectify.makeparser(huge_tree=True)
# #the final steps to print the xml to file for examination in text editor:
# xml_thingy = etree.fromstring(bstring.read(), parser)
# xml_root = xml_thingy.getroottree()
# xml_root.write('header_output.xml', pretty_print=True)
# # without pretty_print everything would be in the single line