# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 12:56:49 2022

@author: Ed
"""
import hyperspy.api as hs
import matplotlib.pylab as pylab
pylab.rcParams['figure.figsize'] = 8, 6  # that's default image size for this interactive session
from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))
#from ipywidgets.widgets import interactive, fixed, interact
from skimage.exposure import equalize_hist, equalize_adapthist, rescale_intensity
from skimage import filters
import matplotlib.pyplot as plt
import os

StackName = 'Stack.bcf'
BinningFactor = 4

def MakeRGBFromMaps(R,G,B, bits=16, equalize=False, cutnoise=False):
    MaxPixelVal = (2**bits)-1
    
    # Make a Signal2D from the three images and copy over the metadata.
    x = hs.stack([R,G,B])
    x = hs.signals.Signal2D(x).T
    x.metadata = R.metadata
    x.original_metadata = R.original_metadata
    x.metadata.General.title = R.metadata.Sample.xray_lines[0] + '(R) ' + G.metadata.Sample.xray_lines[0] + '(G) ' + B.metadata.Sample.xray_lines[0] + '(B)'
    
    rdata = x.data[:,:,0]
    gdata = x.data[:,:,1]
    bdata = x.data[:,:,2]
            
    if cutnoise==True:
        NoiseCutAggressiveness = 1 # 0 means no cut, 1 means full cut.
        rdata[rdata < filters.threshold_triangle(rdata)*NoiseCutAggressiveness] = 0
        gdata[gdata < filters.threshold_triangle(gdata)*NoiseCutAggressiveness] = 0
        bdata[bdata < filters.threshold_triangle(bdata)*NoiseCutAggressiveness] = 0
 
    if equalize==True:
        rdata = equalize_hist(rdata)
        gdata = equalize_hist(gdata)
        bdata = equalize_hist(bdata)

    rdata = rescale_intensity(rdata)*MaxPixelVal
    gdata = rescale_intensity(gdata)*MaxPixelVal
    bdata = rescale_intensity(bdata)*MaxPixelVal

    x.data[:,:,0] = rdata
    x.data[:,:,1] = gdata
    x.data[:,:,2] = bdata

    # Convert to RGB.
    if bits==16:
        x.change_dtype('uint16')
        x.change_dtype('rgb16')
    else:
        x.change_dtype('uint8')
        x.change_dtype('rgb8')    
    
    # Copy the axes info over from the R map.
    x.axes_manager = R.axes_manager.copy()

    return(x)

def EqualizeAxis(fignum=None):
    # Because of a bug in Hyperspi, I sometimes have to set the axis to equal after we rebin.
    if fignum == None:
        fig=plt.gcf()
    else:
        fig = plt.figure(fignum)
    ax = fig.axes
    if type(ax) != list:
        ax.axis('equal')
    else:
        for a in ax:
            a.axis('equal')
    fig.tight_layout()
    
def RemovePtGaContamination(ElementDict):
    '''Make maps without contamination from Pt and Ga in FIB sections.
    The maps in ElementDict must contain Pt and Ga maps (obviously...)
    '''
    ElementDictNoPtGa = dict()
    for ElName, ElMap in ElementDict.items():
        if ElName in ['Pt', 'Ga']:
            ElementDictNoPtGa[ElName] = ElementDict[ElName].deepcopy()
            continue
        ReducedData = ElementDict[ElName].data.astype(float) - ElementDict['Pt'].data.astype(float) - ElementDict['Ga'].data.astype(float)
        ReducedData[ReducedData < 0] = 0
        ReducedElMap = ElMap.deepcopy()
        ReducedElMap.data = ReducedData
        ElementDictNoPtGa[ElName] = ReducedElMap
    return ElementDictNoPtGa

def PlotAndSaveRGBAndRGBNoPtGa(RGBElements, equalizeRGB=True, cutnoiseRGB=True, equalizeNoPtGa=True, cutnoiseNoPtGa=True):
    Img1 = MakeRGBFromMaps(ElementDict[RGBElements[0]],ElementDict[RGBElements[1]],ElementDict[RGBElements[2]], 
                      equalize=equalizeRGB, cutnoise=cutnoiseRGB)
    #Img1.plot()
    Img1.save(os.path.join(f'{BinningFactor}bin', 'Maps', RGBElements[0]+RGBElements[1]+RGBElements[2]+'.tif'), overwrite=True)
    Img2 = MakeRGBFromMaps(ElementDictNoPtGa[RGBElements[0]],ElementDictNoPtGa[RGBElements[1]],ElementDictNoPtGa[RGBElements[2]], 
                          equalize=equalizeNoPtGa, cutnoise=cutnoiseNoPtGa)
    #Img2.plot()
    Img2.save(os.path.join(f'{BinningFactor}bin', 'MapsNoPtGa', RGBElements[0]+RGBElements[1]+RGBElements[2]+'.tif'), overwrite=True)
    hs.plot.plot_images([Img1,Img2], colorbar=None)
    #plt.tight_layout()
    
    
    
    
HAADF, EDS = hs.load(StackName)
print(f'Stack dimensions are: {EDS.data.shape}.')
if BinningFactor > 1:
    EDS = EDS.rebin(scale=(BinningFactor, BinningFactor, 1))
    print(f'New stack dimensions are: {EDS.data.shape} with binning: {BinningFactor}')
EDS.plot()

HAADF.save(os.path.join(f'{BinningFactor}bin', 'Maps', 'HAADF.tif'), overwrite=True)
HAADF.save(os.path.join(f'{BinningFactor}bin', 'MapsNoPtGa', 'HAADF.tif'), overwrite=True)


print(HAADF.metadata)
print(EDS.metadata)




sumspec = EDS.sum()
sumspec.add_lines()
ElLines = ['C_Ka', 'N_Ka', 'O_Ka', 'Na_Ka', 'Pt_La', 'Mg_Ka', 'Al_Ka', 'Si_Ka', 'P_Ka', 'S_Ka', 'Cl_Ka', 'K_Ka', 'Ca_Ka', 'Ti_Ka', 'Cr_Ka', 'Mn_Ka', 'Fe_Ka', 'Ni_Ka', 'Zn_Ka', 'Cu_Ka', 'Ga_Ka' ]
Elements = EDS.get_lines_intensity(ElLines)
Elements = [El.T for El in Elements]
ElementDict = dict()
for El in Elements:
    El.change_dtype('uint16')
    El.save(os.path.join(f'{BinningFactor}bin', 'Maps', El.metadata.Sample.xray_lines[0]+'.tif'), overwrite=True)
    ElementDict[El.metadata.Sample.elements[0]] = El
    #print(f'{El.metadata.Sample.elements[0]}')
    
ElementDictNoPtGa = RemovePtGaContamination(ElementDict)
ElementsNoPtGa = [El for El in ElementDictNoPtGa.values()]
for El in Elements:
    El.save(os.path.join(f'{BinningFactor}bin', 'MapsNoPtGa', El.metadata.Sample.xray_lines[0]+'.tif'), overwrite=True)
fig = plt.figure(figsize=(20,20))
hs.plot.plot_images(Elements, per_row=ceil(sqrt(len(Elements))), colorbar=None, axes_decor='off', fig=fig)
plt.show()
fig.savefig(os.path.join(f'{BinningFactor}bin', 'Maps', 'Mosaic.png'))