# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# createFindingChartsFromDECamFITSImage:
#   creates finding charts from a list of FITS images
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def createFindingChartsFromDECamFITSImage(image_file, fitsType, target_radeg, target_decdeg, target_id, target_name, outputDir):
    
    import numpy as np
    import pandas as pd
    from scipy import interpolate
    import glob
    import math
    import os
    from shutil import copy2
    from datetime import date
    
    from matplotlib import path
    from matplotlib.patches import Arrow
    import matplotlib.pyplot as plt

    import matplotlib.transforms as mtransforms

    from astropy.io import fits
    import astropy.coordinates as coord
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    from astropy.visualization import SqrtStretch
    from astropy.visualization import LinearStretch
    from astropy.visualization import LogStretch

    from astropy.visualization.mpl_normalize import ImageNormalize
    from astropy.visualization import ZScaleInterval
    from astropy.visualization import MinMaxInterval

    from astropy.wcs import WCS
    from astropy.wcs.utils import pixel_to_skycoord
    from astropy.wcs.utils import skycoord_to_pixel

    from astropy.nddata import Cutout2D

    
    # Read in image...
    if fitsType == "DECam":
        # DECam FITS image...
        image_data = fits.getdata(image_file, ext=1)
        image_hdr = fits.getheader(image_file, ext=1)
        arcsec_per_pixel = 0.263
    elif fitsType == "SkyView":
        # Skyview FITS image
        image_data = fits.getdata(image_file)
        image_hdr = fits.getheader(image_file)
        arcsec_per_pixel = 3600.*abs(image_hdr['CDELT1'])
    else:
        # Guess...
        image_data = fits.getdata(image_file)
        image_hdr = fits.getheader(image_file)
        arcsec_per_pixel = 3600.*abs(image_hdr['CDELT1'])


    # We want zscale scaling...
    interval = ZScaleInterval()
    vlimits = interval.get_limits(image_data)
    # We want black stars on white background, so we switch the vmin/vmax limits...
    norm = ImageNormalize(vmin=vlimits[1], vmax=vlimits[0], stretch=LinearStretch())
    
    # Extract WCS from image header, and create a SkyCoord object for the target's coordinates...
    wcs = WCS(image_hdr)
    position = SkyCoord(target_radeg,target_decdeg, unit="deg", frame='icrs')
    
    # The arrow length for the North and East direction arrows...
    arrowLength = 60./arcsec_per_pixel  # 1 arcmin arrow length
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # CCD size finding chart:
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Calculate the x,y position of the target in pixel coordinates...
    xy = skycoord_to_pixel(position, wcs)
    x = xy[0]
    y = xy[1]

    # Create figure...
    fig = plt.figure()
    ax=plt.subplot(1,1,1)
    plt.axis('off')
    
    # DECam FITS image...
    if fitsType == "DECam":

        # To get image to plot at the preferred orientiation (N at top, E at the right)...
        #  we transpose the image_data and set the origin of the plot as "upper"...
        im = ax.imshow(image_data.T, cmap='gray', origin='upper', norm=norm)

        # Circle the target (don't forget to transpose the x,y coordinates!)...
        ax.scatter(y,x, s=30, edgecolor='blue', facecolor='none')

        # Add title...
        title = """ObjID: %d  \n(Name: %s)""" % (target_id, target_name)
        plt.title(title, color='blue', fontsize=10) 

        # Add "compass rose"...
        a =  Arrow(200., 100., 0., arrowLength, edgecolor='black', facecolor='none')
        ax.add_patch(a)
        plt.text(200,100,'N',color='black')

        a =  Arrow(200., 100.+arrowLength, arrowLength, 0., edgecolor='black', facecolor='none')
        ax.add_patch(a)
        plt.text(200+arrowLength,100+arrowLength,'E',color='black')
        plt.text(50,25+0.5*arrowLength,'1.0\'',color='black',rotation=90)

        # Include color bar at the bottom of the plot...
        #fig.colorbar(im, orientation="horizontal")
    
    # Skyview FITS image
    elif fitsType == "SkyView":

        im = ax.imshow(image_data, cmap='gray', origin='lower', norm=norm)        

        # Circle the target...
        ax.scatter(x,y, s=30, edgecolor='blue', facecolor='none')

        # Add title...
        title = """ObjID: %d  \n(Name: %s)""" % (target_id, target_name)
        plt.title(title, color='blue', fontsize=10) 

        # Add "compass rose"...
        a =  Arrow(200., 100., 0., arrowLength, edgecolor='black', facecolor='none')
        ax.add_patch(a)
        plt.text(200,100,'N',color='black')

        a =  Arrow(200., 100., arrowLength, 0., edgecolor='black', facecolor='none')
        ax.add_patch(a)
        plt.text(200+arrowLength,100,'W',color='black')
        plt.text(50,25+0.5*arrowLength,'1.0\'',color='black',rotation=90)

        # Include color bar at the bottom of the plot...
        #fig.colorbar(im, orientation="horizontal")
    

    # Save figure
    outputFile = """objid_%d.%s.ccd_image.png""" % (target_id, target_name)
    outputFile = os.path.join(outputDir,outputFile)
    fig.savefig(outputFile, format='png')
    
    imageNameCCD = outputFile
    
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # 5arcmin x 5arcmin finding chart:
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Create a 5arcmin x 5arcmin "Quantity" object...
    size = u.Quantity((5.0, 5.0), u.arcmin)

    # Create a 5arcmin x 5arcmin cutout from the full image data, 
    #  centered on the target ra, dec
    cutout = Cutout2D(image_data, position, size, mode='partial',wcs=wcs)

    # Calculate the x,y position of the target in pixel coordinates 
    #  within the 5arcmin x 5arcmin cutout...
    xy = skycoord_to_pixel(position, cutout.wcs)
    x = xy[0]
    y = xy[1]

    # Create figure...
    fig = plt.figure()
    ax=plt.subplot(1,1,1)
    plt.axis('off')

    # DECam FITS image...
    if fitsType == "DECam":

        # To get image to plot at the preferred orientiation (N at top, E at the right)...
        ax.imshow(cutout.data.T, cmap='gray', origin='upper', norm=norm)

        # Circle the target (don't forget to transpose the x,y coordinates!)...
        ax.scatter(y,x, s=30, edgecolor='blue', facecolor='none')

        # Add title...
        text1 = """ObjID: %d""" % (target_id)
        text2 = """(Name:  %s)""" % (target_name)
        plt.text(0.5*size.value[1]*60./arcsec_per_pixel,50,text1,color='blue',fontsize=10)
        plt.text(0.5*size.value[1]*60./arcsec_per_pixel,50+0.05*size.value[0]*60./arcsec_per_pixel,text2,color='blue',fontsize=10)

        # Add "compass rose"...
        a =  Arrow(50., 50., 0., arrowLength, edgecolor='black', facecolor='none')
        ax.add_patch(a)
        plt.text(50,50,'N',color='black')

        a =  Arrow(50., 50.+arrowLength, arrowLength, 0., edgecolor='black', facecolor='none')
        ax.add_patch(a)
        plt.text(50+arrowLength,50+arrowLength,'E',color='black')
        plt.text(-10,15+0.5*arrowLength,'1.0\'',color='black',rotation=90)
    
        # Include color bar at the side of the plot...
        #fig.colorbar(im)

    # Skyview FITS image
    elif fitsType == "SkyView":

        im = ax.imshow(cutout.data, cmap='gray', origin='lower', norm=norm)

        # Circle the target...
        ax.scatter(x,y, s=30, edgecolor='blue', facecolor='none')

        # Add title...
        text1 = """ObjID: %d""" % (target_id)
        text2 = """(Name:  %s)""" % (target_name)
        plt.text(0.5*size.value[1]*60./arcsec_per_pixel,50,text1,color='blue',fontsize=10)
        plt.text(0.5*size.value[1]*60./arcsec_per_pixel,50+0.05*size.value[0]*60./arcsec_per_pixel,text2,color='blue',fontsize=10)

        # Add "compass rose"...
        a =  Arrow(50., 50., 0., arrowLength, edgecolor='black', facecolor='none')
        ax.add_patch(a)
        plt.text(50,50+arrowLength,'N',color='black')

        a =  Arrow(50., 50., arrowLength, 0., edgecolor='black', facecolor='none')
        ax.add_patch(a)
        plt.text(50+arrowLength,50,'W',color='black')
        plt.text(30,60+0.5*arrowLength,'1.0\'',color='black',rotation=90)
    
        # Include color bar at the side of the plot...
        #fig.colorbar(im)
    
    # Save figure
    outputFile = """objid_%d.%s.5x5arcmin.png""" % (target_id, target_name)
    outputFile = os.path.join(outputDir,outputFile)
    fig.savefig(outputFile, format='png')
    
    imageName5x5 = outputFile
    
    return (imageNameCCD,imageName5x5)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

