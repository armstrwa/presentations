## --------------------------- ##
#  GEOSPATIAL PYTHON TUTORIAL   #
#  William Armstrong		    #
#  20 Sep 2016					#
## --------------------------- ##

# IMPORT MODULES
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from shapely import geometry
import fiona
import os
from descartes import PolygonPatch
from skimage import exposure

#import fiona
#from shapely import geometry
#shp = geometry.LinearRing([(0,0),(1,0),(0,1)])
#shp.centroid.coords.xy

# DATA LOCATIONS
rgiFn = '/Users/wiar9509/Documents/CU2016-2017/presentations/geospatialPython/01_rgi50_Alaska/01_rgi50_Alaska.shp'
stateFn = '/Users/wiar9509/Documents/CU2016-2017/presentations/geospatialPython/cb_2015_us_state_20m/cb_2015_us_state_20m.shp'
rasterFn = '/Users/wiar9509/Documents/CU2016-2017/presentations/geospatialPython/LC80330322016262LGN00/LC80330322016262LGN00_B8.TIF'

# Lists, tuples, dicts
a = [ 1, 4, 8, 9]
a.append(3)
a.sort()
a.reverse()

t = ( (1,2) , (5,9) )
l = [ 1, 2 , 3 , 4]
d = {'age':( (12,41,19) ), 'name':['brian','claire','joe']}

# FOR LOOPS
for i in range(0,10):
	print i

# OPEN DATA

c = fiona.open(stateFn,'r') # open collection stateFn wirh 'r' for 'read'
#with fiona.open(stateFn,'r') as c:
# Initialize
cxList = []
cyList = []
areaList = []

#for i in range(0,len(c)): # loop over all features in collection
for f in c: # loop over all features in collection
	#print f['properties']['NAME'], f['properties']['STUSPS']
	#f = c.next()
	tmpShp = geometry.shape(f['geometry'])
	tmpX,tmpY = tmpShp.centroid.xy
	cxList.append(tmpX)
	cyList.append(tmpY)
	areaList.append(tmpShp.area)

np.where(areaList==np.max(areaList))
cxList=np.delete(cxList,13)
cyList=np.delete(cyList,13)
areaList=np.delete(areaList,13)
plt.scatter(cxList,cyList,s=areaList*10)
plt.axis('equal')
plt.show()	
plt.close()

# Turn a state's boundaries into a shapely geometry
f = c[0]
#g = f['geometry']['coordinates'][0]
#shp = geometry.Polygon(g)
shp = geometry.shape(f['geometry'])
cent = shp.centroid
hull = shp.convex_hull
buff = shp.buffer(0.5)
ext = shp.exterior
x,y = ext.xy
cx,cy = cent.xy	
hx,hy = hull.exterior.xy

plt.plot(x,y,c='b',lw=2,label='outline')
plt.plot(hx,hy,c='c',lw=2,label='convex hull')
plt.scatter(cx,cy,s=40,marker='o',facecolor='r',label='centroid')
plt.xlabel('Longitude [deg]',fontsize=16)
plt.ylabel('Latitude [deg]',fontsize=16)
plt.legend(loc=1,frameon=False)
plt.axis('tight')
plt.savefig('texas_shapely.png',r=300)
plt.show()
plt.close()

# simple dummy data
square = geometry.Polygon([(0,0),(1,0),(1,1),(0,1)])
triangle = geometry.Polygon([(0.5,0.5),(1.5,0.5),(1,1.5)])
square.intersects(triangle) # Returns true if overlap
isect = square.intersection(triangle) # 'overlap'
bsect = square.intersection(triangle.boundary)
diff = square.difference(triangle) # 'subtracts'
uni = triangle.union(square) # 'combines'

triPatch = PolygonPatch(triangle, fc='b', ec='k', alpha=0.5, zorder=2)
sqPatch = PolygonPatch(square, fc='r', ec='k', alpha=0.5, zorder=2)
iPatch = PolygonPatch(isect, fc='g', ec='k', alpha=0.9, zorder=2,hatch='/')
dPatch = PolygonPatch(dsect, fc='none', ec='k', alpha=0.9, zorder=2,hatch='-')
#uPatch = PolygonPatch(dsect, fc='none', ec='r', lw=3,alpha=0.9, zorder=2)

fig = plt.figure(0)
ax = fig.add_subplot(111)
ax.add_patch(triPatch)
ax.add_patch(sqPatch)
ax.add_patch(iPatch)
ax.add_patch(dPatch)
ax.axis('equal')
ax.set_ylim([-0.5,2])
ax.set_xlim([-0.5,2])
ax.axis('off')   
fig.savefig('test.png')
fig.show()
fig.clf()


### RASTER DATA
# Also getting segmentation fault, but not if just import rasterio & matplotlib
src = rasterio.open(rasterFn) # open raster
ll = src.index(472996,4426147) # approx lower left for boulder subset
ur = src.index(484804,4434040) # approx upper right for boulder subset
#w = ( (ur[1],ll[1]) , (ll[0],ur[0]) )
w = ( (ur[0],ll[0]) , (ll[1],ur[1]) )
sub = src.read(1,window=w)
bLow,bHigh = np.percentile(sub,[1,98])
lims = np.percentile(sub,[0,99.5])
hist,edges = np.histogram(sub,bins=500)
hist = np.insert(hist,0,0)
plt.plot(edges,hist.astype('float')/np.sum(hist),lw=2)
#plt.plot((bLow,bLow),(0,hist.max().astype('float')/np.sum(hist)),ls='--',c='k')
plt.plot((bLow,bLow),(0,0.035),ls='--',c='k')
#plt.plot((bHigh,bHigh),(0,hist.max().astype('float')/np.sum(hist)),ls='--',c='k',label=( ('1%-98%') ))
plt.plot((bHigh,bHigh),(0,0.035),ls='--',c='k',label=( ('1%-98%') ))
plt.xlabel('Brightness [-]',fontsize=16)
plt.ylabel('P(Brightness)',fontsize=16)
plt.legend(loc=1,frameon=False)
plt.xlim(lims)
plt.ylim((0,0.035))
plt.savefig('histogram.pdf')
plt.show()
plt.close()

xmin,ymax = src.affine * (0,0)
xmax,ymin = src.affine * (src.height,src.width)
plt.imshow(sub,cmap='gray',origin='upper',clim=((bLow,bHigh)),extent=((xmin,xmax,ymin,ymax)),aspect=1)
plt.xlabel('Easting [m]',fontsize=16)
plt.ylabel('Northing [m]',fontsize=16)
plt.savefig('boulder.png',r=600)
plt.show()
plt.close()

# WRITING RASTER DATA

stretch=np.copy(sub)
stretch = exposure.rescale_intensity(stretch,(bLow,bHigh))
new_dataset = rasterio.open('boulderSub.tif', 'w', driver='GTiff',
                            height=sub.shape[0], width=sub.shape[1],
                            count=1, dtype=sub.dtype,
                            crs=src.crs, transform=src.transform)
new_dataset.write(stretch,1)                            
new_dataset.close()
