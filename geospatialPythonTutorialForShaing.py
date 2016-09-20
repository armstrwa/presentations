## --------------------------- ##
#  GEOSPATIAL PYTHON TUTORIAL   #
#  William Armstrong		    #
#  20 Sep 2016					#
## --------------------------- ##

### IMPORT MODULES ###
import matplotlib.pyplot as plt
import numpy as np
from shapely import geometry
import fiona
import rasterio
import os
from descartes import PolygonPatch
from skimage import exposure


### DATA LOCATIONS ###
rgiFn = '/Users/wiar9509/Documents/CU2016-2017/presentations/geospatialPython/01_rgi50_Alaska/01_rgi50_Alaska.shp'
stateFn = '/Users/wiar9509/Documents/CU2016-2017/presentations/geospatialPython/cb_2015_us_state_20m/cb_2015_us_state_20m.shp'
rasterFn = '/Users/wiar9509/Documents/CU2016-2017/presentations/geospatialPython/LC80330322016262LGN00/LC80330322016262LGN00_B8.TIF'


### PYTHON BASICS ###

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

### VECTOR DATA ###
# OPEN DATA
c = fiona.open(stateFn,'r') # open collection stateFn wirh 'r' for 'read'

# Initialize empty lists
cxList = []
cyList = []
areaList = []

for f in c: # loop over all features in collection
	print f['properties']['NAME'], f['properties']['STUSPS'] # print name and abbrev
	tmpShp = geometry.shape(f['geometry']) # get geometry
	tmpX,tmpY = tmpShp.centroid.xy # get centroid coords
	cxList.append(tmpX) # append to list
	cyList.append(tmpY)
	areaList.append(tmpShp.area) # append area to list

# Plotting, first bit deletes Alaska to get better detail in CONUS
np.where(areaList==np.max(areaList))
cxList=np.delete(cxList,13)
cyList=np.delete(cyList,13)
areaList=np.delete(areaList,13)
plt.scatter(cxList,cyList,s=areaList*10)
plt.axis('equal')
plt.show()	
plt.close()

## MANIPULATING A SINGLE GEOMETRY
# Turn a state's boundaries into a shapely geometry
f = c[0] # get first feature
shp = geometry.shape(f['geometry']) # turn it into a geometry
cent = shp.centroid # get its centroid
hull = shp.convex_hull # get convex hull
ext = shp.exterior # get exterior ring
x,y = ext.xy # get coords of ring
cx,cy = cent.xy	 # get coords of centroid
hx,hy = hull.exterior.xy # get coords of convex hull

# plot
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

# MAKING NEW GEOMETRY
# simple dummy data
square = geometry.Polygon([(0,0),(1,0),(1,1),(0,1)])
triangle = geometry.Polygon([(0.5,0.5),(1.5,0.5),(1,1.5)])
square.intersects(triangle) # Returns true if overlap
isect = square.intersection(triangle) # 'overlap'
bsect = square.intersection(triangle.boundary)
diff = square.difference(triangle) # 'subtracts'
uni = triangle.union(square) # 'combines'

# Plotting using Descartes
triPatch = PolygonPatch(triangle, fc='b', ec='k', alpha=0.5, zorder=2)
sqPatch = PolygonPatch(square, fc='r', ec='k', alpha=0.5, zorder=2)
iPatch = PolygonPatch(isect, fc='g', ec='k', alpha=0.9, zorder=2,hatch='/')
dPatch = PolygonPatch(dsect, fc='none', ec='k', alpha=0.9, zorder=2,hatch='-')

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

# WRITING VECTOR DATA
# WHA: This doesn't work quite right as of 20 spe 2016
# Create schema
if 0: # if 0 turns of section because broken
	mySchema = { 'geometry':'Polygon',
				'properties': {
					'id':int
					}
				}		
	if 1:
		geoms = geometry.MultiPolygon([square,triangle,isect])
		sink = fiona.open('out.shp', 'w',driver=c.driver,crs=crs.from_epsg(4326), schema = mySchema)
	
		for i in range(0,len(geoms)):
			id = i
			print id
			geom = geometry.shape(feat['geometry'])		
			prop = {'id':id}	
			sink.write({'geometry': geometry.mapping(geom), 'properties':prop}) 
			write(feat)
		
		
	if 0:			
		# Open file to write; define driver, crs, and schema
		with fiona.open('out.shp', 'w',driver=c.driver,crs=crs.from_epsg(4326),
					schema = source.schema) as sink:
			# Iterate over features to write
			for feat in geoms:
				# get geometry using shapely           
				geom = geometry.shape(feat['geometry'])
				# map geometry to GeoJSON format
				feat['geometry'] = geometry.mapping(geom)
				sink.write(feat) # write geometry to file
		
		

### RASTER DATA
# Also getting segmentation fault, but not if just import rasterio & matplotlib
src = rasterio.open(rasterFn) # open raster
ll = src.index(472996,4426147) # approx lower left for boulder subset
ur = src.index(484804,4434040) # approx upper right for boulder subset
w = ( (ur[0],ll[0]) , (ll[1],ur[1]) ) # defining window for subset. (row start, row end) (col start, col end)
sub = src.read(1,window=w) # subset image
bLow,bHigh = np.percentile(sub,[1,98]) # calculate percentile of brightnesses
lims = np.percentile(sub,[0,99.5]) # use percentiles for plot limits
hist,edges = np.histogram(sub,bins=500) # make histogram of image brightness
hist = np.insert(hist,0,0) # stick a 0 a the start so len(edges)==len(hist)

# Plot histogram and label
plt.plot(edges,hist.astype('float')/np.sum(hist),lw=2)
plt.plot((bLow,bLow),(0,0.035),ls='--',c='k')
plt.plot((bHigh,bHigh),(0,0.035),ls='--',c='k',label=( ('1%-98%') ))
plt.xlabel('Brightness [-]',fontsize=16)
plt.ylabel('P(Brightness)',fontsize=16)
plt.legend(loc=1,frameon=False)
plt.xlim(lims)
plt.ylim((0,0.035))
plt.savefig('histogram.pdf')
plt.show()
plt.close()

# Plot raster image
xmin,ymax = src.affine * (0,0) # get geographic coordinates of raster
xmax,ymin = src.affine * (src.height,src.width) # get geographic coordinates of raster
plt.imshow(sub,cmap='gray',origin='upper',clim=((bLow,bHigh)),extent=((xmin,xmax,ymin,ymax)),aspect=1)
plt.xlabel('Easting [m]',fontsize=16)
plt.ylabel('Northing [m]',fontsize=16)
plt.savefig('boulder.png',r=600)
plt.show()
plt.close()

# WRITING RASTER DATA
from skimage import exposure
stretch=np.copy(sub)
stretch = exposure.rescale_intensity(stretch,(bLow,bHigh))
new_dataset = rasterio.open('boulderSub.tif', 'w', driver='GTiff',
                            height=sub.shape[0], width=sub.shape[1],
                            count=1, dtype=sub.dtype,
                            crs=src.crs, transform=src.transform)
new_dataset.write(stretch,1)                            
new_dataset.close()
