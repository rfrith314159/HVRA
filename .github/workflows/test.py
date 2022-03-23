import geopandas
import argparse
import glob
import os
import pandas
import h3
import sys
from shapely.geometry import Polygon
import rasterio
import exposure_idx.shared as share
import pandas as pd
from rasterio.warp import calculate_default_transform, reproject, Resampling
from pathlib import Path
from matplotlib import pyplot as plt
from shapely.geometry import shape
import fiona
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.pandas2ri import py2rpy, rpy2py
import json
import re
import numpy
import ast


polygonise = lambda hex_id: Polygon(
    h3.h3_to_geo_boundary(
        hex_id, geo_json=True)
)

def reproject_raster(in_path, out_path, crs):

    """
    """
    # reproject raster to project crs
    dst_crs = 'EPSG:4326'

    with rasterio.open(in_path) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        with rasterio.open(out_path, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=crs,
                    resampling=Resampling.nearest)

    return(out_path)

def process_raster1(raster):
    rasterFile = os.path.basename(raster) 
    rasterDir = os.path.dirname(raster)
    projectedRaster = rasterDir + "/projected_" + rasterFile
    #print(projectedRaster)
    reproject_raster(raster,projectedRaster,None)
    whp_file = rasterio.open(projectedRaster)
    return whp_file

def build_aoiHexGrid(swsp, aoi, hs): 
    
    stateShp = aoi
    
    shapefile = fiona.open(stateShp)
    
    first = shapefile.next()
    
    hexGridLst = []
    idx = 0
    while first:
        hexs = h3.polyfill(first['geometry'], hs, geo_json_conformant = True)
        hex_grid = geopandas.GeoSeries(list(map(polygonise, hexs)), \
                                     index=hexs, \
                                     crs="EPSG:4326" \
                                    )
        hexGridLst.append(hex_grid)
        idx = idx + 1
        try:
            first = shapefile.next()
        except:
            break

    gdf0 = geopandas.GeoDataFrame(hexGridLst[0])

    gdf0.columns = ['geometry']
    
    for gs in hexGridLst:
        g = geopandas.GeoDataFrame(gs)
        g.columns = ['geometry']
        gdf0 = gdf0.append(g)
    
    
    gdf0 = gdf0.to_crs(epsg=4326)
    #gdf0.to_file(swsp + "//" + 'gdf0.shp' )
    return gdf0

def build_aoiHexGrid1(sourceFolder, aoi, hs): 
    
    stateShp = sourceFolder + "/" + aoi
    
    shapefile = fiona.open(stateShp)
    
    first = shapefile.next()
    
    hexGridLst = []
    idx = 0
    while first:
        hexs = h3.polyfill(first['geometry'], hs, geo_json_conformant = True)
        hex_grid = geopandas.GeoSeries(list(map(polygonise, hexs)), \
                                     index=hexs, \
                                     crs="EPSG:4326" \
                                    )
        hexGridLst.append(hex_grid)
        idx = idx + 1
        try:
            first = shapefile.next()
        except:
            break

    gdf0 = geopandas.GeoDataFrame(hexGridLst[0])

    gdf0.columns = ['geometry']
    
    for gs in hexGridLst:
        g = geopandas.GeoDataFrame(gs)
        g.columns = ['geometry']
        gdf0 = gdf0.append(g)
    
    
    gdf0 = gdf0.to_crs(epsg=4326)
    return gdf0

def build_hexGrid(hs, buffDict,swsp):
   
    gdf = geopandas.GeoDataFrame(pandas.concat(buffDict[key] for key in buffDict.keys())) 
    gdf.columns = ['geometry']        
    gdf['dissolveColumn'] = 1
    dissolved = gdf.dissolve('dissolveColumn')
    #dissolved.to_file(swsp + "//" + 'compiledBuffer.shp' )
    #hexagonal tesselation
    gdf = dissolved   
    gdf.crs = "EPSG:3857"
    gdf = gdf.to_crs(epsg=4326)
    exploded = gdf.explode(ignore_index=True)
    area = 0
    idx = 0
    for x in range(len(gdf) + 1):
        try:
            if(exploded.geometry[x].area >= area):
                area = exploded.geometry[x].area
                idx = x
        except:
            pass
    
    hexs = h3.polyfill(exploded.geometry[idx].__geo_interface__, hs, geo_json_conformant = True)
    
    hex_grid = geopandas.GeoSeries(list(map(polygonise, hexs)), \
                                      index=hexs, \
                                      crs="EPSG:4326" \
                                     )
    
    return hex_grid, dissolved

def build_buffers_from_file(inputs):
   
    hs = 0
    sourceRaster = ''
    aoi = ''
    scratchWorkspace = ''
    unbufferedDict = {}
    unbufferedSubtypeDict = {}
    buffDict = {}
    
    for item in inputs['configuration']:
        #print(item['configtype'])
        if(item['configtype'] == 'parameter' and item['designation'] == 'scratch workspace'):
            scratchWorkspace = item['source']
            print('scratch workspace : ' + scratchWorkspace)
        if(item['configtype'] == 'parameter' and item['designation'] == 'hex scale'):
            try:
                hs = int(item['scale'])
                if hs < 0:  # if not a positive int print message and ask for input again
                    print("Sorry, hex scale be a positive integer <= 15, try again")
                    return
                else:
                    print('hex scale : ' + str(hs))    
            except ValueError:
                print("Input a valid numeric hexagonal scale.")
                return
    
    for item in inputs['systems']:
        if(item['designation'] == 'infrastructure'):
            #print(item['source'] + ' ' + str(item['buffer_distance']))
            if(item['subtypes'] == '{}'):
                unbufferedDict[item['source']] = item['buffer_distance']
            else:
                unbufferedSubtypeDict[item['source']] = item['subtypes']
        if(item['geotype'] == 'polygon' and item['designation'] == 'aoi'):
            aoi = item['source']
            print('area of interest : ' + aoi)
        if(item['geotype'] == 'raster'):
            sourceRaster = item['source']
            print('source raster : ' + sourceRaster)
       
    if(len(unbufferedDict) > 0):
        for key, value in unbufferedDict.items():
            shp = key
            bd = value
            print('key : ' + key + ', buffer distance : ' + str(value))
            file = os.path.basename(shp)
            if "Buffer" not in file: 
                dfPoly = geopandas.read_file(shp)
                buff = geopandas.GeoSeries(dfPoly["geometry"])
                buffName = file.split('.')[0] + '_Buffer.shp'           
                buffDict[buffName] = buff.buffer(bd)

    if(len(unbufferedSubtypeDict) > 0):
        for key,value in unbufferedSubtypeDict.items():
            shp = key
            file = os.path.basename(shp)
            print(shp)
            if "Buffer" not in file:
                dfPoly = geopandas.read_file(shp)
                #dfPoly = geopandas.read_file(shp)
                #print(dfPoly.head(10))
                jsonObj = json.loads(value)
                print(type(jsonObj))
                for k,v in jsonObj.items():
                    print('--' + k)
                    v = v.replace('\'','\"')
                    res = json.loads(v)
                    print(type(res))
                    for k1,v1 in res.items():
                        print('key : ' + k1 + ' value : ' + str(v1))
                        try:
                            
                            dfPoly[dfPoly[k] == k1].to_file(scratchWorkspace + "//" + file.split('.')[0] + '_' + k1 + ".shp")
                            tmpPoly = dfPoly[dfPoly[k] == k1]
                            buff = geopandas.GeoSeries(tmpPoly["geometry"])
                            buffName = file.split('.')[0] + '_' + k1 + '_Buffer.shp'
                            buffDict[buffName] = buff.buffer(v1)
                        except Exception as ex:
                            print(ex)
                        
   
   
    hex_grid, compiledBuffer = build_hexGrid(hs,buffDict,scratchWorkspace)
    return scratchWorkspace, hs, sourceRaster, aoi, buffDict, hex_grid, compiledBuffer

def area_intersection(polyShp,dfHexPoly,title,refField):
    dfHexPolyTemp = dfHexPoly
    dfPoly = geopandas.read_file(polyShp)
    dfPoly = dfPoly.to_crs(epsg=3857)
    dfHexPolyTemp = dfHexPolyTemp.to_crs(epsg=3857)
    dfHexPolyTemp['doesIntersect'] = 0
    dfHexPolyTemp['area'] = dfHexPolyTemp['geometry'].area/10**6  #square km
    dfHexPolyTemp['index'] = dfHexPolyTemp.index
    dfHexPolyTemp[refField + 'area'] = 0
    
   
    for index, row in dfHexPolyTemp.iterrows():
        polyClip = geopandas.clip(dfPoly,row['geometry'])
        polyClip = polyClip[~polyClip.is_empty]
        lst = polyClip.area/10**6
        s = lst.to_frame()
        
        s.columns = ['geometry']
        if not s.empty:
            dfHexPolyTemp.loc[index, refField + 'area'] = s['geometry'].sum()
    
            dfHexPolyTemp.loc[index, 'doesIntersect'] = 1
    
    dfHexPolyTemp['index'] = dfHexPolyTemp.index  
    return dfHexPolyTemp
    
def area_intersection_compiledBuffer(dfHexPoly,compiledBuffer):
    
    dfHexPolyTemp = dfHexPoly
    dfPoly = compiledBuffer
    dfPoly = dfPoly.to_crs(epsg=3857)
    dfHexPolyTemp = dfHexPolyTemp.to_crs(epsg=3857)
    dfHexPolyTemp['polyArea'] = 0
    dfHexPolyTemp['doesIntersect'] = 0
    dfHexPolyTemp['area'] = dfHexPolyTemp['geometry'].area/10**6  #square km
    
    dfHexPolyTemp['catArea'] = ''
    
   
    for index, row in dfHexPolyTemp.iterrows():
        polyClip = geopandas.clip(dfPoly,row['geometry'])
        polyClip = polyClip[~polyClip.is_empty]
        lst = polyClip.area/10**6
        s = lst.to_frame()
        #dfHexPolyTemp.loc[index, 'catArea'] = '{ "fromArea" :'
        s.columns = ['geometry']
        if not s.empty:
            dfHexPolyTemp.loc[index, 'polyArea'] = s['geometry'].sum()
            dfHexPolyTemp.loc[index, 'doesIntersect'] = 1
          

    dfHexPolyTemp['index'] = dfHexPolyTemp.index  
    return dfHexPolyTemp
    
#############################
#https://gis.stackexchange.com/questions/332167/length-of-intersections-from-a-linestring-and-a-grid-shapefile-by-using-python-g
#https://gis.stackexchange.com/questions/340058/geopandas-cut-lines-with-polygon

    
# plt.ion()
# fig, ax = plt.subplots()
# dfHexPoly.plot(ax=ax, alpha=0.1)    
# plt.savefig(sourceFolder + "//" + 'testfig.png' )
def line_intersection1(lineShp,dfHexPoly,title,refField):
    
    dfHexPolyTemp = dfHexPoly
    dfLine = geopandas.read_file(lineShp)   
    dfLine = dfLine.to_crs(epsg=3857)
    dfHexPolyTemp = dfHexPolyTemp.to_crs(epsg=3857)
    dfHexPolyTemp['index'] = dfHexPolyTemp.index
    
   
    dfHexPolyTemp[title + 'Length'] = 0
    for index, row in dfHexPolyTemp.iterrows():
        lineClip = geopandas.clip(dfLine,row['geometry'])
        lineClip = lineClip[~lineClip.is_empty]
        lineSet = set(lineClip[refField].to_list())
        if(len(lineSet) == 0): 
            continue

        for line in lineSet:
            dfLin = lineClip[lineClip[refField] == line]
           
            lst = dfLin.length
            #print(lst)
            s = lst.to_frame()
            s.columns = ['geometry']
           
            dfHexPolyTemp.loc[index,title + 'Length'] += s['geometry'].sum()
   
    return dfHexPolyTemp

def point_intersection1(ptShp,dfHexPoly,title,refField):
    #run pip uninstall rtree
    #https://automating-gis-processes.github.io/site/notebooks/L3/spatial_index.html
    #https://gis.stackexchange.com/questions/306674/geopandas-spatial-join-and-count
    #dfHexPtPoly['type'] = 'ptPp'
    dfTmpHexPoly = dfHexPoly
    #dfTmpHexPoly[title + 'Cnt'] = 0
    dfPt = geopandas.read_file(ptShp)   
    dfPt = dfPt.to_crs(epsg=3857)
    dfTmpHexPoly = dfTmpHexPoly.to_crs(epsg = 3857)
    dfTmpHexPoly['index'] = dfTmpHexPoly.index

    ptLst = geopandas.sjoin(dfTmpHexPoly, dfPt, how='left').groupby('index')[refField].apply(set).to_frame()
    #ptLst.to_csv('C:/Users/russ.frith/P1/Data/buffers/check.csv')
    ptLst['ptIndex'] = ptLst.index
   
    ptLst[title + 'Cnt'] = 0 
    ptLst[title + 'Cnt'] = ptLst[refField].apply(lambda x: len(x) if str(x) != '{nan}' else 0)
    
   
    dfTmpHexPoly['ptIndex'] = dfTmpHexPoly['index']
    dfPtTotal = pandas.merge(dfTmpHexPoly, ptLst, how='inner', on=['ptIndex'])
    
    df = dfPtTotal[["geometry","index",title + 'Cnt']]
    return df

def combineHexPtPoly(scratch,dfHexPtPolyLst,ptCategories):
    dfHexPtPoly = dfHexPtPolyLst.pop()
    
    
    for ptPoly in dfHexPtPolyLst:
        try:
            dfHexPtPoly = dfHexPtPoly.merge(ptPoly,on='index')           
        except Exception as ex:
            print('cannot merge : ',ex)
    print(type(dfHexPtPoly))
    print(dfHexPtPoly.head(10))
    delColLst = []
    for col in dfHexPtPoly.columns:
        if(col.endswith("_y")):
            delColLst.append(col)

    
    dfHexPtPoly.drop(delColLst, axis=1, inplace=True)
    dfHexPtPoly.columns = dfHexPtPoly.columns.str.replace('_x','')
    dfHexPtPoly = dfHexPtPoly.loc[:,~dfHexPtPoly.columns.duplicated()]
    print(type(dfHexPtPoly))
    print(dfHexPtPoly.head(10))

   
    gdfPt = geopandas.GeoDataFrame(dfHexPtPoly)
    gdfPt.to_file(scratch + "//" + 'dfHexPtPoly.shp' )
    
    return gdfPt

def combineHexLinePoly(scratch, dfHexLinePolyLst,lineCategories):
    dfHexLinePoly = dfHexLinePolyLst.pop()
    for linePoly in dfHexLinePolyLst:
        dfHexLinePoly = dfHexLinePoly.merge(linePoly,on='index')

    delColLst = []
    for col in dfHexLinePoly.columns:
        if(col.endswith("_y")):
            delColLst.append(col)

    dfHexLinePoly.drop(delColLst, axis=1, inplace=True)
    dfHexLinePoly.columns = dfHexLinePoly.columns.str.replace('_x','')
    dfHexLinePoly = dfHexLinePoly.loc[:,~dfHexLinePoly.columns.duplicated()]
   
    #dfHexLinePoly.to_csv(scratch + "//" + 'dfHexLinePoly.csv')
    gdfLine = geopandas.GeoDataFrame(dfHexLinePoly)
    gdfLine.to_file(scratch + "//" + 'dfHexLinePoly.shp' )
    return dfHexLinePoly

def combineHexAreaPoly(scratch, dfHexAreaPolyLst,areaCategories):
    dfHexAreaPoly = dfHexAreaPolyLst.pop()
    for areaPoly in dfHexAreaPolyLst:
        dfHexAreaPoly = dfHexAreaPoly.merge(areaPoly,on='index')

    delColLst = []
    for col in dfHexAreaPoly.columns:
        if(col.endswith("_y")):
            delColLst.append(col)

    dfHexAreaPoly.drop(delColLst, axis=1, inplace=True)
    dfHexAreaPoly.columns = dfHexAreaPoly.columns.str.replace('_x','')
    dfHexAreaPoly = dfHexAreaPoly.loc[:,~dfHexAreaPoly.columns.duplicated()]
    #dfHexAreaPoly = merge_area_categories(dfHexAreaPoly,areaCategories)
    #dfHexAreaPoly.to_csv(scratch + "//" + 'dfHexAreaPoly.csv')
    gdfArea = geopandas.GeoDataFrame(dfHexAreaPoly)
    gdfArea.to_file(scratch + "//" + 'dfHexAreaPoly.shp' )
    return dfHexAreaPoly

def combineHexPoly(scratch,dfHexPtPoly,dfHexLinePoly,dfHexAreaPoly,whipPoly):
    whipPoly['index'] = whipPoly.index
    print(type(dfHexPtPoly))
    dfHexPoly = dfHexPtPoly
    dfHexPoly = dfHexPoly.merge(dfHexLinePoly,on='index')
    print(dfHexPoly.head(10))
    print(type(dfHexPoly))
    dfHexPoly = dfHexPoly.merge(dfHexAreaPoly,on='index')
    dfHexPoly =  dfHexPoly.merge(whipPoly,on='index')

    delColLst = []
    for col in dfHexPoly.columns:
        if(col.endswith("_y")):
            delColLst.append(col)
        
    
    dfHexPoly.drop(delColLst, axis=1, inplace=True)
    dfHexPoly.columns = dfHexPoly.columns.str.replace('_x','')
    dfHexPoly = dfHexPoly.loc[:,~dfHexPoly.columns.duplicated()]
    print(dfHexPoly.columns)
    # #dfHexPoly = merge_categories(dfHexPoly)
    # dfHexPoly.to_csv(scratch + "//" + 'dfHexPoly.csv')
    gdf = geopandas.GeoDataFrame(dfHexPoly)
    gdf.to_file(scratch + "//" + 'dfHexPoly.shp' )
    # # normal_df = r_processing(scratch + "//" + 'dfHexPoly.shp')
    normal_df = normal_processing(scratch + "//" + 'dfHexPoly.shp')
    print(normal_df.head(10))
    for col in normal_df.columns:
        if col.startswith('norm_'):
            gdf[col] = normal_df[col]
        if(col.startswith('tot_')):
            gdf[col] = normal_df[col]
        if(col.startswith('final_')):
            gdf[col] = normal_df[col]
    
    gdf.to_file(scratch + "//" + "normalHexPoly.shp")

    return gdf
    
def normalize(feature):
  normal_feature = (feature - feature.min()) / (feature.max() - feature.min())
  return normal_feature

def rank2weight(x):
  mx = max(x)
  rev = [(mx + 1) - i for i in x]
  weights = [i / sum(rev) for i in rev]
  return weights

def weighted_score(df, features, weights):
  length = len(features)
  tmp_df = pd.DataFrame()
  for i in range(length):
    tmp_var_name = 'weighted_feat_' + str(i + 1)
    w_vals = df[features[i]] * weights[i]
    tmp_df[tmp_var_name] = w_vals
    
  df['wght_score'] = tmp_df.sum(axis=1)
  
  print(df.head(10))
  return df
  
def normal_processing(shapeFile):
    #features = ['schlNameCn', 'interstate', 'polyArea']
    features = ['schlNameCn', 'I1Length', 'p2area']
    data = geopandas.read_file(shapeFile)
    for f in features:
        norm_feat_name = 'norm_' + f
        
        norm_feat = normalize(data[f])
        data[norm_feat_name] = norm_feat
        print(data.head(10))
        print('----')

    normal_features = ['norm_schlNameCn', 'norm_I1Length', 'norm_p2area']
    test_weights = rank2weight([1,2,3])

    return weighted_score(data, normal_features, test_weights)

def r_processing(shapeFile):
    r = robjects.r
    
    print('sourcing R script')
    r['source']('HVRA_scoring.R')
    # Loading the function we have defined in R.
    print('Loading R defined function')
    normalize_function_r = robjects.globalenv['normalizeData']
    
    
    #Invoking the R function and getting the result
    print('getting result from R function')
    try:
        df_result_r = normalize_function_r(shapeFile)
        
        itm = df_result_r[df_result_r.colnames.index('index')]
        list_of_tuples = list(zip(itm))
        df = pd.DataFrame(list_of_tuples,columns=['index'])   

        for col in list(tuple(df_result_r.names)):
            if('index' not in col):
                itm = df_result_r[df_result_r.colnames.index(col)]
                list_of_tuples = list(zip(itm))
                dftmp = pd.DataFrame(list_of_tuples,columns=[col])
                dftmp = dftmp.reset_index(drop=True)
                df = pd.concat([df,dftmp],axis=1)
                
            
        #print(df.head(10))
        return df
        
    except Exception as ex:
        print(ex)
        return None
    
def main(inputs):
    
    
    scratch, hs, sourceRaster, aoi, buffDict, hex_grid, compiledBuffer = build_buffers_from_file(inputs)
    
    compiledBuffer.to_file(scratch + "//" + 'compiledBuffer.shp' )

    dfHexPoly = build_aoiHexGrid(scratch, aoi, hs)
    if dfHexPoly.empty:
        print('empty dfHexPoly')
        return -1
    # Calc median WHP    
    dfHexPoly[share.whp_score] = dfHexPoly.apply(
      share.raster_stats_median, intersection_data=dfHexPoly, raster_file=process_raster1(sourceRaster), normalize=False, axis=1)
    dfHexPoly.to_file(scratch + "//" + 'dfHexPoly.shp' )

    ptPolys = []
    ptCategories = []
    linePolys = []
    areaPolys = []
    lineCategories = []
    areaCategories = []

    ptCatIdx_people_and_property = {}
    lineCatIdx_people_and_property = {}
    polyCatIdx_people_and_property = {}
    ptCatIdx_infrastructure = {}
    lineCatIdx_infrastructure = {}
    polyCatIdx_infrastructure = {}

    catDict = {}

    for item in inputs['systems']:
        if(item['geotype'] == 'point' and item['designation'] == 'infrastructure'):
            #print('item catIdx = ',item['catIdx'])
            ptCatIdx_infrastructure[item['catIdx']] = item['name'] + 'Cnt'
            #print('item[catIdx] : ',item['catIdx'],' ptCatIdx[item[catIdx]] = : ',ptCatIdx[item['catIdx']])
            dfHexPtPoly = point_intersection1(item['source'],dfHexPoly,item['name'],item['referenceField'])
            dfHexPtPoly = dfHexPtPoly.drop_duplicates(keep='first', ignore_index=True)
            key = item['name'] + 'Cnt'
            catDict[key] = 'ptInfra'
            dfHexPtPoly[item['name'] + 'Type'] = 'ptInfra'
            ptPolys.append(dfHexPtPoly)
            #ptCategories.append('cat' + item['referenceField'])
            ptCategories.append(item['name'] + 'Cnt')
            dfHexPtPoly.to_file(scratch + "//" + 'dfHexPt_' + item['referenceField'] + '_Poly.shp' )
        if(item['geotype'] == 'point' and item['designation'] == 'people_and_property'):
            #print('item catIdx = ',item['catIdx'])
            ptCatIdx_people_and_property[item['catIdx']] = item['name'] + 'Cnt'
            #print('item[catIdx] : ',item['catIdx'],' ptCatIdx[item[catIdx]] = : ',ptCatIdx[item['catIdx']])            
            dfHexPtPoly = point_intersection1(item['source'],dfHexPoly,item['name'],item['referenceField'])
            dfHexPtPoly = dfHexPtPoly.drop_duplicates(keep='first', ignore_index=True)
            dfHexPtPoly[item['name'] + 'Type'] = 'ptPp'
            key = item['name'] + 'Cnt'
            catDict[key] = 'ptPp'
            ptPolys.append(dfHexPtPoly)
            ptCategories.append(item['name'] + 'Cnt')
            dfHexPtPoly.to_file(scratch + "//" + 'dfHexPt_' + item['referenceField'] + '_Poly.shp' )
        if(item['geotype'] == 'line' and item['designation'] == 'infrastructure'):
            lineCatIdx_infrastructure[item['catIdx']] = item['name'] + 'Length'
            dfHexLinePoly = line_intersection1(item['source'],dfHexPoly,item['name'],item['referenceField'])
            dfHexLinePoly = dfHexLinePoly.drop_duplicates(keep='first', ignore_index=True)
            dfHexLinePoly[item['name'] + 'Type'] = 'lineInfra'
            key = item['name'] + 'Length'
            catDict[key] = 'lineInfra'
            linePolys.append(dfHexLinePoly)
            #lineCategories.append('cat' + item['referenceField'])
            dfHexLinePoly.to_file(scratch + "//" + 'dfHexLine_' + item['referenceField'] + '_Poly.shp' )
        if(item['geotype'] == 'line' and item['designation'] == 'people_and_property'):
            lineCatIdx_people_and_property[item['catIdx']] = item['name'] + 'Length'
            dfHexLinePoly = line_intersection1(item['source'],dfHexPoly,item['name'],item['referenceField'])
            dfHexLinePoly = dfHexLinePoly.drop_duplicates(keep='first', ignore_index=True)
            dfHexLinePoly[item['name'] + 'Type'] = 'linePp' 
            linePolys.append(dfHexLinePoly)
            key = item['name'] + 'Length'
            catDict[key] = 'linePp'
            #lineCategories.append('cat' + item['referenceField'])
            dfHexLinePoly.to_file(scratch + "//" + 'dfHexLine_' + item['referenceField'] + '_Poly.shp' )
        if(item['geotype'] == 'polygon' and item['designation'] == 'infrastructure'):
            polyCatIdx_infrastructure[item['catIdx']] = 'cat' + item['referenceField'] 
            dfHexAreaPoly = area_intersection(item['source'],dfHexPoly,item['name'],item['referenceField']) 
            dfHexAreaPoly = dfHexAreaPoly.drop_duplicates(keep='first', ignore_index=True)
            dfHexAreaPoly[item['name'] + 'Type'] = 'polyInfra'
            key = item['name'] + 'area'
            catDict[key] = 'polyInfra'
            areaPolys.append(dfHexAreaPoly)
            #areaCategories.append('cat' + item['referenceField'])           
            dfHexAreaPoly.to_file(scratch + "//" + 'dfHexArea_' + item['referenceField'] + '_Poly.shp' )
        if(item['geotype'] == 'polygon' and item['designation'] == 'people_and_property'):
            polyCatIdx_people_and_property[item['catIdx']] = 'cat' + item['referenceField'] 
            dfHexAreaPoly = area_intersection(item['source'],dfHexPoly,item['name'],item['referenceField']) 
            dfHexAreaPoly = dfHexAreaPoly.drop_duplicates(keep='first', ignore_index=True)
            dfHexAreaPoly[item['name'] + 'Type'] = 'polyPp'
            key = item['name'] + 'area'
            catDict[key] = 'polyPp'
            areaPolys.append(dfHexAreaPoly)
            #areaCategories.append('cat' + item['referenceField'])           
            dfHexAreaPoly.to_file(scratch + "//" + 'dfHexArea_' + item['referenceField'] + '_Poly.shp' )
            
    dfHexPtPoly = combineHexPtPoly(scratch,ptPolys,ptCategories)
    dfHexLinePoly = combineHexLinePoly(scratch,linePolys,lineCategories) 
    dfHexAreaPoly = combineHexAreaPoly(scratch,areaPolys,areaCategories)
    dfHexCombinedAreaPoly = area_intersection_compiledBuffer(dfHexPoly, compiledBuffer)
    dfHexCombinedAreaPoly = dfHexCombinedAreaPoly.drop_duplicates(keep='first', ignore_index=True) 
    #dfHexCombinedAreaPoly.to_file(scratch + "//" + 'dfHexCombined_Poly.shp' )
    df = combineHexPoly(scratch,dfHexPtPoly,dfHexLinePoly,dfHexAreaPoly,dfHexPoly)
    
    
    # print(dfHexPtPoly.columns)
    # print(dfHexLinePoly.columns)
    # print(dfHexAreaPoly.columns)
    
    cat_field = {
        "people_and_property": {
        "point_count": 0,
        "line_dist": 0,
        "poly_area": 0,
        "sub_types": []
        },
        "infrastructure": {
        "point_count": 0,
        "line_dist": 0,
        "poly_area": 0,
        "sub_types": []
        },
        "cultural": {
            "point_count": 0,
            "line_dist": 0,
            "poly_area": 0,
            "sub_types": []             
        },
        "drinking_water": {
            "point_count": 0,
            "line_dist": 0,
            "poly_area": 0,
            "sub_types": []
        }
    }

    cat_field1 = {
        "people_and_property": {
        "point_count": 0,
        "line_dist": 0,
        "poly_area": 0,
        "sub_types": []
        },
        "infrastructure": {
        "point_count": 0,
        "line_dist": 0,
        "poly_area": 0,
        "sub_types": []
        },
        "cultural": {
            "point_count": 0,
            "line_dist": 0,
            "poly_area": 0,
            "sub_types": []
        },
        "drinking_water": {
            "point_count": 0,
            "line_dist": 0,
            "poly_area": 0,
            "sub_types": []
        }
    }

    df['category'] =  json.dumps(cat_field)
    print(df.columns)
    df = populate_catField(df,catDict)
    
    df.to_file(scratch + "//" + "normalHexPoly.shp")
    df.to_file(scratch + "//" + "normalHexPoly.geojson", driver='GeoJSON')
    df.to_csv(scratch + "//" + "normalHexPoly.csv")

def populate_catField(df,catDict):
    for index, row in df.iterrows():
        cfj = ast.literal_eval(df.loc[index,"category"])
        for key,value in catDict.items():
            #print('key = ',key,' value = ',value)
            if(value == 'ptPp'):
                tmpLst = cfj["people_and_property"]["sub_types"]
                tmpDict = {}
                tmpDict['name'] = key
                tmpDict['point_count'] = df.loc[index,key]
                tmpDict['line_dist'] = 0
                tmpDict['poly_area'] = 0
                tmpLst.append(tmpDict)
                cfj["people_and_property"]["sub_types"] = tmpLst
                cfj["people_and_property"]["point_count"] += df.loc[index,key]
                df.loc[index,"category"] = str(cfj)
                #print(df.loc[index,"category"])
            if(value == 'ptInfra'):
                tmpLst = cfj["infrastructure"]["sub_types"]
                tmpDict = {}
                tmpDict['name'] = key
                tmpDict['point_count'] = df.loc[index,key]
                tmpDict['line_dist'] = 0
                tmpDict['poly_area'] = 0
                tmpLst.append(tmpDict)
                cfj["infrastructure"]["sub_types"] = tmpLst
                cfj["infrastructure"]["point_count"] += df.loc[index,key]
                df.loc[index,"category"] = str(cfj)
                #print(df.loc[index,"category"])
            if(value == 'linePp'):
                tmpLst = cfj["people_and_property"]["sub_types"]
                tmpDict = {}
                tmpDict['name'] = key
                tmpDict['point_count'] = 0
                tmpDict['line_dist'] = df.loc[index,key]
                tmpDict['poly_area'] = 0
                tmpLst.append(tmpDict)
                cfj["people_and_property"]["sub_types"] = tmpLst
                cfj["people_and_property"]["line_dist"] += df.loc[index,key]
                df.loc[index,"category"] = str(cfj)
                #print(df.loc[index,"category"])
            if(value == 'lineInfra'):
                tmpLst = cfj["infrastructure"]["sub_types"]
                tmpDict = {}
                tmpDict['name'] = key
                tmpDict['point_count'] = 0
                tmpDict['line_dist'] = df.loc[index,key]
                tmpDict['poly_area'] = 0
                tmpLst.append(tmpDict)
                cfj["infrastructure"]["sub_types"] = tmpLst
                cfj["infrastructure"]["line_dist"] += df.loc[index,key]
                df.loc[index,"category"] = str(cfj)
            if(value == 'linePp'):
                tmpLst = cfj["people_and_property"]["sub_types"]
                tmpDict = {}
                tmpDict['name'] = key
                tmpDict['point_count'] = 0
                tmpDict['line_dist'] = df.loc[index,key]
                tmpDict['poly_area'] = 0
                tmpLst.append(tmpDict)
                cfj["people_and_property"]["sub_types"] = tmpLst
                cfj["people_and_property"]["line_dist"] += df.loc[index,key]
                df.loc[index,"category"] = str(cfj)
                #print(df.loc[index,"category"])
            if(value == 'polyInfra'):
                tmpLst = cfj["infrastructure"]["sub_types"]
                tmpDict = {}
                tmpDict['name'] = key
                tmpDict['point_count'] = 0
                tmpDict['line_dist'] = 0
                tmpDict['poly_area'] = df.loc[index,key]
                tmpLst.append(tmpDict)
                cfj["infrastructure"]["sub_types"] = tmpLst
                cfj["infrastructure"]["poly_area"] += df.loc[index,key]
                df.loc[index,"category"] = str(cfj)
                #print(df.loc[index,"category"])
            if(value == 'polyPp'):
                tmpLst = cfj["people_and_property"]["sub_types"]
                tmpDict = {}
                tmpDict['name'] = key
                tmpDict['point_count'] = 0
                tmpDict['line_dist'] = 0
                tmpDict['poly_area'] = df.loc[index,key]
                tmpLst.append(tmpDict)
                cfj["people_and_property"]["sub_types"] = tmpLst
                cfj["people_and_property"]["poly_area"] += df.loc[index,key]
                df.loc[index,"category"] = str(cfj)
                #print(df.loc[index,"category"])
    return df

# Threat - Color ramp is based on this

# Raw WHP median 
# ((People and Property - Avg Normalized 0.25) + 
# (Infrastructure - Avg Normalized 0.25) + 
# (Cultural - Avg Normalized 0.25) + 
# (Drinking Water - Avg Normalized * 0.25))

if __name__ == "__main__":
    sysDict = {}
    f = open('systems.json')
    inputs = json.load(f)
    f.close()
    main(inputs)


    # Notes:
    # 1. https://github.com/twisted/pydoctor#simple-usage
    # 2. pydoctor --make-html --html-output=docs test.py
    