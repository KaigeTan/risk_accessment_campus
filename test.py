from shapely.geometry import Point
import geopandas as gpd
 
pnt1 = Point(18.069835099071902, 59.34909685637183)
pnt2 = Point(18.072124382629877, 59.34654725898077)
points_df = gpd.GeoDataFrame({'geometry': [pnt1, pnt2]}, crs='EPSG:4326')
points_df = points_df.to_crs('EPSG:3854')
points_df2 = points_df.shift() #We shift the dataframe by 1 to align pnt1 with pnt2
dis = points_df.distance(points_df2)