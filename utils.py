import geopandas as gpd
from math import radians, sin, cos, sqrt, atan2
from shapely.geometry import Point, MultiLineString, LineString

def convert_point_crs(point, target_crs='EPSG:3854'):
    """
    Convert a Shapely Point to a specified Coordinate Reference System (CRS).
    
    Parameters:
        point (shapely.geometry.Point or tuple): The input point geometry as a Point object or a tuple (lon, lat).
        target_crs (str or dict): The target CRS to which the point will be converted.
            Can be specified as an EPSG code (e.g., 'EPSG:4326') or as a dictionary
            representing the CRS definition.
    
    Returns:
        geopandas.GeoDataFrame: A GeoDataFrame containing the point geometry in the target CRS.
    """
    # Convert tuple to Point if necessary
    if isinstance(point, tuple):
        point = Point(point)

    # Create a GeoDataFrame with the point
    gdf_point = gpd.GeoDataFrame(geometry=[point], crs='EPSG:4326')  # Assume WGS84 (EPSG:4326) for initial CRS
    
    # Convert the CRS
    gdf_point.to_crs(crs=target_crs, inplace=True)
    
    return gdf_point



def haversine_distance(point1, point2):
    """
    Calculate the distance between two points on the Earth's surface
    using the Haversine formula.

    Parameters:
        point1 (tuple): Longitude and latitude of the first point in degrees.
        point2 (tuple): Longitude and latitude of the second point in degrees.

    Returns:
        distance (float): Distance between the two points in meters.
    """
    # Convert latitude and longitude from degrees to radians
    lon1, lat1 = radians(point1[0]), radians(point1[1])
    lon2, lat2 = radians(point2[0]), radians(point2[1])

    # Radius of the Earth in meters
    R = 6371000  # Earth radius in meters

    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    distance = R * c

    return distance



def filter_lines_in_bbox(multilinestring, bbox_polygon):
    """
    Function to filter out line segments not in the bounding box
    """
    filtered_lines = []
    for line in multilinestring.geoms:  # Iterate over the individual LineStrings in the MultiLineString
        if line.intersects(bbox_polygon):
            intersection = line.intersection(bbox_polygon)
            if intersection.geom_type == 'LineString':
                filtered_lines.append(intersection)
            elif intersection.geom_type == 'MultiLineString':
                filtered_lines.extend(intersection.geoms)  # Iterate over the individual LineStrings in the MultiLineString
    return MultiLineString(filtered_lines)


