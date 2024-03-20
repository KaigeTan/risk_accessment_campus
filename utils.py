import geopandas as gpd
from shapely.geometry import Point, MultiLineString
import math
from pyproj import Proj, transform

def convert_point_crs(point, target_crs='EPSG:3854'):
    """
    Convert a Shapely Point to a specified Coordinate Reference System (CRS).
    
    Parameters:
        point (shapely.geometry.Point or tuple): The input point geometry as a Point object or a tuple (lon, lat).
        target_crs (str or dict): The target CRS to which the point will be converted.
            Can be specified as an EPSG code (e.g., 'EPSG:4326') or as a dictionary
            representing the CRS definition.
        
        *** Note: EPSG:3854 is for Stockholm/Sweden, check on https://epsg.io/ ***
    
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



def filter_lines_in_bbox(multilinestring, bbox_polygon):
    """
    Function to filter out line segments not in the bounding box
    
    Parameters:
        multilinestring (shapely.geometry.MultiLineString): The MultiLineString containing line segments.
        bbox_polygon (shapely.geometry.Polygon): The bounding box polygon defining the area of interest.

    Returns:
        shapely.geometry.MultiLineString: A MultiLineString containing line segments within the bounding box.
        
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



def is_point_crossing_segment(point, segment_start, segment_end):
    """
    Check if the point is crossing a line segment defined by segment_start and segment_end.
    
    Parameters:
        point (tuple): The coordinates of the point as a tuple (x, y).
        segment_start (tuple): The coordinates of the start point of the line segment as a tuple (x, y).
        segment_end (tuple): The coordinates of the end point of the line segment as a tuple (x, y).
    
    Returns:
        bool: True if the point crosses the line segment, o.w. False.
        
    """
    # Calculate vectors representing the line segment and the point to segment_start
    vec_segment = (segment_end[0] - segment_start[0], segment_end[1] - segment_start[1])
    vec_point_to_start = (segment_start[0] - point[0], segment_start[1] - point[1])

    # Calculate the cross product
    cross_product = vec_segment[0] * vec_point_to_start[1] - vec_segment[1] * vec_point_to_start[0]

    # If the cross product is zero and the dot product is negative, the point is on the line segment
    if abs(cross_product) < 0.05:
        dot_product = vec_segment[0] * vec_point_to_start[0] + vec_segment[1] * vec_point_to_start[1]
        if dot_product < 0:
            return True
    return False



def time_to_collision(distance, vel_a, vel_b, acc_a=0, acc_b=0):
    """
    Function to calculate the time to collision between two objects.
    *** Note: If the vehicle is driving towards the intersection, the velocity is positive, o.w. negative. ***

    Parameters:
        distance (float): Initial distance between the objects.
        vel_a (float): Initial velocity of object A.
        vel_b (float): Initial velocity of object B.
        acc_a (float, optional): Acceleration of object A. Default is 0.
        acc_b (float, optional): Acceleration of object B. Default is 0.
        
    Returns:
        float: Time to collision in seconds.

    """
    relative_velocity = vel_a + vel_b
    relative_acceleration = acc_a + acc_b

    if relative_acceleration == 0:
        # Prevent division by zero
        if relative_velocity <= 0:
            return None  # Objects are not approaching each other
        else:
            time_to_collision = distance/relative_velocity
    else:
        # Calculate time to collision using the quadratic formula
        # Discriminant
        discriminant = relative_velocity**2 + 2*relative_acceleration*distance

        # Check if roots are real
        if discriminant < 0:
            return None  # No real roots, no collision possible
        # Calculate the positive root
        time_to_collision = (-relative_velocity + math.sqrt(discriminant))/relative_acceleration
        # Ensure time to collision is positive
    if time_to_collision < 0:
        time_to_collision = None  # Negative time, no collision possible

    return time_to_collision



# Define a function to map values to colors
def value_to_color(value):
    """
    
    """
    # Normalize the value to range [0, 1]
    normalized_value = (value - 0) / (40 - 0)
    
    # Use the normalized value to interpolate between red and green
    red = 1 - normalized_value
    green = normalized_value
    blue = 0  # No blue component
    
    # Return the color as an RGB tuple
    return (red, green, blue)