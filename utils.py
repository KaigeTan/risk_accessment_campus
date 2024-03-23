import geopandas as gpd
import matplotlib.colors as mcolors
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
    if isinstance(point, tuple) or isinstance(point, list):
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



def value_to_color(value):
    """
    Convert a numerical value (TTC) to a corresponding RGB color.

    Parameters:
        value (float): The numerical value to be converted.

    Returns:
        tuple: A tuple representing the RGB color corresponding to the value.
    """
    if value >= 15:
        return (0, 1, 0)  # Bright green
    
    # Define the colors for the gradient (red to yellow to green)
    red = (1, 0, 0)  # Bright red
    yellow = (1, 1, 0)  # Yellow
    green = (0, 1, 0)  # Bright green
    
    # Normalize the value to range [0, 1] for values below 10
    normalized_value = value / 15
    
    # Interpolate between red and yellow (0 to 5)
    if normalized_value <= 0.5:
        # Calculate the interpolation factor
        factor = normalized_value / 0.5
        # Interpolate between red and yellow
        color = (red[0] * (1 - factor) + yellow[0] * factor,
                 red[1] * (1 - factor) + yellow[1] * factor,
                 red[2] * (1 - factor) + yellow[2] * factor)
    # Interpolate between yellow and green (5 to 10)
    else:
        # Calculate the interpolation factor
        factor = (normalized_value - 0.5) / 0.5
        # Interpolate between yellow and green
        color = (yellow[0] * (1 - factor) + green[0] * factor,
                 yellow[1] * (1 - factor) + green[1] * factor,
                 yellow[2] * (1 - factor) + green[2] * factor)
    
    return color



def on_same_line(point1, point2, point3):
    """
    Check if three points are collinear (lie on the same line).

    Args:
    point1, point2, point3: Tuples representing the (x, y) coordinates of the points.

    Returns:
    True if the points are collinear, False otherwise.
    """
    
    # If two points have the same coordinates, they are collinear by default
    if point1 == point2 or point2 == point3 or point1 == point3:
        return True
    
    x1, y1 = point1
    x2, y2 = point2
    x3, y3 = point3

    # Using the formula for slope of a line: (y2 - y1) / (x2 - x1)
    # If the slopes of the lines connecting any two pairs of points are equal, the points are collinear.
    # Avoiding division by zero by checking if the denominators are not equal.
    if (x2 - x1) != 0 and (x3 - x2) != 0:
        slope1 = (y2 - y1) / (x2 - x1)
        slope2 = (y3 - y2) / (x3 - x2)
        return slope1 == slope2
    # If two of the points have the same x-coordinate, they must be collinear with the third point.
    elif (x2 - x1) == 0 and (x3 - x2) == 0:
        return True
    else:
        return False