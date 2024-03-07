import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import geopandas as gpd
from shapely.geometry import Point, LineString, MultiLineString, box
import rasterio
import networkx as nx
import numpy as np
from utils import convert_point_crs, haversine_distance, filter_lines_in_bbox, is_point_crossing_segment

class RoadNet:
    def __init__(self, road_data, road_img):
        # Read the road data from GeoJSON using geopandas
        roads = gpd.read_file(road_data)
        # in KTH campus, our case only consider unclassified type and residential type (check on OpenStreetMap)
        roads = roads[roads['highway'].isin(['unclassified', 'residential'])]
        # Read the .tif file
        self.image = rasterio.open(road_img)
        self.roads = roads[(roads.geometry.type == 'MultiLineString')] # only care about lines
        
        
        
    def trim_road(self, min_x=18.0655, min_y=59.3475, max_x=18.0715, max_y=59.354):
        """
        Trim the road network based on the defined ROI.
        """
        # roi coordinates
        self.min_x = min_x
        self.min_y = min_y
        self.max_x = max_x
        self.max_y = max_y
        
        # Create a bounding box polygon
        bbox_polygon = box(self.min_x, self.min_y, self.max_x, self.max_y)
        # create roads_in_bbox as the trimmed network within bbox
        roads_in_bbox = self.roads.copy()
        roads_in_bbox['geometry'] = self.roads['geometry'].apply(filter_lines_in_bbox, 
                                                                 bbox_polygon=bbox_polygon)
        # Clip road segments to the boundary of the ROI
        self.clipped_roads_GEO = roads_in_bbox
        # Convert  EPSG:4326 to EPSG:3854, lag/lat --> meter      
        self.clipped_roads_CRS = self.clipped_roads_GEO.to_crs(epsg=3854) # epsg=3854 for sweden https://epsg.io/3854
        
        
        
    def draw_road_fig(self):
        """
        Draw the KTH geo-image and motor roads
        
        """
        # Read the region of interest from the raster image
        window = rasterio.windows.from_bounds(self.min_x, self.min_y, self.max_x, self.max_y, self.image.transform)
        roi = self.image.read(window=window)

        # Plot the ROI
        plt.figure(figsize=(5, 10))
        plt.axis('off')
        plt.imshow(roi.transpose(1, 2, 0), extent=(self.min_x, self.max_x, self.min_y, self.max_y))

        # Plot the road data
        self.clipped_roads_GEO.plot(ax=plt.gca(), linewidth=1.5, color='red')
        
        plt.title("KTH geo-image and motor roads")
        plt.show()
    
    
    
    def draw_network_fig(self, original_points=None, shortest_path=None):
        """
        Draw the KTH graph network with nodes and edges, optionally highlighting the shortest path between original_points.
    
        Parameters:
            original_points (list of tuples): List of original points to be drawn. 
                            Each point should be a tuple of (x, y) coordinates. Default is None.
            shortest_path (list): List of nodes representing the shortest path. Default is None.

        """
        positions = {n: [n[0], n[1]] for n in self.G.nodes}
        
        fig, ax = plt.subplots(figsize=(5, 10))
        
        # Draw the network without highlighting the shortest path
        nx.draw_networkx(self.G, positions, with_labels=False, node_size=20, ax=ax)
        
        # # Get edge labels (weights)
        # edge_labels = nx.get_edge_attributes(self.G, 'weight')
        # formatted_edge_labels = {edge: f'{label:.3f}' for edge, label in edge_labels.items()}
        # # Draw edge labels
        # nx.draw_networkx_edge_labels(self.G, positions, edge_labels=formatted_edge_labels, ax=ax)
        
        # Highlight the shortest path, if provided
        if shortest_path:
            shortest_path_edges = [(shortest_path[i], shortest_path[i+1]) for i in range(len(shortest_path)-1)]
            nx.draw_networkx_edges(self.G, positions, edgelist=shortest_path_edges, edge_color='g', width=3, ax=ax)
    
        # Draw the original points, if provided
        if original_points:
            x, y = zip(*original_points)
            ax.scatter(x, y, color='r', marker='x', label='Original Points')
        
        # Set x-axis formatter to display in 10e3 format
        ax.xaxis.set_visible(False)
        # Add x and y axes
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
        plt.title("KTH graph network with nodes and edges")
        plt.ylabel('in CRS units [m]')
        plt.show()
    
    
    def _build_network_graph(self):
        """
        Create a network graph from road segments within the ROI.
        If road is a MultiLineString, iterate over its LineString components and add them as edges to the graph.
        
        """
        G = nx.Graph()
        
        for index, road in self.clipped_roads_CRS.iterrows():
            if isinstance(road.geometry, MultiLineString):
                for line in road.geometry.geoms:  # Iterate over the individual LineStrings
                    for i in range(len(line.coords) - 1):
                        coord1, coord2 = line.coords[i], line.coords[i + 1]
                        G.add_edge(coord1, coord2, weight=LineString([coord1, coord2]).length)
            else:
                raise ValueError("Unexpected geometry type")
        self.G = G
        
        
            
    def shortest_distance_along_roads(self, point1, point2):
        """
        Given two points, calculate the shortest path and the length of the shortest path
        
        """
        if not hasattr(self, 'G'):
            self._build_network_graph()
        
        # find the nearest touching points (nodes or points on edges) on the road network to the given points
        distance1, touching_point1 = self.find_nearest_distance_and_point_to_graph(point1)
        distance2, touching_point2 = self.find_nearest_distance_and_point_to_graph(point2)
        
        # Find the nearest nodes on the road network to the given points
        nearest_node1_idx = np.argmin([touching_point1.distance(Point(node)) for node in self.G.nodes])
        nearest_node2_idx = np.argmin([touching_point2.distance(Point(node)) for node in self.G.nodes])
        nearest_dis1 = np.min([touching_point1.distance(Point(node)) for node in self.G.nodes])
        nearest_dis2 = np.min([touching_point2.distance(Point(node)) for node in self.G.nodes])
        
        # Retrieve the corresponding node identifiers
        nearest_node1 = list(self.G.nodes())[nearest_node1_idx]
        nearest_node2 = list(self.G.nodes())[nearest_node2_idx]
        # Calculate the shortest path, node-to-node
        shortest_path = nx.shortest_path(self.G, source=nearest_node1, target=nearest_node2, weight='weight')
        shortest_distance = nx.shortest_path_length(self.G, source=nearest_node1, target=nearest_node2, weight='weight')
        
        # compensate the distance between nodes and touching points
        if self.is_path_crossing_point(shortest_path, (touching_point1.x, touching_point1.y)):
            shortest_distance = shortest_distance - nearest_dis1
        else:
            shortest_distance = shortest_distance + nearest_dis1
        
        if self.is_path_crossing_point(shortest_path, (touching_point2.x, touching_point2.y)):
            shortest_distance = shortest_distance - nearest_dis2
        else:
            shortest_distance = shortest_distance + nearest_dis2
        
        # draw the shortest path and points on the road network
        touch_points = [(touching_point1.x, touching_point1.y), (touching_point2.x, touching_point2.y)]
        self.draw_network_fig(original_points=touch_points, shortest_path=shortest_path)

        return shortest_distance, shortest_path

        
        
    def find_nearest_distance_and_point_to_graph(self, point):
        """
        Find the closest points in the road network to the given point

        """
        if not hasattr(self, 'G'):
            self._build_network_graph()
        
        
        point_crs = convert_point_crs(point)
        
        # Initialize variables for storing nearest distance and touching point
        nearest_distance = np.inf
        touching_point = None
        
        # Iterate through nodes and edges in the graph
        for node in self.G.nodes:
            distance_to_node = point_crs.distance(Point(node))[0]
            if distance_to_node < nearest_distance:
                nearest_distance = distance_to_node
                touching_point = node
        
        for u, v, data in self.G.edges(data=True):
            edge = LineString([Point(u), Point(v)])
            distance_to_edge = point_crs.distance(edge)[0]
            if distance_to_edge < nearest_distance:
                nearest_distance = distance_to_edge
                touching_point = edge.interpolate(edge.project(point_crs))
        
        # if touching_point is a GeoDataFrame, convert it to a point
        if not isinstance(touching_point, tuple):
            touching_point = Point(touching_point['geometry'].iloc[0].x, touching_point['geometry'].iloc[0].y)
        
        return nearest_distance, touching_point
    
    
    
    def is_path_crossing_point(self, path, point):
        """
        Check if the point is across a line segment defined by the path

        """
        for u_coords, v_coords in zip(path[:-1], path[1:]):            
            if is_point_crossing_segment(point, u_coords, v_coords):
                return True
        return False
    
    
    