import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point, LineString, MultiLineString, box
import rasterio
import networkx as nx
import numpy as np
from utils import convert_point_crs, filter_lines_in_bbox, is_point_crossing_segment, time_to_collision, value_to_color

class RoadNet:
    def __init__(self, road_data, road_img=None):
        # Read the road data from GeoJSON using geopandas
        roads = gpd.read_file(road_data)
        # in KTH campus, our case only consider unclassified type and residential type (check on OpenStreetMap)
        roads = roads[roads['highway'].isin(['unclassified', 'residential'])]
        # Read the .tif file
        self.image = rasterio.open(road_img)
        self.roads = roads[(roads.geometry.type == 'MultiLineString')] # only care about lines
        self.image_cnt = 0
        
        # trim the road
        self.trim_road()
        # build the road network
        self._build_network_graph()
        
        
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
        window = rasterio.windows.from_bounds(self.min_x, self.min_y, self.max_x, self.max_y, 
                                              self.image.transform)
        roi = self.image.read(window=window)

        # Plot the ROI
        plt.figure(figsize=(5, 15))
        plt.axis('off')
        plt.imshow(roi.transpose(1, 2, 0), extent=(self.min_x, self.max_x, self.min_y, self.max_y))

        # Plot the road data
        self.clipped_roads_GEO.plot(ax=plt.gca(), linewidth=1, color='red')
        
        plt.title("KTH geo-image and motor roads")
        plt.show()
    
    
    
    def draw_network_fig(self, original_points=None, shortest_path=None, generated_points=None, color=None):
        """
        Draw the KTH graph network with nodes and edges, and also visualize generated points.
        
        Parameters:
            original_points (list of tuples): List of original points to be drawn. 
                            Each point should be a tuple of (x, y) coordinates. Default is None.
            shortest_path (list): List of nodes representing the shortest path. Default is None.
            generated_points (list of tuples): List of generated points to be drawn.
                            Each point should be a tuple of (x, y) coordinates. Default is None.
            color (): Color for the cloest route to show. It indicates the risk value.
        """
        positions = {n: [n[0], n[1]] for n in self.G.nodes}
        
        fig, ax = plt.subplots(figsize=(5, 10))
        
        # Draw the network without highlighting the shortest path
        nx.draw_networkx(self.G, positions, with_labels=False, node_size=20, ax=ax)
        
        # Highlight the shortest path, if provided
        if shortest_path:
            # Define an empty list to store nodes where line segments do not cross original points
            valid_nodes = []
            
            # Iterate over each pair of nodes in the shortest path
            for i in range(len(shortest_path) - 1):
                # Get the current node and the next node in the path
                current_node = shortest_path[i]
                next_node = shortest_path[i + 1]
                edge = [current_node, next_node]
                
                # Draw the edge between the current node and the next node
                if not any(self.is_path_crossing_point(edge, pt) for pt in original_points): # check if line segment cross the points
                    nx.draw_networkx_edges(self.G, positions, edgelist=[(current_node, next_node)], edge_color=color, width=3, ax=ax)
                    valid_nodes.append(current_node)
                # else:
                #     crossing_point = next(pt for pt in original_points if self.is_path_crossing_point(edge, pt))
                #     next_node_coords = positions[next_node]
                #     ax.plot([crossing_point[0], next_node_coords[0]], [crossing_point[1], next_node_coords[1]], color='g', linewidth=3)
            
            # TODO: double check this part, this is a temperal solution right now
            # Find the closest node to the original points among the valid nodes
            closest_node1 = min(valid_nodes, key=lambda node: Point(node).distance(Point(original_points[0])))
            closest_node2 = min(valid_nodes, key=lambda node: Point(node).distance(Point(original_points[1])))
            ax.plot([original_points[0][0], closest_node1[0]], [original_points[0][1], closest_node1[1]], color=color, linewidth=3)
            
            if 100773.2 < original_points[1][0] < 100811.8 and 80980.7 < original_points[1][1] < 81041.6:
                ax.plot([original_points[1][0], 100810.7035522659], [original_points[1][1], 80982.33937331103], color=color, linewidth=3)
            elif 100811.7 < original_points[1][0] and original_points[1][1] < 80980.7:
                ax.plot([original_points[1][0], 100863.51895834938], [original_points[1][1], 80899.63354030065], color=color, linewidth=3)
            else:
                ax.plot([original_points[1][0], closest_node2[0]], [original_points[1][1], closest_node2[1]], color=color, linewidth=3)
            
            
        # Draw the original points, if provided
        box_marker = 's'
        if original_points:
            x, y = zip(*original_points)
            ax.scatter(x, y, marker=box_marker, color='grey', s=75, label='Car')
        
        # Draw the generated points, if provided
        if generated_points:
            x_gen, y_gen = zip(*generated_points)
            ax.scatter(x_gen, y_gen, color='b', marker='o', label='Generated Points')
        
        # Set x-axis formatter to display in 10e3 format
        ax.xaxis.set_visible(False)
        # Add x and y axes
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
        plt.title("KTH graph network with nodes and edges")
        plt.ylabel('in CRS units [m]')
        # Save the plot as an image
        self.image_cnt += 1
        # plt.savefig(f'image_{self.image_cnt}.png')
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
        
        
            
    def shortest_distance_along_roads(self, point1, point2, is_crs):
        """
        Given two points with geo loc. and spd, calculate the shortest path and the length of the shortest path
        
        """
        if not hasattr(self, 'G'):
            self._build_network_graph()
        
        # find the nearest touching points (nodes or points on edges) on the road network to the given points
        distance1, touching_point1 = self.find_nearest_distance_and_point_to_graph(point1, is_crs)
        distance2, touching_point2 = self.find_nearest_distance_and_point_to_graph(point2, is_crs)
        
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
        
        return shortest_distance, shortest_path, touching_point1, touching_point2

        
        
    def find_nearest_distance_and_point_to_graph(self, point, is_crs):
        """
        Find the closest points in the road network to the given point

        """
        if not hasattr(self, 'G'):
            self._build_network_graph()
        
        if is_crs:
            point_crs = point
        else:
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
        if not isinstance(touching_point, Point):
            if isinstance(touching_point, tuple):
                touching_point = Point(touching_point)
            else:
                touching_point = Point(touching_point['geometry'].iloc[0].x, touching_point['geometry'].iloc[0].y)
        
        return nearest_distance, touching_point
    
    
    
    def is_path_crossing_point(self, path, point):
        """
        Check if the point is across the path

        """
        for u_coords, v_coords in zip(path[:-1], path[1:]):            
            if is_point_crossing_segment(point, u_coords, v_coords):
                return True
        return False
    
    
    
    def generate_points_along_path(self, pt1, pt2):
        """
        Generate points evenly distributed along the shortest path between two nodes with a point every meter.
    
        Parameters:
            pt1 (tuple): The coordinates of the starting point as a tuple (x, y).
            pt2 (tuple): The coordinates of the ending point as a tuple (x, y).
    
        Returns:
            list: A list of generated points as (x, y) tuples.
        """
        points = []
        distance, path, _, _ = self.shortest_distance_along_roads(pt1, pt2, is_crs=0)
        num_points = max(round(distance), 1)  # Ensure at least 1 point
        
        if len(path) > 1:
            for i in range(len(path) - 1):
                segment_length = self.G[path[i]][path[i + 1]]['weight']
                num_segment_points = round(segment_length)  # Points every meter
                
                # Distribute points along the segment at regular intervals
                segment_points = np.linspace(0, segment_length, num_segment_points)
                for t in segment_points:
                    # Interpolate between node coordinates based on the distance
                    x = (1 - t / segment_length) * path[i][0] + (t / segment_length) * path[i + 1][0]
                    y = (1 - t / segment_length) * path[i][1] + (t / segment_length) * path[i + 1][1]
                    points.append((x, y))
        
        # Plot the generated points
        # self.draw_network_fig(generated_points=points)
        
        return points
    
    
    
    def calculate_distances_and_ttcs(self, pts_trace1, pts_trace2, vel1, vel2, acc1, acc2, is_draw=1):
        """
        Calculate distances, paths, and time-to-collision (TTC) between pairs of points in pts_trace1 and pts_trace2.
    
        Parameters:
            pts_trace1 (list of tuples): List of points for trace 1.
            pts_trace2 (list of tuples): List of points for trace 2.
            vel1 (float): Velocity for trace 1.
            vel2 (float): Velocity for trace 2.
            acc1 (float): Acceleration for trace 1.
            acc2 (float): Acceleration for trace 2.
    
        Returns:
            list of tuples: List of tuples containing distances, paths, and time-to-collision (TTC) for each pair of points.
        """

        # Sample points along the traces at specified density
        sampled_pts_trace1 = pts_trace1[::vel1]
        sampled_pts_trace2 = pts_trace2[::vel2]
        # Initialize list to store results for each pair of points
        results = []
    
        # Iterate over each pair of points in pts_trace1 and pts_trace2
        for pt1, pt2 in zip(sampled_pts_trace1, sampled_pts_trace2):
            gdf_pt1 = gpd.GeoDataFrame(geometry=[Point(pt1)], crs='EPSG:3854')
            gdf_pt2 = gpd.GeoDataFrame(geometry=[Point(pt2)], crs='EPSG:3854')
            
            # Calculate distance and path between the current pair of points
            distance, shortest_path, touch_pt1, touch_pt2 = self.shortest_distance_along_roads(gdf_pt1, gdf_pt2, is_crs=1) # TODO: revise is draw
            
            # Calculate time-to-collision using the calculated distance and provided velocities and accelerations
            t_collision = time_to_collision(distance, vel1, vel2, acc1, acc2)
            
            # draw the shortest path and points on the road network
            if is_draw:
                line_color = value_to_color(t_collision)
                touch_points = [(touch_pt1.x, touch_pt1.y), (touch_pt2.x, touch_pt2.y)]
                self.draw_network_fig(original_points=touch_points, shortest_path=shortest_path, color=line_color)
            
            # Append the results to the list
            results.append((pt1, pt2, distance, t_collision))
    
        return results


    def risk_access(self, point1, point2, prev_point1, prev_point2, case, is_crs=0, is_draw=0):
        """
        
        
        """
        
        if case == 'with_time':
            # estimate the speed of two vehicles
            dis1, _, _, _ = self.shortest_distance_along_roads(prev_point1[0], point1[0])
            spd1 = dis1/(point1[1] - prev_point1[1])
            dis2, _, _, _ = self.shortest_distance_along_roads(prev_point2[0], point2[0])
            spd2 = dis2/(point2[1] - prev_point2[1])
        elif case == 'with_spd':
            spd1 = point1[1]
            spd2 = point2[1]
        else:
            print('error!')
        
        # calculate the distance between two vehicles
        shortest_distance, shortest_path = self.shortest_distance_along_roads(point1[0], point2[0])
        
        # TODO: revise this part, right now direction judgement
        t_collision = time_to_collision(shortest_distance, spd1, spd2)
        # print(t_collision)
        
        return shortest_distance, t_collision