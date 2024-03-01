from RoadNetClass import RoadNet


road_img = 'KTH_geofig.tif'
road_data = 'KTH_allroad.geojson'
#road_data = 'KTH_path_test.geojson'

KTHRoadNet = RoadNet(road_data, road_img)
KTHRoadNet.trim_road()
KTHRoadNet.draw_road_fig()
KTHRoadNet._build_network_graph()
KTHRoadNet.draw_network_fig()

# 18.06574449138022, 59.35305240139645
# 18.07019912436768, 59.34851018180592
# 18.070000661082304, 59.34911730809832
point1 = (18.070000661082304, 59.34911730809832)
point2 = (18.07019912436768, 59.34851018180592)
shortest_distance, path = KTHRoadNet.shortest_distance_along_roads(point1, point2)
nearest_distance1, touching_point1 = KTHRoadNet.find_nearest_distance_and_point_to_graph(point1)
nearest_distance2, touching_point2 = KTHRoadNet.find_nearest_distance_and_point_to_graph(point2)


distance = touching_point1.distance(touching_point2)
print('distance is :', distance)