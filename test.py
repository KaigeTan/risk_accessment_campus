from RoadNetClass import RoadNet
from utils import time_to_collision


road_img = 'KTH_geofig.tif'
road_data = 'KTH_allroad.geojson'

KTHRoadNet = RoadNet(road_data, road_img)
KTHRoadNet.draw_road_fig()
KTHRoadNet.draw_network_fig()

# pt1: 18.06574449138022, 59.35305240139645
# pt2: 18.07019912436768, 59.34851018180592
# pt3: 18.070000661082304, 59.34911730809832
point1 = (18.06574449138022, 59.35305240139645)
point2 = (18.07019912436768, 59.34851018180592)
vel1 = 10
vel2 = 5
acc1 = 1
acc2 = 1
distance, path, _, _ = KTHRoadNet.shortest_distance_along_roads(point1, point2, is_crs=0)
t_collision = time_to_collision(distance, vel1, vel2, acc1, acc2)

if t_collision is not None:
    print('distance is: {:.3f}, and time-to-collision is: {:.3f}'.format(distance, t_collision))
else:
    print("No collision possible")

