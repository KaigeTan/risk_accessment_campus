from RoadNetClass import RoadNet
from utils import time_to_collision


road_img = 'KTH_geofig.tif'
road_data = 'KTH_allroad.geojson'

KTHRoadNet = RoadNet(road_data, road_img)
# KTHRoadNet.draw_road_fig()
# KTHRoadNet.draw_network_fig()

# pt1: 18.06574449138022, 59.35305240139645
# pt2: 18.07019912436768, 59.34851018180592
# pt3: 18.070000661082304, 59.34911730809832
point1 = (18.070000661082304, 59.34911730809832)
point2 = (18.07019912436768, 59.34851018180592)

ori_pt1_trace1 = (18.067205655966372, 59.35265551322761)
# ori_pt2_trace1 = (18.069985681266857, 59.34961492199544)
ori_pt2_trace1 = (18.070189367461456, 59.349340449537394)
ori_pt1_trace2 = (18.06752642321192, 59.3508789146111)
ori_pt2_trace2 = (18.069776592417707, 59.349110427936)

pts_trace1 = KTHRoadNet.generate_points_along_path(ori_pt1_trace1, ori_pt2_trace1)
pts_trace2 = KTHRoadNet.generate_points_along_path(ori_pt1_trace2, ori_pt2_trace2)

vel1 = 14
vel2 = 8
acc1 = 0
acc2 = 0

results = KTHRoadNet.calculate_distances_and_ttcs(pts_trace1, pts_trace2, vel1, vel2, acc1, acc2)
