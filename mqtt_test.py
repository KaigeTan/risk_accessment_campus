from RoadNetClass import RoadNet


# initialization step
road_img = 'KTH_geofig.tif'
road_data = 'KTH_allroad.geojson'
KTHRoadNet = RoadNet(road_data)

# %% when only position and time step are received

# receive and parse data from mqtt message
# data_ego = recv_message(Volvo_Car)    --> receive information from ego vehicle
# data_ext = recv_message(XSense)       --> receive information from Xsense sensor
# veh_pos1 = parse_data_ego(data_ego)   --> get veh_pos1 = [lon, lat, time]
# veh_pos2 = parse_data_ext(data_ext)   --> get veh_pos2 = [lon, lat, time]

# calculate the velocity
# TODO: write vel calculation
case = 'with_time'
point1 = (18.06574449138022, 59.35305240139645)
prev_point1 = point1
point2 = (18.07019912436768, 59.34851018180592)
prev_point2 = point2

while True:
    # point1 = recv_message(Volvo_Car)    --> receive information from ego vehicle, point1 = [(lon, lat), time]
    # point2 = recv_message(XSense)       --> receive information from Xsense sensor, point2 = [(lon, lat), time]
    distance, _, t_collision = KTHRoadNet.risk_access(point1, point2, prev_point1, prev_point2, case)
    # update the points
    prev_point1 = point1
    prev_point2 = point2


# %% when only position and velocity are received