from RoadNetClass import RoadNet


# initialization step
road_data = 'KTH_allroad.geojson'
KTHRoadNet = RoadNet(road_data)

# %% when only position and time step are received
point1 = [18.06574449138022, 59.35305240139645, 1]
prev_point1 = point1
point2 = [18.07019912436768, 59.34851018180592, 1]
prev_point2 = point2

while True:
    # point1 = recv_message(Volvo_Car)
    #        --> receive information from ego vehicle, point1 = [lon, lat, time]
    # point2 = recv_message(XSense)       
    #        --> receive information from Xsense sensor, here I assume only one vehicle point is received, 
    #        --> point2 = [lon, lat, time]
    
    # calculate risk with ttc (t_collision)
    distance, t_collision = KTHRoadNet.risk_access(point1, point2, case='with_time')


# %% when only position and velocity are received
point1 = [18.06574449138022, 59.35305240139645, 10]
point2 = [18.07019912436768, 59.34851018180592, 10]

while True:
    # point1 = recv_message(Volvo_Car)
    #        --> receive information from ego vehicle
    #        --> point1 = [lon, lat, spd]
    # point2 = recv_message(XSense)       
    #        --> receive information from Xsense sensor, here I assume only one vehicle point is received, 
    #        --> point2 = [lon, lat, spd]
    
    # calculate risk with ttc (t_collision)
    distance, t_collision = KTHRoadNet.risk_access(point1, point2, case='with_speed')