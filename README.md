# Unscented Kalman Filter

Tools.cpp contains several simple tools to calculate RMSE and to send and receive data.
UKF.cpp provides values to the variables. These ones here work. You might get better results by tweeking them more. I got to a point of diminishing returns and stopped.

This code is compiled using cmake and make. It was tested by running against obj_pose_laser-radar-synthetic-intput.txt located in the Mercedes directory. This is a set of lidar and radar measurements of a bike that is 'orbiting' a car. (Think of a hostile bicyclist harassing a driver.) 

Much of this project was about getting the right values for the input noise variables. The goal was to get the final RMSE values below certain targets. My outputs were:
px: 0.0687 < 0.09
py: 0.0827 < 0.10
vx: 0.3368 < 0.40
vy: 0.2185 < 0.30


