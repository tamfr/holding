# coding: utf-8
# Holding entries via parametric equations

import numpy as np
import matplotlib.pyplot as plt

def normalize_deg(deg_to_norm):
	""" Ensures degrees are between 0 and 360. """
	add_360 = (deg_to_norm < 0) * 360
	sub_360 = (deg_to_norm > 360) * 360
	return deg_to_norm + add_360 - sub_360 

def plot_entry_orbit(
	hdg, 
	inbnd_crs,
	R, 
	v, 
	turn_rate, 
	wind_direction, 
	wind_vel
	):
	
	adjusted_hdg = normalize_deg(hdg-inbnd_crs+180)
	
	if 0 + 291*R <= adjusted_hdg <= 69 + 291*R:
		entry_type = 'Teardrop Entry'
		teardrop_hdg = normalize_deg(inbnd_crs + 180 + 30*(-1)**R)
		
		if hdg > teardrop_hdg:
			turn = hdg - (180 + 30*(-1)**R) - inbnd_crs
			R = (not R)*1*R		
		else:
			turn = 180 + 30*(-1)**R + inbnd_crs - hdg
			R = (not R)*1 + R
			
	elif 69 + 40*R < adjusted_hdg < 251 + 40*R:
		entry_type = 'Direct Entry'
		turn = (hdg - inbnd_crs)*(-1)**R + 180
		
	else:
		entry_type = 'Parallel Entry'
		turn = (180 - hdg)*(-1)**R - inbnd_crs
		R = (not R)*1
	
	hdg_xform = - hdg * np.pi / 180
	
	angular_vel =  turn_rate * np.pi / 180 # rad/s
	wind_xform = - wind_direction + 3 * np.pi / 2
	
	r = v / angular_vel
	turn = normalize_deg(turn)
	time_to_turn = (turn / turn_rate)
	t = np.arange(0, time_to_turn, np.pi / 50)

	wind_x = wind_vel * np.cos(wind_xform)
	wind_y = wind_vel * np.sin(wind_xform)
	
	phase = hdg_xform + R * np.pi
	
	# Shift starting point to origin
	origin_x = r * np.cos(phase)
	origin_y = r * np.sin(phase)
	
	theta = angular_vel * t * (-1)**R
	
	x = r * np.cos(theta + phase) + wind_x * t - origin_x
	y = r * np.sin(theta + phase) + wind_y * t - origin_y
	 
	
	inbnd_crs_xform = (-inbnd_crs + 180) * np.pi / 180 + np.pi / 2
	
	r = np.arange(0, 6000, 10)
	
	x_crs = r*np.cos(inbnd_crs_xform)
	y_crs = r*np.sin(inbnd_crs_xform)
	
	ax = plt.subplot(111, polar=False)
	ax.plot(x,y)
	ax.plot(x_crs,y_crs)
	ax.grid(True)
	ax.set(aspect=1, adjustable='datalim')
	ax.set_title(entry_type, va='bottom')
	plt.show()


hdg = 30
R = 1
inbnd_crs = 360

v = 128 # m/s
turn_rate = 3 # degrees/s

wind_vel = 0
wind_direction = 180 * np.pi / 180 # rad

entry_orbit = plot_entry_orbit(
	hdg, 
	inbnd_crs,
	R,
	v,
	turn_rate,
	wind_direction,
	wind_vel
	)
	

