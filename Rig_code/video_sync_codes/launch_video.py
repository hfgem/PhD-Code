#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 11:23:58 2023

This program will synchronize time to a given server (or computer) time and 
launch video recording software. The time of video recording will be stored to
a csv file

@author: Hannah Germaine
"""

import os, easygui, subprocess, csv
import ntplib
from time import ctime
import time
from datetime import datetime

# Ask the user for the directory to save the video files in
directory = easygui.diropenbox(msg='Select the directory to save the videos from this experiment', title='Select directory')
# Change to that directory
os.chdir(directory)
# Ask the user for the video filename
ask_loop = 1
while ask_loop == 1:
	filename = input('What is the name of your video recording? [any spaces will be removed]')
	filename = filename.replace(' ','_')
	try:
		len_filename = len(filename)
		if len_filename > 0:
			print("Selected name: " + filename)
			ask_loop = 0
		else:
			print("Please try again. Incorrect entry.")
	except:
		print("Please try again. Incorrect entry.")

#Grab current date/time from Raspberry Pi
raspberry_pi_ip = input('What is the IP address of the Raspberry Pi running deliveries?')
# ask_loop = 1
# while ask_loop == 1:
# 	raspberry_pi_ip = input('What is the IP address of the Raspberry Pi running deliveries?')
# 	response = os.system("ping -c 1 " + raspberry_pi_ip)
# 	#and then check the response...
# 	if response == 0:
# 		print(f"{raspberry_pi_ip} is reachable!")
# 		ask_loop = 0
# 	else:
# 		print(f"{raspberry_pi_ip} is down! Check why and try again.")

c = ntplib.NTPClient() #Set up ntp client
#Create a .csv file in the video location and add current date/time
csv_name = filename + '_times.csv'
with open(csv_name, 'w', newline='') as csvfile:
	spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
	response = c.request(raspberry_pi_ip)
	t = datetime.fromtimestamp(response.tx_time)
	cur_time = t.strftime("%a %b %d %H:%M:%S.%f")
	call_time = cur_time
	spamwriter.writerow(call_time)

subprocess.run(record_video)

with open(csv_name, 'a+', newline='') as csvfile:
	spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
	response = c.request(raspberry_pi_ip)
	t = datetime.fromtimestamp(response.tx_time)
	cur_time = t.strftime("%a %b %d %H:%M:%S.%f")
	end_time = cur_time
	spamwriter.writerow(end_time)
