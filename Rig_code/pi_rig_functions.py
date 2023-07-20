#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 12:50:49 2023

@author: hannahgermaine
Functions to be used in Raspberry Pi delivery code.
"""

# Import things for running pi codes
import time
from math import floor
import random
import RPi.GPIO as GPIO
import numpy as np
import easygui
import os

# Import other things for video
from subprocess import Popen
GPIO.setwarnings(False)
GPIO.cleanup()
GPIO.setmode(GPIO.BOARD)	
	

# To empty taste lines
def clearout(outports, dur):

    # Setup pi board GPIO ports
	GPIO.setmode(GPIO.BOARD)
	for i in outports:
		GPIO.setup([int(i)], GPIO.OUT)
	#Deliver for given amount of time
	for i in outports:
		GPIO.output([int(i)], 1)
	time.sleep(dur)
	for i in outports:
		GPIO.output([int(i)], 0)

	print('Tastant line clearing complete.')
	
	
# To calibrate taste lines
def calibrate(outport, opentime, repeats):

        # Setup pi board GPIO ports
	GPIO.setmode(GPIO.BOARD)
	GPIO.setup([int(outport)], GPIO.OUT)

        # Open ports  
	for rep in range(repeats):
		GPIO.output([int(outport)], 1)
		time.sleep(opentime)
		GPIO.output([int(outport)], 0)
		time.sleep(1)

	print('Calibration procedure complete.')
	
	
def passive(outports, intaninputs, opentimes, itimin, itimax, trials, taste_names):
	"""
	INPUTS:
		- outports: a list of integers of tastant delivery lines
		- intaninputs: a list of integers of intan input lines associated with outports
		- opentimes: a list of floats of open times for the delivery lines
		- itimin: minimum number of seconds between deliveries
		- itimax: maximum number of seconds between deliveries
		- trials: a list of integers of how many trials per taste
	OUTPUT: tastant deliveries according to the input information.
	"""
	#Passive tastant deliveries
	num_outports = len(outports)

	# Setup pi board GPIO ports
	GPIO.setmode(GPIO.BOARD)
	for i in outports:
		GPIO.setup([int(i)], GPIO.OUT)
	for i in intaninputs:
		GPIO.setup([int(i)], GPIO.OUT)

	# Randomize trial order
	tot_trials = np.sum(trials)
	count = 0
	trial_array = []
	for i in range(len(outports)):
		trial_array.extend([i for j in range(trials[i])])
	random.shuffle(trial_array)

	time.sleep(15)
	iti = random.randint(itimin, itimax)
	
	# Loop through trials
	trial_counter = np.zeros(len(outports))
	for i in range(len(trial_array)):
		t_i = trial_array[i]
		GPIO.output([int(outports[t_i])], 1)
		GPIO.output([int(intaninputs[t_i])], 1)
		time.sleep(opentimes[t_i])
		GPIO.output([int(outports[t_i])], 0)
		GPIO.output([int(intaninputs[t_i])], 0)
		count += 1
		trial_counter[t_i] += 1
		print('Trial '+str(count)+' of '+str(tot_trials)+' completed. ITI = '+str(iti)+' sec.')
		iti = random.randint(itimin, itimax)
		time.sleep(iti)

	print("Total deliveries per tastant:")
	[print("\t" + taste_names[i] + ": " + str(int(trial_counter[i]))) for i in range(len(outports))]

# Passive deliveries with video recordings for Christina's CTA test protcol with 40 trials each of water, saccharin, quinine
def passive_with_video(outports, intaninputs, opentimes, itimin, itimax, trials, taste_names, video_cue, directory, segment_num):
	"""
	INPUTS:
		- outports: a list of integers of tastant delivery lines
		- intaninputs: a list of integers of intan input lines associated with outports
		- opentimes: a list of floats of open times for the delivery lines
		- itimin: minimum number of seconds between deliveries
		- itimax: maximum number of seconds between deliveries
		- trials: a list of integers of how many trials per taste
		- taste_names: a list of strings of each outport's taste name
		- video_cue: single integer of the video port
		- directory: directory to store video recordings
	OUTPUT: tastant deliveries according to the input information.
	"""

	# Set the outports to outputs
	GPIO.setmode(GPIO.BOARD)
	for i in outports:
		GPIO.setup([int(i)], GPIO.OUT)

	# Set the input lines for Intan to outputs
	for i in intaninputs:
		GPIO.setup([int(i)], GPIO.OUT)
		GPIO.output([int(i)], 0)

	# Define the port for the video cue light, and set it as output
	GPIO.setup(video_cue, GPIO.OUT)
	
	# Randomize trial order
	tot_trials = np.sum(trials)
	count = 0
	trial_array = []
	for i in range(len(outports)):
		trial_array.extend([i for j in range(trials[i])])
	random.shuffle(trial_array)

    # Change to that directory
	os.chdir(directory)

    # A 10 sec wait before things start
	time.sleep(10)
	iti = random.randint(itimin, itimax)

    # Session one: deliver water and CS only in random order
	trial_counter = np.zeros(len(outports))
	for i in trial_array:
        # Make filename, and start the video in a separate process
		process = Popen('sudo streamer -q -c /dev/video0 -s 1280x720 -f jpeg -t 180 -r 30 -j 75 -w 0 -o ' + 
				  '_segment_' + str(int(segment_num)) + '_' + taste_names[i] + '_trial_' + str(int(trial_counter[i])) + 
				  '.avi', shell = True, stdout = None, stdin = None, stderr = None, close_fds = True)

        # Switch on the cue light
		GPIO.output(video_cue, 1)

        # Deliver the taste, and send outputs to Intan
		t_i = trial_array[i]
		GPIO.output([int(outports[t_i])], 1)
		GPIO.output([int(intaninputs[t_i])], 1)
		time.sleep(opentimes[t_i])	
		GPIO.output([int(outports[t_i])], 0)
		GPIO.output([int(intaninputs[t_i])], 0)

        # Switch the light off after 50 ms
		time.sleep(0.050)
		GPIO.output(video_cue, 0)

                
                # Increment the trial counter for the taste by 1
		trial_counter[t_i] += 1    
                # Print number of trials completed
		print("Trial " + str(np.sum(trial_counter)) + " of " + str(tot_trials) + " completed.")
		iti = random.randint(itimin, itimax)

        # Wait for the iti before delivering next taste
		time.sleep(iti)
		
	print("Total deliveries per tastant:")
	[print("\t" + taste_names[i] + ": " + str(trial_counter[i])) for i in range(len(outports))]
	
# Clear all pi board GPIO settings
def clearall(outports,inports,video_port):
	
	# Set all ports to default/low state
	for i in outports:
		try: #In case in a rig without certain value
			GPIO.setup([int(i)], GPIO.OUT)
			GPIO.output([int(i)], 0)
		except:
			"do nothing"
	for i in inports:
		try: #In case in a rig without certain value
			GPIO.setup([int(i)], GPIO.OUT)
			GPIO.output([int(i)], 0)
		except:
			"do nothing"
	GPIO.setup(video_port,GPIO.OUT)
	GPIO.output(video_port,0)
		