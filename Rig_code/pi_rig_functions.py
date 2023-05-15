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
		GPIO.setup(i, GPIO.OUT)
 
	for i in outports:
		GPIO.output(i, 1)
	time.sleep(dur)
	for i in outports:
		GPIO.output(i, 0)

	print('Tastant line clearing complete.')
	
	
# To calibrate taste lines
def calibrate(outport, opentime, repeats):

        # Setup pi board GPIO ports
	GPIO.setmode(GPIO.BOARD)
	GPIO.setup(outport, GPIO.OUT)

        # Open ports  
	for rep in range(repeats):
		GPIO.output(outport, 1)
		time.sleep(opentime)
		GPIO.output(outport, 0)
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
		GPIO.setup(i, GPIO.OUT)
	for i in intaninputs:
		GPIO.setup(i, GPIO.OUT)

	# Randomize trial order
	tot_trials = np.sum(trials)
	count = 0
	trial_array = []
	for i in range(len(outports)):
		trial_array.extend([i for j in range(trials[i])])
	random.shuffle(trial_array)

	time.sleep(15)
	
	# Loop through trials
	trial_counter = np.zeros(len(outports))
	for i in range(len(trial_array)):
		GPIO.output(outports[i], 1)
		GPIO.output(intaninputs[i], 1)
		time.sleep(opentimes[i])
		GPIO.output(outports[i], 0)
		GPIO.output(intaninputs[i], 0)
		count += 1
		iti = random.randint(itimin, itimax)
		print('Trial '+str(count)+' of '+str(tot_trials)+' completed. ITI = '+str(iti)+' sec.')
		time.sleep(iti)

	print("Total deliveries per tastant:")
	[print("\t" + taste_names[i] + ": " + str(trial_counter[i])) for i in range(len(outports))]

# Passive deliveries with video recordings for Christina's CTA test protcol with 40 trials each of water, saccharin, quinine
def passive_with_video(outports, intaninputs, opentimes, itimin, itimax, trials, taste_names):
	"""
	INPUTS:
		- outports: a list of integers of tastant delivery lines
		- intaninputs: a list of integers of intan input lines associated with outports
		- opentimes: a list of floats of open times for the delivery lines
		- itimin: minimum number of seconds between deliveries
		- itimax: maximum number of seconds between deliveries
		- trials: a list of integers of how many trials per taste
		- taste_names: a list of strings of each outport's taste name
	OUTPUT: tastant deliveries according to the input information.
	"""

	# Set the outports to outputs
	GPIO.setmode(GPIO.BOARD)
	for i in outports:
		GPIO.setup(i, GPIO.OUT)

	# Set the input lines for Intan to outputs
	for i in intaninputs:
		GPIO.setup(i, GPIO.OUT)
		GPIO.output(i, 0)

	# Define the port for the video cue light, and set it as output
	video_cue = 16
	GPIO.setup(video_cue, GPIO.OUT)
	
	# Randomize trial order
	tot_trials = np.sum(trials)
	count = 0
	trial_array = []
	for i in range(len(outports)):
		trial_array.extend([i for j in range(trials[i])])
	random.shuffle(trial_array)

    # Ask the user for the directory to save the video files in	
	directory = easygui.diropenbox(msg = 'Select the directory to save the videos from this experiment', title = 'Select directory')
    # Change to that directory
	os.chdir(directory)

    # A 10 sec wait before things start
	time.sleep(10)

    # Session one: deliver water and CS only in random order
	trial_counter = np.zeros(len(outports))
	for i in trial_array:
        # Make filename, and start the video in a separate process
		process = Popen('sudo streamer -q -c /dev/video0 -s 1280x720 -f jpeg -t 180 -r 30 -j 75 -w 0 -o ' + taste_names[i] + '_trial_' + str(trial_counter[i]) + '.avi', shell = True, stdout = None, stdin = None, stderr = None, close_fds = True)
    
        # Wait for 2 sec, before delivering tastes
		time.sleep(2)

        # Switch on the cue light
		GPIO.output(video_cue, 1)

        # Deliver the taste, and send outputs to Intan
		GPIO.output(outports[i], 1)
		GPIO.output(intaninputs[i], 1)
		time.sleep(opentimes[i])	
		GPIO.output(outports[i], 0)
		GPIO.output(intaninputs[i], 0)

        # Switch the light off after 50 ms
		time.sleep(0.050)
		GPIO.output(video_cue, 0)

                
                # Increment the trial counter for the taste by 1
		trial_counter[i] += 1    
		iti = random.randint(itimin, itimax)
                # Print number of trials completed
		print("Trial " + str(np.sum(trial_counter)) + " of " + str(tot_trials) + " completed.")

        # Wait for the iti before delivering next taste
		time.sleep(iti)
		
	print("Total deliveries per tastant:")
	[print("\t" + taste_names[i] + ": " + str(trial_counter[i])) for i in range(len(outports))]
	
# Clear all pi board GPIO settings
def clearall():

	# Pi ports to be cleared
	outports = [23, 29, 31, 33, 35, 37]
	inports = [24, 26, 32, 36, 38, 40]
	pokelights = [36, 38, 40]
	houselight = 18
	lasers = [12, 22, 16]
	intan = [8, 10, 24, 26, 19, 21]
	
	# Set all ports to default/low state
	for i in intan:
		try: #In case in a rig without certain value
			GPIO.setup(i, GPIO.OUT)
			GPIO.output(i, 0)
		except:
			"do nothing"
	
	for i in outports:
		try: #In case in a rig without certain value
			GPIO.setup(i, GPIO.OUT)
			GPIO.output(i, 0)
		except:
			"do nothing"
		
	for i in inports:
		try: #In case in a rig without certain value
			GPIO.setup(i, GPIO.IN, GPIO.PUD_UP)
		except:
			"do nothing"
		
	for i in pokelights:
		try: #In case in a rig without certain value
			GPIO.setup(i, GPIO.OUT)
			GPIO.output(i, 0)
		except:
			"do nothing"
		
	for i in lasers:
		try: #In case in a rig without certain value
			GPIO.setup(i, GPIO.OUT)
			GPIO.output(i, 0)
		except:
			"do nothing"
		
	GPIO.setup(houselight, GPIO.OUT)
	GPIO.output(houselight, 0)	