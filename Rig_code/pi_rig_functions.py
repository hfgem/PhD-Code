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


def single_option_loop(question, options):
	"""This function asks the user for input on a single selection from a list"""
	ask_loop = 1
	while ask_loop == 1:
		print(question)
		[print('\t' + str(i) + ': ' + options[i] + '\n') for i in range(len(options))]
		answer = input("Enter index of selection: ")
		try:
			int_answer = int(answer)
			if 0 <= int_answer <= len(options):
				print("Selected option: " + options[answer])
				ask_loop = 0
			else:
				print("Please try again. Incorrect entry.")
		except:
			print("Please try again. Incorrect entry.")
	
	return int_answer

def list_loop(question, options):
	"""This function asks the user for a list of inputs"""
	ask_loop = 1
	while ask_loop == 1:
		print(question)
		[print('\t' + str(i) + ': ' +  options[i] + '\n') for i in range(len(options))]
		answer = input("Enter indices from above list, comma separated: ")
		try:
			individual_selections = [int(o) for o in answer.split(',')]
			if len(individual_selections) > 0:
				print("Given answer: " + individual_selections)
				ask_loop = 0
			else:
				print("Please try again. Incorrect entry.")
		except:
			print("Please try again. Incorrect entry.")
	
	return answer

def single_int_loop(question):
	"""This function asks the user for an input of a float"""
	ask_loop = 1
	while ask_loop == 1:
		print(question)
		answer = input("Enter integer value: ")
		try:
			int_answer = int(answer)
			print("Given value: " + int_answer)
			ask_loop = 0
		except:
			print("Please try again. Incorrect entry.")
	
	return int_answer

def single_float_loop(question):
	"""This function asks the user for an input of a float"""
	ask_loop = 1
	while ask_loop == 1:
		print(question)
		answer = input("Enter float value: ")
		try:
			float_answer = float(answer)
			print("Given value: " + float_answer)
			ask_loop = 0
		except:
			print("Please try again. Incorrect entry.")
	
	return float_answer

def multi_float_loop(question):
	"""This function asks the user for a list of floats"""
	ask_loop = 1
	while ask_loop == 1:
		print(question)
		answer = input("Enter float values, comma separated: ")
		try:
			individual_selections = [float(o) for o in answer.split(',')]
			if len(individual_selections) > 0:
				print("Given answer: " + individual_selections)
				ask_loop = 0
			else:
				print("Please try again. Incorrect entry.")
		except:
			print("Please try again. Incorrect entry.")
	
	return individual_selections

def more_stuff_loop():
	"""This function asks the user if they want to keep running things or end 
	the session"""
	ask_loop = 1
	while ask_loop == 1:
		repeat_input = input("Would you like to perform any other actions [y/n]? ")
		if repeat_input != 'y' and repeat_input != 'n':
			print("Incorrect entry, please try again.")
		elif repeat_input == 'y':
			answer = 1
		elif repeat_input == 'n':
			answer = 0
			
	return answer
	
	

# To empty taste lines
def clearout(outports = [23, 29, 31, 33, 35, 37], dur = 5):

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
def calibrate(outport, opentime = 0.05, repeats = 5):

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
	
	
# Passive H2O deliveries
def passive(outports = [23, 29, 31, 33, 35, 37], intaninputs = [24, 26, 32, 36, 38, 40], opentimes = [0.01], itimin = 10, itimax = 30, trials = 150):


        # Setup pi board GPIO ports
	GPIO.setmode(GPIO.BOARD)
	for i in outports:
		GPIO.setup(i, GPIO.OUT)
	for i in intaninputs:
		GPIO.setup(i, GPIO.OUT)

        # Set and radomize trial order
		tot_trials = len(outports) * trials
		count = 0
		trial_array = trials * range(len(outports))
		random.shuffle(trial_array)

	time.sleep(15)
	
	# Loop through trials
	for i in trial_array:
		GPIO.output(outports[i], 1)
		GPIO.output(intaninputs[i], 1)
		time.sleep(opentimes[i])
		GPIO.output(outports[i], 0)
		GPIO.output(intaninputs[i], 0)
		count += 1
		iti = random.randint(itimin, itimax)
		print('Trial '+str(count)+' of '+str(tot_trials)+' completed. ITI = '+str(iti)+' sec.')
		time.sleep(iti)

	print('Passive deliveries completed')
	
# Passive deliveries with video recordings for Christina's CTA test protcol
def passive_with_video_CM(
	outports = [31, 33, 35], 
	intan_inports = [24, 26, 19], 
	tastes = ['water', 'sucrose', 'quinine'], 
	opentimes = [0.015, 0.015, 0.015], 
	iti = 20,
	session_one_indices = [0,1],
	session_two_indices = [0,1,2],
	session_three_indices = [0,2],
	two_taste_repeats = 10,
	three_taste_repeats = 10):

	# Set the outports to outputs
	GPIO.setmode(GPIO.BOARD)
	for i in outports:
		GPIO.setup(i, GPIO.OUT)

	# Set the input lines for Intan to outputs
	for i in intan_inports:
		GPIO.setup(i, GPIO.OUT)
		GPIO.output(i, 0)


	# Define the port for the video cue light, and set it as output
	video_cue = 16
	GPIO.setup(video_cue, GPIO.OUT)

	# Make an ordered array of session one tastes and randomize it
	session_one_taste_array = []
	for i in range(two_taste_repeats):
		session_one_taste_array.append(np.random.choice(session_one_indices, replace = False, size = len(session_one_indices)))
    
	session_one_taste_array = np.concatenate(session_one_taste_array)

	# Make an ordered array of session two tastes and randomize it
	session_two_taste_array = []
	for i in range(three_taste_repeats):
		session_two_taste_array.append(np.random.choice(session_two_indices, replace = False, size = len(session_two_indices)))
    
	session_two_taste_array = np.concatenate(session_two_taste_array)
        
    # Make an ordered array for session three tastes and randomize it
	session_three_taste_array = []
	for i in range(two_taste_repeats):
		session_three_taste_array.append(np.random.choice(session_three_indices, replace = False, size = len(session_three_indices)))
    
	session_three_taste_array = np.concatenate(session_three_taste_array)
    
    #create array that is length of all trials
	total_trials = two_taste_repeats * len(session_one_indices) + three_taste_repeats * len(session_two_indices) + two_taste_repeats * len(session_three_indices)
 
    # Ask the user for the directory to save the video files in	
	directory = easygui.diropenbox(msg = 'Select the directory to save the videos from this experiment', title = 'Select directory')
    # Change to that directory
	os.chdir(directory)

    # A 10 sec wait before things start
	time.sleep(10)

    # Session one: deliver water and CS only in random order
	trial_counter = [0 for i in range(len(outports))]
	for taste in session_one_taste_array:
        # Make filename, and start the video in a separate process
		process = Popen('sudo streamer -q -c /dev/video0 -s 1280x720 -f jpeg -t 180 -r 30 -j 75 -w 0 -o ' + tastes[taste] + '_trial_' + str(trial_counter[taste]) + '.avi', shell = True, stdout = None, stdin = None, stderr = None, close_fds = True)
    
        # Wait for 2 sec, before delivering tastes
		time.sleep(2)

        # Switch on the cue light
		GPIO.output(video_cue, 1)

        # Deliver the taste, and send outputs to Intan
		GPIO.output(outports[taste], 1)
		GPIO.output(intan_inports[taste], 1)
		time.sleep(opentimes[taste])	
		GPIO.output(outports[taste], 0)
		GPIO.output(intan_inports[taste], 0)

        # Switch the light off after 50 ms
		time.sleep(0.050)
		GPIO.output(video_cue, 0)

                
                # Increment the trial counter for the taste by 1
		trial_counter[taste] += 1                  
                                   
                
                # Print number of trials completed
		print("Trial " + str(np.sum(trial_counter)) + " of " + str(total_trials) + " completed.")

        # Wait for the iti before delivering next taste
		time.sleep(iti)
    
    # Session two: deliver all 3 tastes in random order
	for taste in session_two_taste_array:
        # Make filename, and start the video in a separate process
		process = Popen('sudo streamer -q -c /dev/video0 -s 1280x720 -f jpeg -t 180 -r 30 -j 75 -w 0 -o ' + tastes[taste] + '_trial_' + str(trial_counter[taste]) + '.avi', shell = True, stdout = None, stdin = None, stderr = None, close_fds = True)

                # Wait for 2 sec, before delivering tastes
		time.sleep(2)

		# Switch on the cue light
		GPIO.output(video_cue, 1)

		# Deliver the taste, and send outputs to Intan
		GPIO.output(outports[taste], 1)
		GPIO.output(intan_inports[taste], 1)
		time.sleep(opentimes[taste])	
		GPIO.output(outports[taste], 0)
		GPIO.output(intan_inports[taste], 0)

		# Switch the light off after 50 ms
		time.sleep(0.050)
		GPIO.output(video_cue, 0)

                
                # Increment the trial counter for the taste by 1
		trial_counter[taste] += 1                  
                                   
                
                # Print number of trials completed
		print("Trial " + str(np.sum(trial_counter)) + " of " + str(total_trials) + " completed.")

		# Wait for the iti before delivering next taste
		time.sleep(iti)

    # Session three: deliver water and quinine only in random order
	for taste in session_three_taste_array:
        # Make filename, and start the video in a separate process
		process = Popen('sudo streamer -q -c /dev/video0 -s 1280x720 -f jpeg -t 180 -r 30 -j 75 -w 0 -o ' + tastes[taste] + '_trial_' + str(trial_counter[taste]) + '.avi', shell = True, stdout = None, stdin = None, stderr = None, close_fds = True)

                # Wait for 2 sec, before delivering tastes
		time.sleep(2)

		# Switch on the cue light
		GPIO.output(video_cue, 1)

		# Deliver the taste, and send outputs to Intan
		GPIO.output(outports[taste], 1)
		GPIO.output(intan_inports[taste], 1)
		time.sleep(opentimes[taste])	
		GPIO.output(outports[taste], 0)
		GPIO.output(intan_inports[taste], 0)

		# Switch the light off after 50 ms
		time.sleep(0.050)
		GPIO.output(video_cue, 0)

                
                # Increment the trial counter for the taste by 1
		trial_counter[taste] += 1                  
                                   
                
                # Print number of trials completed
		print("Trial " + str(np.sum(trial_counter)) + " of " + str(total_trials) + " completed.")

		# Wait for the iti before delivering next taste
		time.sleep(iti)


# Passive deliveries with video recordings for Christina's CTA test protcol with 40 trials each of water, saccharin, quinine
def passive_with_video_CM_equal_trials(
	outports = [31, 33, 35], 
	intan_inports = [24, 26, 19], 
	tastes = ['water', 'sucrose', 'quinine'], 
	opentimes = [0.015, 0.015, 0.015], 
	iti = 20):

	# Set the outports to outputs
	GPIO.setmode(GPIO.BOARD)
	for i in outports:
		GPIO.setup(i, GPIO.OUT)

	# Set the input lines for Intan to outputs
	for i in intan_inports:
		GPIO.setup(i, GPIO.OUT)
		GPIO.output(i, 0)


	# Define the port for the video cue light, and set it as output
	video_cue = 16
	GPIO.setup(video_cue, GPIO.OUT)

	# Make an ordered array of session one tastes and randomize it
	session_one_tastes = [0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
	session_one_taste_array = np.random.choice(session_one_tastes, replace = False, size = len(session_one_tastes))
    

	# Make an ordered array of session two tastes and randomize it
	session_two_tastes = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
	session_two_taste_array = np.random.choice(session_two_tastes, replace = False, size = len(session_two_tastes))
    
        
    # Make an ordered array for session three tastes and randomize it
	session_three_tastes = [0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
	session_three_taste_array = np.random.choice(session_three_tastes, replace = False, size = len(session_three_tastes))
    
    
    #create array that is length of all trials
	total_trials = len(session_one_taste_array) + len(session_two_taste_array) + len(session_three_taste_array)
 
    # Ask the user for the directory to save the video files in	
	directory = easygui.diropenbox(msg = 'Select the directory to save the videos from this experiment', title = 'Select directory')
    # Change to that directory
	os.chdir(directory)

    # A 10 sec wait before things start
	time.sleep(10)

    # Session one: deliver water and CS only in random order
	trial_counter = [0 for i in range(len(outports))]
	for taste in session_one_taste_array:
        # Make filename, and start the video in a separate process
		process = Popen('sudo streamer -q -c /dev/video0 -s 1280x720 -f jpeg -t 180 -r 30 -j 75 -w 0 -o ' + tastes[taste] + '_trial_' + str(trial_counter[taste]) + '.avi', shell = True, stdout = None, stdin = None, stderr = None, close_fds = True)
    
        # Wait for 2 sec, before delivering tastes
		time.sleep(2)

        # Switch on the cue light
		GPIO.output(video_cue, 1)

        # Deliver the taste, and send outputs to Intan
		GPIO.output(outports[taste], 1)
		GPIO.output(intan_inports[taste], 1)
		time.sleep(opentimes[taste])	
		GPIO.output(outports[taste], 0)
		GPIO.output(intan_inports[taste], 0)

        # Switch the light off after 50 ms
		time.sleep(0.050)
		GPIO.output(video_cue, 0)

                
                # Increment the trial counter for the taste by 1
		trial_counter[taste] += 1                  
                                   
                
                # Print number of trials completed
		print("Trial " + str(np.sum(trial_counter)) + " of " + str(total_trials) + " completed.")

        # Wait for the iti before delivering next taste
		time.sleep(iti)
    
    # Session two: deliver all 3 tastes in random order
	for taste in session_two_taste_array:
        # Make filename, and start the video in a separate process
		process = Popen('sudo streamer -q -c /dev/video0 -s 1280x720 -f jpeg -t 180 -r 30 -j 75 -w 0 -o ' + tastes[taste] + '_trial_' + str(trial_counter[taste]) + '.avi', shell = True, stdout = None, stdin = None, stderr = None, close_fds = True)

                # Wait for 2 sec, before delivering tastes
		time.sleep(2)

		# Switch on the cue light
		GPIO.output(video_cue, 1)

		# Deliver the taste, and send outputs to Intan
		GPIO.output(outports[taste], 1)
		GPIO.output(intan_inports[taste], 1)
		time.sleep(opentimes[taste])	
		GPIO.output(outports[taste], 0)
		GPIO.output(intan_inports[taste], 0)

		# Switch the light off after 50 ms
		time.sleep(0.050)
		GPIO.output(video_cue, 0)

                
                # Increment the trial counter for the taste by 1
		trial_counter[taste] += 1                  
                                   
                
                # Print number of trials completed
		print("Trial " + str(np.sum(trial_counter)) + " of " + str(total_trials) + " completed.")

		# Wait for the iti before delivering next taste
		time.sleep(iti)

    # Session three: deliver water and quinine only in random order
	for taste in session_three_taste_array:
        # Make filename, and start the video in a separate process
		process = Popen('sudo streamer -q -c /dev/video0 -s 1280x720 -f jpeg -t 180 -r 30 -j 75 -w 0 -o ' + tastes[taste] + '_trial_' + str(trial_counter[taste]) + '.avi', shell = True, stdout = None, stdin = None, stderr = None, close_fds = True)

                # Wait for 2 sec, before delivering tastes
		time.sleep(2)

		# Switch on the cue light
		GPIO.output(video_cue, 1)

		# Deliver the taste, and send outputs to Intan
		GPIO.output(outports[taste], 1)
		GPIO.output(intan_inports[taste], 1)
		time.sleep(opentimes[taste])	
		GPIO.output(outports[taste], 0)
		GPIO.output(intan_inports[taste], 0)

		# Switch the light off after 50 ms
		time.sleep(0.050)
		GPIO.output(video_cue, 0)

                
                # Increment the trial counter for the taste by 1
		trial_counter[taste] += 1                  
                                   
                
                # Print number of trials completed
		print("Trial " + str(np.sum(trial_counter)) + " of " + str(total_trials) + " completed.")

		# Wait for the iti before delivering next taste
		time.sleep(iti)



	
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