'''
pi_rig_CTA_Test.py contrains a user-friendly protocol for calling tastant deliveries 
in a CTA-induction and CTA-testing protocol.
'''

import pi_rig_functions as pf
import user_input_functions as ui
import numpy as np
import time, tqdm, os, easygui
file_path = ('/').join(os.path.abspath(__file__).split('/')[0:-1])
os.chdir(file_path)

#Set up variables
purpose_options = ['Clear lines', 'Calibrate open times', 'Run deliveries']
delivery_options = ['Normal Deliveries', 'CTA induction day', 'CTA test day']

#Confirm rig information
print("=========================")
outport_options = ui.multi_int_loop('\nWhat are the GPIO outport indices associated with the tastant lines?')
intan_options = ui.multi_int_loop('\nWhat are the GPIO intan indices associated with the ' + str(len(outport_options)) + ' tastant lines?')
video = ui.single_option_loop('\nDoes this rig have video set up? / Do you wish to record video?',['Yes','No'])
if video == 0:
	video_port = ui.single_int_loop('\nWhat is the video cue light GPIO port?')
	# Ask the user for the directory to save the video files in	
	directory = easygui.diropenbox(msg = 'Select the directory to save the videos from this experiment', title = 'Select directory')
else:
	video_port = outport_options[0] #Just a placeholder value

repeat_loop = 1

while repeat_loop == 1: #Returns to main menu asking purpose
	#Ask the user the purpose of the current session
	print("=========================")
	purpose_answer = ui.single_option_loop('\nWhat is the purpose of your running this code?',purpose_options)
			
	#Run correct pipeline depending on selection
	if purpose_answer == 0:
		print("=========================")
		print("\nFor clearout, each outport selected will be individually cleared, \
		in order, for the given number of seconds. \n")
		outports_selected = ui.list_loop('\nWhich tastant lines (outports) are you using? ', outport_options)
		dur_selected = ui.single_float_loop('\nWhat opentime should be used for clearout? ')
		pf.clearout(np.array(outport_options)[outports_selected],dur_selected)
		#Check if user wants to perform any more actions
		repeat_loop = ui.more_stuff_loop()
		
	elif purpose_answer == 1:
		print("=========================")
		print("\nFor calibrating open times, you will be asked which line is being calibrated, \
		how long the open time is, and how many deliveries are being calibrated. \n")
		outport_selected = ui.single_option_loop('\nWhich single outport would you like to calibrate? ',outport_options)
		dur_selected = ui.single_float_loop('\nWhat is the calibration opentime? ')
		repeat_options = [1,3,5]
		num_repeats = ui.single_option_loop('\nHow many deliveries? ',repeat_options)
		pf.calibrate(np.array(outport_options)[outport_selected], opentime = dur_selected, repeats = repeat_options[num_repeats])
		#Check if user wants to perform any more actions
		repeat_loop = ui.more_stuff_loop()
	
	elif purpose_answer == 2: #If running deliveries, confirm which kind
		delivery_answer = ui.single_option_loop('\nWhat kind of delivery? ',delivery_options)
		outports_selected = ui.list_loop('\nWhich outports would you like to use? ', outport_options)
		taste_names = [input("\nWhat is the name of the tastant in outport " + str(outport_options[i]) + "? ") for i in outports_selected]
		dur_selected = ui.multi_float_loop('\nWhat are the opentimes?\n'+'Enter comma-separated values for each of these outports:\n'+str(np.array(outport_options)[outports_selected]))
		iti_min = ui.single_int_loop('\nWhat is the iti minimum? (integer value)')
		iti_max = ui.single_int_loop('\nWhat is the iti maximum? (integer value)')
		
		if delivery_answer == 0: #Testing delivery
			print("=========================")
			print("\nTesting delivery selected. If multiple outports are provided, all will be tested (in order)")
			num_trials = ui.multi_int_loop('Outport options: '+str(np.array(outport_options)[outports_selected])+'\nAssociated tastant options: '+str(np.array(taste_names)[outports_selected])+'\nHow many deliveries per outport?')
			outports_i = np.array(outport_options)[outports_selected]
			intan_i = np.array(intan_options)[outports_selected]
			tastants_i = np.array(taste_names)[outports_selected]
			trials_i = np.ones(len(outports_selected))*num_trials
			input("Press enter when ready to begin.")
			print("Beginning deliveries.")
			if video == 0:
				pf.passive_with_video(outports_i, intan_i, dur_selected, iti_min, iti_max, num_trials, tastants_i, video_port, directory, segment_num=0)
			else:
				pf.passive(outports_i, intan_i, dur_selected, iti_min, iti_max, num_trials, tastants_i)
		
		elif delivery_answer == 1: #CTA induction day
			print("=========================")
			print("CTA induction day selected. Tastant deliveries will begin after a specified wait time.")
			pre_wait_time = ui.single_int_loop("\nWhat is the wait time PRIOR TO tastant delivery start (minutes)? ")
			post_wait_time = ui.single_int_loop("\nWhat is the wait time FOLLOWING tastant delivery (minutes)? ")
			num_trials = ui.multi_int_loop('Outport options: '+str(np.array(outport_options)[outports_selected])+'\nAssociated tastant options: '+str(np.array(taste_names)[outports_selected])+'\nHow many deliveries per outport?')
			outports_i = np.array(outport_options)[outports_selected]
			intan_i = np.array(intan_options)[outports_selected]
			input("Press enter when ready to begin.")
			#Begin delivery protocol - pre-wait + deliveries + post-wait
			print("=========================")
			print("Since you selected CTA induction day, there will be a pre- wait time, tastant delivery, and a post- wait time.")
			print('Beginning Pre-Delivery Wait Time.')
			for i in tqdm.tqdm(range(pre_wait_time*6)):
				time.sleep(10)
			print("Beginning deliveries.")
			if video == 0:
				pf.passive_with_video(outports_i, intan_i, dur_selected, iti_min, iti_max, num_trials, taste_names, video_port, directory, segment_num=0)
			else:
				pf.passive(outports_i, intan_i, dur_selected, iti_min, iti_max, num_trials, taste_names)
			print('Beginning Post-Delivery Wait Time.')
			for i in tqdm.tqdm(range(post_wait_time*6)):
				time.sleep(10)
		
		elif delivery_answer == 2: #CTA test day
			print("=========================")
			print("CTA test day selected. Tastant deliveries will begin after a specified wait time.")
			pre_wait_time = ui.single_int_loop("\nWhat is the wait time PRIOR TO tastant delivery start (minutes)? ")
			post_wait_time = ui.single_int_loop("\nWhat is the wait time FOLLOWING tastant delivery (minutes)? ")
			print("Since you selected CTA test day, there will be three delivery segments.")
			segment_tastants = [ui.list_loop("\nWhich tastants are being delivered in part " + str(i+1) + "? ",taste_names) for i in range(3)]
			segment_num_trials = [ui.multi_int_loop('\nHow many deliveries per outport in part ' + str(i+1) + ' for tastants ' + str(np.array(taste_names)[segment_tastants[i]]) + '?') for i in range(3)]
			input("Press enter when ready to begin.")
			#Begin delivery protocol - pre-wait + deliveries + post-wait
			print("=========================")
			print('Beginning Pre-Delivery Wait Time.')
			for i in tqdm.tqdm(range(pre_wait_time*6)):
				time.sleep(10)
			for i in range(3):
				lines_i = segment_tastants[i]
				tastants_i = np.array(taste_names)[lines_i]
				outports_i = np.array(outport_options)[lines_i]
				intan_i = np.array(intan_options)[lines_i]
				dur_i = np.array(dur_selected)[lines_i]
				trial_num_i = segment_num_trials[i]
				#Print info for user
				print('\nPart ' + str(i+1))
				print('\toutports = ' + str(outports_i) + '\n\tintans = ' + str(intan_i) + 
		  '\n\tdur_selected = ' + str(dur_i) + '\n\titi_min = ' + str(iti_min) + 
		  '\n\titi_max = ' + str(iti_max) + '\n\tsegment_num_trials = ' + str(trial_num_i) + 
		  '\n\ttastants = ' + str(tastants_i))
				#Run deliveries
				print("Beginning deliveries.")
				if video == 0:
					pf.passive_with_video(outports_i, intan_i, dur_i, iti_min, iti_max, trial_num_i, tastants_i, video_port, directory, segment_num=i)
				else:
					pf.passive(outports_i, intan_i, dur_i, iti_min, iti_max, trial_num_i, tastants_i)
			print('Beginning Post-Delivery Wait Time.')
			for i in tqdm.tqdm(range(post_wait_time*6)):
				time.sleep(10)
		#Check if user wants to perform any more actions
		repeat_loop = ui.more_stuff_loop()

#After performing actions, ensure everything is cleared
pf.clearall(outport_options,intan_options,video_port)