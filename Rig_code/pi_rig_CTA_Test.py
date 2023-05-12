'''
pi_rig_CTA_Test.py contrains a user-friendly protocol for calling tastant deliveries 
in a CTA-induction and CTA-testing protocol.
'''

import pi_rig_functions as pf

#Set up variables
purpose_options = ['Clear lines', 'Calibrate open times', 'Run deliveries']
delivery_options = ['Test', 'CTA induction day', 'CTA test day']
outport_options = [23, 29, 31, 33, 35, 37]
intan_options = [24, 26, 32, 36, 38, 40]

repeat_loop = 1

while repeat_loop == 1: #Returns to main menu asking purpose
	#Ask the user the purpose of the current session
	purpose_answer = pf.single_option_loop('What is the purpose of your running this code?',purpose_options)
			
	#Run correct pipeline depending on selection
	if purpose_answer == 0:
		print("For clearout, each outport selected will be individually cleared, \
		in order, for the given number of seconds. \n")
		outports_selected = pf.list_loop('Which outports are you using? ', outport_options)
		dur_selected = pf.single_float_loop('How long should lines be cleared for? ')
		pf.clearout(outport_options[outports_selected],dur_selected)
		#Check if user wants to perform any more actions
		repeat_loop = pf.more_stuff_loop()
		
	elif purpose_answer == 1:
		print("For calibrating open times, you will be asked which line is being calibrated, \
		how long the open time is, and how many deliveries are being calibrated. \n")
		outport_selected = pf.single_option_loop('Which outport would you like to calibrate? ',outport_options)
		dur_selected = pf.single_float_loop('How long would you like to calibrate for? ')
		repeat_options = [1,3,5]
		num_repeats = pf.single_option_loop('How many deliveries? ',repeat_options)
		pf.calibrate(outport_options[outport_selected], opentime = dur_selected, repeats = repeat_options[num_repeats])
		#Check if user wants to perform any more actions
		repeat_loop = pf.more_stuff_loop()
	
	elif purpose_answer == 2: #If running deliveries, confirm which kind
		delivery_answer = pf.single_option_loop('What kind of delivery? ',delivery_options)
		
		if delivery_answer == 0: #Testing delivery
			print("Testing delivery selected. If multiple outports are provided, all will be tested (in order)")
			outports_selected = pf.list_loop('Which outports would you like to test? ', outport_options)
			dur_selected = pf.multi_float_loop('How long should lines be cleared for?\n'+'Enter comma-separated values for each of these outports:\n'+outport_options[outports_selected])
			iti_min = pf.single_int_loop('What is the iti minimum? (integer value)')
			iti_max = pf.single_int_loop('What is the iti maximum? (integer value)')
			num_trials = pf.single_int_loop('How many deliveries per outport? (single integer value)')
			pf.passive(outport_options[outports_selected], intan_options[outports_selected], dur_selected, iti_min, iti_max, num_trials)
		
		elif delivery_answer == 1: #CTA induction day
		
		
		elif delivery_answer == 2: #CTA test day
		
			

#After performing actions, ensure everything is cleared
pf.clearall()