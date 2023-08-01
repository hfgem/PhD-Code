'''
pi_rig contrains basic functions for using the raspberry pi behavior and electrophysiology rig in the Katz Lab

These functions can be used directly via ipython in a terminal window or called by other codes
'''

# Import things for running pi codes
import time, random, easygui, os, csv
from math import floor
import RPi.GPIO as GPIO

# Import other things for video
from subprocess import Popen
import numpy as np

# Setup pi board
GPIO.setwarnings(False)
GPIO.cleanup()
GPIO.setmode(GPIO.BOARD)

# To empty taste lines


def clearout(outports=[7, 11, 12, 13, 16, 31, 32, 33, 35, 36, 37, 38, 40, 15, 5, 3, 29], dur=5):

#7 and 11 are for ortho odor vacuum, 32, 36, 38, 40 are rig 2 outports, 13 and 16 are ortho odor input
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
def calibrate(outports=[7, 11, 13, 12, 16, 23, 29, 31, 32, 33, 35, 36, 37, 38, 40], opentime=0.05, repeats=5):

    # Setup pi board GPIO ports
    GPIO.setmode(GPIO.BOARD)
    for i in outports:
        GPIO.setup(i, GPIO.OUT)

    # Open ports
    for rep in range(repeats):
        for i in outports:
            GPIO.output(i, 1)
        time.sleep(opentime)
        for i in outports:
            GPIO.output(i, 0)
        time.sleep(1)

    print('Calibration procedure complete.')


# Passive deliveries
def passive(outports=[37, 36, 38, 40, 32, 16, 18],
    intaninputs=[15, 19, 21, 23, 11, 12, 13], 
    opentimes=[0.01, 0.01, 0.01, 0.01, 0.01, 0.01], 
    itimin=22, itimax=22, trials=30):
	
	# Ask the user for the directory to save the video files in
    directory = easygui.diropenbox(msg='Select the directory to save the delivery times from this experiment.', title='Select directory')
	# Change to that directory
    os.chdir(directory)
	
    # Setup pi board GPIO ports
    GPIO.setmode(GPIO.BOARD)
    for i in outports:
        GPIO.setup(i, GPIO.OUT)
    for i in intaninputs:
        GPIO.setup(i, GPIO.OUT)

    # Set and radomize trial order
    tot_trials = len(outports) * trials
    count = 0
    trial_array = trials * list(np.arange(len(outports)))
    time_array = [] #Store delivery times
    random.shuffle(trial_array)

    time.sleep(15)
    print(trial_array)
    # Loop through trials
    for i in trial_array:
        time_array.append(time.ctime())
        GPIO.output(outports[i], 1) #opens the solenoid
        GPIO.output(intaninputs[i], 1) #changes dig_in signal to 1
        time.sleep(opentimes[i]) #waits for opentime
        GPIO.output(outports[i], 0) #closes the solenoid
        GPIO.output(intaninputs[i], 0) #changes dig_in signal to 0
        count += 1
        iti = random.randint(itimin, itimax)
        print('Trial '+str(count)+' of '+str(tot_trials) +
              ' completed. ITI = '+str(iti)+' sec.')
        time.sleep(iti)

    print('Passive deliveries completed')
	
	#Save csv of trial delivery times
    csv_name = 'output_times.csv'
    with open(csv_name, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=' ',
	                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for r_i in range(len(time_array)):
            spamwriter.writerow(time_array[r_i])
			
    print('Delivery times .csv saved.')
    

def passive_cue(
    outports=[7, 11, 13, 16, 31, 32, 33, 35, 36, 37, 38, 40], 
    intaninputs=[24, 26, 19, 21], 
    opentimes=[0.01], itimin=10, itimax=30, trials=150,
    cue_input = 40):

    # Setup pi board GPIO ports
    GPIO.setmode(GPIO.BOARD)
    for i in outports:
        GPIO.setup(i, GPIO.OUT)
    for i in intaninputs:
        GPIO.setup(i, GPIO.OUT)
    GPIO.setup(cue_input, GPIO.OUT)

    # Set and radomize trial order
    tot_trials = len(outports) * trials
    count = 0
    #trial_array = trials * range(len(outports))
    trial_array = trials * list(np.arange(len(outports)))
    time_array = [] #Store delivery times
    random.shuffle(trial_array)

    time.sleep(3)

    # Loop through trials
    for i in trial_array:
        time_array.append(time.ctime())
        GPIO.output(cue_input, 1)
        #time.sleep(1) #if you want cue on before taste delivery, uncomment
        #GPIO.output(cue_input, 0)
        #time.sleep(1)
        GPIO.output(outports[i], 1)
        GPIO.output(intaninputs[i], 1)
        time.sleep(opentimes[i])

        GPIO.output(outports[i], 0)
        GPIO.output(intaninputs[i], 0)
        time.sleep(1)
        GPIO.output(cue_input, 0)
        count += 1
        iti = random.randint(itimin, itimax)
        print('Trial '+str(count)+' of '+str(tot_trials) +
              ' completed. ITI = '+str(iti)+' sec.')
        time.sleep(iti)

    print('Passive deliveries completed')
	
	#Save csv of trial delivery times
    csv_name = 'output_times.csv'
    with open(csv_name, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=' ',
	                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for r_i in range(len(time_array)):
            spamwriter.writerow(time_array[r_i])
			
    print('Delivery times .csv saved.')


# Basic nose poking procedure to train poking for discrimination 2-AFC task
def basic_np(outport=40, opentime=0.012, iti=[.4, 1, 2], trials=200, outtime=0):

    intaninput = 8
    trial = 1
    inport = 13
    pokelight = 15 #edit
    houselight = 22 #edit
    lights = 0
    maxtime = 60

    # Setup pi board GPIO ports
    GPIO.setmode(GPIO.BOARD)
    GPIO.setup(pokelight, GPIO.OUT)
    GPIO.setup(houselight, GPIO.OUT)
    GPIO.setup(inport, GPIO.IN)
    GPIO.setup(outport, GPIO.OUT)
    GPIO.setup(intaninput, GPIO.OUT)

    time.sleep(15)
    starttime = time.time()

    while trial <= trials:

        # Timer to stop experiment if over 60 mins
        curtime = time.time()
        elapsedtime = round((curtime - starttime)/60, 2)
        if elapsedtime > maxtime:
            GPIO.output(pokelight, 0)
            GPIO.output(houselight, 0)
            break

        if lights == 0:
            GPIO.output(pokelight, 1)
            GPIO.output(houselight, 1)
            lights = 1

        # Check for pokes
        if GPIO.input(inport) == 0:
            poketime = time.time()
            curtime = poketime

            # Make rat remove nose from nose poke to receive reward
            while (curtime - poketime) <= outtime:
                if GPIO.input(inport) == 0:
                    poketime = time.time()
                curtime = time.time()

            # Taste delivery and switch off lights
            GPIO.output(outport, 1)
            GPIO.output(intaninput, 1)
            time.sleep(opentime)
            GPIO.output(outport, 0)
            GPIO.output(intaninput, 1)
            GPIO.output(pokelight, 0)
            GPIO.output(houselight, 0)
            print('Trial '+str(trial)+' of '+str(trials)+' completed.')
            trial += 1
            lights = 0

            # Calculate and execute ITI delay.  Pokes during ITI reset ITI timer.
            if trial <= trials/2:
                delay = floor((random.random()*(iti[1]-iti[0]))*100)/100+iti[0]
            else:
                delay = floor((random.random()*(iti[2]-iti[0]))*100)/100+iti[0]

            poketime = time.time()
            curtime = poketime

            while (curtime - poketime) <= delay:
                if GPIO.input(inport) == 0:
                    poketime = time.time()
                curtime = time.time()

    print('Basic nose poking has been completed.')

# Passive H2O deliveries


def affective(intaninputs=[24], tim_dur=1200):

    # Setup pi board GPIO ports
    GPIO.setmode(GPIO.BOARD)
    for i in intaninputs:
        GPIO.setup(i, GPIO.OUT)

    # Loop through trials

    GPIO.output(intaninputs[0], 1)
    time.sleep(0.1)
    GPIO.output(intaninputs[0], 0)

    time.sleep(tim_dur)

    GPIO.output(intaninputs[0], 1)
    time.sleep(0.1)
    GPIO.output(intaninputs[0], 0)

    print('Test completed')


# Clear all pi board GPIO settings
def clearall():

    # Pi ports to be cleared
    outports = [31, 33, 35, 37, 32, 36, 16, 18, 38, 40]
    inports = [11, 13, 15]
    pokelights = [15]
    houselight = 22
    lasers = [12, 18, 23, 8]
    intan = [12, 15, 19, 21, 23, 11, 13]
    ortho = [7, 11, 13, 16]

    # Set all ports to default/low state
    for i in intan:
        GPIO.setup(i, GPIO.OUT)
        GPIO.output(i, 0)

    for i in outports:
        GPIO.setup(i, GPIO.OUT)
        GPIO.output(i, 0)

    for i in inports:
        GPIO.setup(i, GPIO.IN, GPIO.PUD_UP)

    for i in pokelights:
        GPIO.setup(i, GPIO.OUT)
        GPIO.output(i, 0)

    for i in lasers:
        GPIO.setup(i, GPIO.OUT)
        GPIO.output(i, 0)

    GPIO.setup(houselight, GPIO.OUT)
    GPIO.output(houselight, 0)
