# Rig_code
This repository contains a collection of python scripts for running tastant deliveries on the raspberry pi.
 
## long_paradigm_codes:
### pi_rig_CTA_Test.py
Code to run deliveries and CTA induction/testing protocols.
### pi_rig_functions.py
Functions used by pi_rig_CTA_Test.py in delivering tastants
### user_input_functions.py
Functions used by pi_rig_CTA_Test.py in getting user input.

## video_sync_codes:
### pi_rig.py
Functions to run the rig.
### launch_video.py
Code to grab the pi time, trigger video recording, and record the start time of the video.

NOTE: Will need to run <pip install ntplib> in the python environment on the terminal of the computer that is performing the video recording.