# SETUP INSTRUCTIONS:
On the Raspberry pi:
- Update local repository index: sudo apt-get update
- Install NTP on the pi: sudo apt-get install ntp
- Verify ntp installation/version number: sntp --version
- (OPTIONAL) You can modify the NTP configuration file to fetch time from a local server: sudo nano /etc/ntp.conf
- (OPTIONAL cont...) In the open NTP file, modify the pool ... iburst lines with a server from the NTP server list: https://support.ntp.org/bin/view/Servers/NTPPoolServers
- Restart the NTP server: sudo service ntp restart
- Check the status of the NTP service: sudo service ntp status
- (IF THERE'S A FIREWALL) Configure the firewall to allow the other computer to access the NTP server (opens port 123 for incoming traffic): sudo ufw allow from any to any port 123 proto udp
    
On the video recording computer:
To set the computer time to the Raspberry Pi's:
- Install ntpdate to be able to access Raspberry Pi's server time: sudo apt-get install ntpdate
- Modify the NTP server host file:
-- Open the file: sudo nano /etc/hosts
-- Add the Raspberry Pi's IP address and specify a hostname ex. 000.000.000.0  NTP-server-name
-- Exit the file with Ctrl+X and save by entery y.
- Check if the time is synchronized with the server (you will need to replace NTP-server-name with the name you specified in the host file): sudo ntpdate NTP-server-name
- Disable the systemd timesyncd service on the client: sudo timedatectl set-ntp off
- Install NTP: sudo apt-get install ntp
- Configure the ntp config file to use the Raspberry Pi as the new time server.
-- Open the file: sudo nano /etc/ntp.conf
-- Add the following line (again, replace NTP-server-name with the name you specified above): server NTP-server-name prefer iburst
-- Exit the file with Ctrl+X and save by entery y.
- Restart the NTP server: sudo service ntp restart

To set up Python to grab Raspberry Pi's time:
- Install ntplib for python: pip3 install ntplib
- Install easygui for python: pip3 install easygui

# CONTENTS:
## pi_rig.py
Functions to run the rig.
## launch_video.py
Code to grab the pi time, trigger video recording, and record the start time of the video.

#Need .mp4 format for above ^