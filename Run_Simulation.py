#This script is meant to be run on HPC3
#Clark Wey 10/18/20

import os
import shutil
from datetime import datetime
import sys
from pathlib import Path

#checks to see if user inputted number of cells to simulate
if len(sys.argv) != 2:
    print('ERROR: Did not input number of cells')
    sys.exit()

#determine which system the scripts are running on (HPC3 or Windows)
#first finds the absolute path to BioNetGen (Should be included in project package)
#then finds the correct path to where sim data will be saved
os_name = sys.platform
if os_name == 'win32': #running on Clark's Laptop
   project_path = Path(__file__).parent.absolute()
   path_to_bng = os.path.join(project_path, 'BioNetGen-2.5.1')
   path_to_data_save = 'C:/Users/clark/OneDrive/Documents/Files/Read_Research/Research_Data/BioNetGen_Sim/Sims'
   
elif os_name == 'linux': #running on HPC3
   project_path = Path(__file__).parent.absolute()
   path_to_bng = os.path.join(project_path, 'BioNetGen-2.5.1')
   path_to_data_save = '/data/homezvol2/cwey/Research_Data/BioNetGen_Sim/Sims'

else:
   print("ERROR: Run type not recognized.")
   sys.exit()


sim_name = 'Pluripotency_Mouse_Change.bngl' #name of the BNG simulation
path_to_model = os.path.join(project_path, 'Models', sim_name)
current_wd = os.getcwd()


#establish new directory to store data
now = datetime.now()
dt_string = now.strftime("%H_%M_%S_%d-%m-%Y")
simulation_data_folder = 'Sim_' + dt_string

#temp foldername
simulation_data_folder = 'PM_10_1_0.1'

os.makedirs(os.path.join(path_to_data_save, simulation_data_folder)) #Estalish new data directory
os.mkdir(os.path.join(path_to_data_save, simulation_data_folder, "cdat")) #Estalish cdat directory under data directory
os.mkdir(os.path.join(path_to_data_save, simulation_data_folder, "gdat")) #Estalish gdat directory under data directory
os.mkdir(os.path.join(path_to_data_save, simulation_data_folder, "net")) #Estalish net directory under data directory


#number of cells wanted in simulation
num_cells = sys.argv[1]
#concatenates the command to be run 
command = r'perl ' + str(project_path) + r'/BioNetGen-2.5.1/BNG2.pl ' + str(path_to_model) #command that is to be run under terminal

#print('Location of gdat files are: ' + path_to_data_save + '/' + simulation_data_folder + '/gdat')
#suppress_output = r' > /dev/null'

#dynamically obtain raw file name
sim_name_raw = os.path.splitext(sim_name)[0]


for i in range(int(num_cells)):
   os.system(command) #Runs the simulation
   #move data files created and rename them
   shutil.move(current_wd + '/' + sim_name_raw + '.cdat', path_to_data_save + '/' + simulation_data_folder + '/cdat/' + str(i) + '.cdat')
   shutil.move(current_wd + '/' + sim_name_raw + '.gdat', path_to_data_save + '/' + simulation_data_folder + '/gdat/' + str(i) + '.gdat')
   shutil.move(current_wd + '/' + sim_name_raw + '.net', path_to_data_save + '/' + simulation_data_folder +  '/net/'  + str(i) + '.net')




