import pandas as pd
#from tkinter import filedialog
import sys
from Data_formatter import *
import numpy as np
from datetime import datetime

#This script takes in an upper bound time and lower bound time, and selects n random timepoints

# Extracts limits from the command line arguements
lower_time_limit = int(sys.argv[1])
upper_time_limit = int(sys.argv[2])
num_timepoints = int(sys.argv[3])

#In case user inputs a lower limit that is higher than the upper limit
while lower_time_limit > upper_time_limit:
    print('\nERROR: Your lower limit was higher than the upper limit!')
    sys.exit()

#prompt user to select which simulation run to choose from
# print('Select the folder with the gdat files')
# gdat_dir = filedialog.askdirectory()
# print('\nExtracting files ...')
gdat_dir = sys.argv[4]

#calls the extract_cell_data function which returns a list of all the df
list_of_df = extract_cell_data(gdat_dir)

#calls the function to extract timepoint 'i' from each cell
columns = list(list_of_df[0].columns.values)
final_cell_data = []

for i in np.linspace(lower_time_limit, upper_time_limit, num_timepoints):
    data_at_timepoint_i = cell_timepoint_extract(i, list_of_df)
    final_cell_data.extend(data_at_timepoint_i)

final_cell_data_df = pd.DataFrame(final_cell_data, columns=columns)

#swtich the order of the dataframe columns so that time appears last
cols = final_cell_data_df.columns.values.tolist()
cols.append(cols.pop(cols.index('time')))
final_cell_data_df = final_cell_data_df.reindex(columns=cols)

#pwd to where data should be saved : /home/clarkwey/Documents/Read_Research/SINCERITIES_toolbox/BNG_Data
print('\nSaving Data to csv')
now = datetime.now()
dt_string = now.strftime("%d-%m-%Y_%H:%M:%S")
full_dt_string = 'Sim_' + dt_string
file_name = full_dt_string + '_' + str(len(list_of_df)) + '_cells_' + str(num_timepoints) + '_timepoints_from_' + str(lower_time_limit) + '_to_' + str(upper_time_limit) + '.csv'
save_directory = r'C:\Users\clark\OneDrive\Documents\Files\Read_Research\Research_Data\BioNetGen_Sim\Processed_Data'
final_cell_data_df.to_csv(save_directory + '/' + file_name, index=False)
