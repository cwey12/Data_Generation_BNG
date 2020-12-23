# this function takes an input of 'gene' and 'number of cells' and 'directory' and outputs a df with the time resolutions x data for the genes in the cells
import pandas as pd
import os
import math
import random
import matplotlib.pyplot as plt

#this function takes in a directory of gdat files and returns one list of dataframes, with each cell having their own df
def extract_cell_data(dir):

    list_of_cells = os.listdir(dir)  # dir is your directory path
    num_cells = len(list_of_cells)
    list_of_df = []

    for cell in list_of_cells:
        #reads the table of data
        df = pd.read_table(dir + '/' + cell, sep="\s+", engine='python')

        # post process to remove the initial '#' as a column and reformat data
        df_modified = df[df.columns[:-1]]
        df_modified.columns = df.columns[1:]

        #appends it to growing list of all cells
        list_of_df.append(df_modified)

    return list_of_df

#this function takes in the specific gene the user wants to visualize, and the fraction of cells from the raw data that are to be visualized.
# It returns a list of df for  randomly selected cells that follow the fraction of cells requested
def gene_data_extractor(gene, fraction, list_of_df):

    #calculate which cells will be taken into account
    num_cells = len(list_of_df)
    fraction_cells = math.ceil(fraction*num_cells)

    #generates list of which files will be read
    cells_selected = random.sample(range(0, num_cells), fraction_cells)
    cells_selected.sort()

    #initilize list of extracted cells; elements are df
    extracted_cells = []
    flag = True
    for i in cells_selected:

        #initializes temp df with each df selected
        full_cell_df = list_of_df[i]

        #checks to make sure gene is inside of the file
        genes_present_in_cell = list(full_cell_df.columns)
        if gene not in genes_present_in_cell:
            flag = False
            return [], [], flag

        extracted_cell_df = full_cell_df[['time', gene]].copy()
        extracted_cells.append(extracted_cell_df)

    return extracted_cells, cells_selected, flag

def plot_cells(list_of_df, gene, cells_selected):
    #cycles through the list of df and plots each one
    fig, ax = plt.subplots(figsize=(30,10))
    for cell in list_of_df:
        cell.plot(x='time', y=gene, ax=ax)
    ax.set_ylabel(gene + ' Expression')
    ax.set_xlabel('Time')
    cells_selected_for_legend = ['Cell ' + str(cell_num) for cell_num in cells_selected]
    ax.legend(cells_selected_for_legend, loc = 'best')
    ax.set_title(gene + ' for ' + str(len(cells_selected)) + ' cells vs time')
    plt.show()

def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper


#this function takes in an input of time point i and returns a data frame with that timepoint for all cells from simulation
def cell_timepoint_extract(timepoint, list_of_df):
    cell_data_list = []
    print('\nExtracting data for timepoint: ' + str(truncate(timepoint, 2)))
    for cell in list_of_df:
        cell_data_temp_list = []
        cell_data = cell.loc[cell['time'].astype(str).str[:4].astype(float) == truncate(timepoint, 2)]
        if cell_data.dropna().empty:
            print('empty!')
        else:
            cell_data_temp_list = list(cell_data.iloc[0])
        cell_data_list.append(cell_data_temp_list)

    return cell_data_list


