# Author: Thijs van 't Riet
# Date: 11-5-2024
# Description: Data analysis script for BCA assay data performed in a 96-well plate format. 
# The script reads an Excel file containing the raw absorbance values, 
# calculates the protein concentration based on a BSA calibration curve, 
# and plots the results for different sample types (e.g., pure protein, cell-free extract, load).

# The script does not take into account if the absorbance exceeds the linear range of the BCA assay, or the range of the spectrophotometer.

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

# Filepath to the excel file of your BCA measurement
file_path = "H:\\file\\project\\Enzymes\\ADH\\20240521_BCA_ADH.xlsx"

# Enter your own parameters
direction = "columns" # Are the samples in columns or rows?
wavelength = 562 # Wavelength used for the BCA assay
Multiplicity = 2 # Number of samples per time point, so if it is in dulpicates, Multiplicity = 2
Sample_types = ["Pure", "CFE", "Pure1","CFE1"]
BSA_start_conc = 2 # Concentration of BSA calibration used for the BCA assay in g/L
BSA_dilution = 2 # Dilution factor of the BSA calibration used for the BCA assay
Dilution_factor = 2 # Dilution factor of the samples

# Read the Excel file
df = pd.read_excel(file_path)

# Slice the data
start = df[df.iloc[:, 0] == wavelength].index
end = df[df.iloc[:, 0] == "Results"].index
df =  df.iloc[(start[0]+2):(end[0]-1), 1:]
df.columns = df.iloc[0]
df = df.drop(df.index[0])

# Create dictionaries to hold the lists
abs_avg_dict = {}
abs_std_dict = {}
conc_avg_dict = {}
total_conc_avg_dict = {}

t = Multiplicity+1   # This is to initialize the column count for the regex pattern, it skips the first columns with BSA samples in it

# Define a function to calculate BSA calibration curve when samples are in columns 
def calculate_calibration_curve_columns(row, start_conc, dilution):
    BSA_abs_std = []
    BSA_conc = []
    BSA_abs_list = []

    for i in range(1, Multiplicity+1):
        pattern = r'[A-Z]' + str(i) + r'\b'
        filtered_df = row.filter(regex=pattern)
        filtered_df = filtered_df.dropna(axis=0)
        BSA_abs = filtered_df.values.tolist()
        BSA_abs_list.append(BSA_abs)
        
    BSA_abs_avg = np.mean(BSA_abs_list, axis=0)
    BSA_abs_std = np.std(BSA_abs_list, axis=0)
    BSA_conc = [start_conc / (dilution ** i) for i in range(len(BSA_abs_avg) - 1)] + [0] # The last value is 0 because the last BSA concentration is 0, -1 to remove the last value
    return BSA_abs_avg, BSA_abs_std, BSA_conc

# Define a function to calculate BSA calibration curve when the samples are in rows 
def calculate_calibration_curve_rows(row, start_conc, dilution):
    BSA_abs_std = []
    BSA_conc = []
    BSA_abs_list = []

    for i in range(1, Multiplicity+1):
        letter = chr(64 + i)  # Convert the number to its respective letter in the alphabet
        pattern = r'' + letter + r'\d+\b'
        filtered_df = row.filter(regex=pattern)
        filtered_df = filtered_df.dropna(axis=0)
        BSA_abs = filtered_df.values.tolist()
        BSA_abs_list.append(BSA_abs)
        
    BSA_abs_avg = np.mean(BSA_abs_list, axis=0)
    BSA_abs_std = np.std(BSA_abs_list, axis=0)
    BSA_conc = [start_conc / (dilution ** i) for i in range(len(BSA_abs_avg) - 1)] + [0] # The last value is 0 because the last BSA concentration is 0, -1 to remove the last value
    return BSA_abs_avg, BSA_abs_std, BSA_conc

# Initialize lists for each sample type
for sample_type in Sample_types:
    abs_avg_dict[sample_type] = []
    abs_std_dict[sample_type] = []
    conc_avg_dict[sample_type] = []
    total_conc_avg_dict[sample_type] = []
    BSA_abs_avg_list = []
    BSA_abs_std_list = []
    time_list = []  # List to store the time values

    for index, row in df.iterrows():
        if direction == "columns":
            # Getting the values for the BSA calibration curve
            BSA_abs_avg, BSA_abs_std, BSA_conc = calculate_calibration_curve_columns(row, BSA_start_conc, BSA_dilution)
            BSA_abs_avg_list.append(BSA_abs_avg)
            BSA_abs_std_list.append(BSA_abs_std)
            time_list.append(row.iloc[0])  # Append the time value to the list

            # Fit a linear regression to the calibration curve
            regression_model = LinearRegression()
            regression_model.fit(np.array(BSA_abs_avg).reshape(-1, 1), np.array(BSA_conc).reshape(-1, 1))

            # Getting the values for the protein sample
            Sample_abs_list = []

            for i in range(t, Multiplicity+t):
                pattern = r'[A-Z]' + str(i) + r'\b'
                filtered_df = row.filter(regex=pattern)
                filtered_df = filtered_df.dropna(axis=0)
                Sample_abs = filtered_df.values.tolist()
                Sample_abs_list.append(Sample_abs)

        if direction == "rows":
            # Getting the values for the BSA calibration curve
            BSA_abs_avg, BSA_abs_std, BSA_conc = calculate_calibration_curve_rows(row, BSA_start_conc, BSA_dilution)
            BSA_abs_avg_list.append(BSA_abs_avg)
            BSA_abs_std_list.append(BSA_abs_std)
            time_list.append(row.iloc[0])  # Append the time value to the list

            # Fit a linear regression to the calibration curve
            regression_model = LinearRegression()
            regression_model.fit(np.array(BSA_abs_avg).reshape(-1, 1), np.array(BSA_conc).reshape(-1, 1))

            # Getting the values for the protein sample
            Sample_abs_list = []

            for i in range(t, Multiplicity+t):
                letter = chr(64 + i)  # Convert the number to its respective letter in the alphabet
                pattern = r'' + letter + r'\d+\b'
                filtered_df = row.filter(regex=pattern)
                filtered_df = filtered_df.dropna(axis=0)
                Sample_abs = filtered_df.values.tolist()
                Sample_abs_list.append(Sample_abs)

        Sample_abs_avg = np.mean(Sample_abs_list, axis=0)
        Sample_abs_std = np.std(Sample_abs_list, axis=0)
        Sample_conc = regression_model.predict(np.array(Sample_abs_avg).reshape(-1, 1))

        Purity = [100 / (Dilution_factor ** i) for i in range(len(Sample_abs_avg) - 1)] + [0]
        
        Sample_conc_flat = Sample_conc.flatten().tolist()
        total_conc_avg = [x * (100/y) for x, y in zip(Sample_conc_flat[:-1], Purity[:-1])]

        total_conc_avg_dict[sample_type].append(total_conc_avg)
        abs_avg_dict[sample_type].append(Sample_abs_avg)
        abs_std_dict[sample_type].append(Sample_abs_std)
        conc_avg_dict[sample_type].append(Sample_conc)
    t += Multiplicity # This is to go to the columns of the next sample

# A nice function where you can choose what to plot

def plot_data(plot_type):
    # Plot the BSA calibration curve
    if plot_type == "BSA":
        plt.style.use('ggplot')
        for i, BSA_abs_avg in enumerate(BSA_abs_avg_list):
            plt.errorbar(BSA_conc, BSA_abs_avg, yerr=BSA_abs_std_list[i], fmt='-o', label=time_list[i])  # Use the time value as the label
        plt.xlabel("BSA Concentration (g/L)")
        plt.ylabel("Absorbance")
        plt.title("BSA Calibration Curve")
        plt.legend()  # Add the legend
        plt.show()
    # Plot the protein concentration for different sample types
    elif plot_type == "protein":
        plt.style.use('ggplot')
        fig, axs = plt.subplots(1, len(Sample_types), figsize=(6 * len(Sample_types),8))

        for i, sample_type in enumerate(Sample_types):
            for j, abs_avg in enumerate(abs_avg_dict[sample_type]):
                axs[i].errorbar(conc_avg_dict[sample_type][j], abs_avg, yerr=abs_std_dict[sample_type][j], fmt='-o', label=time_list[j])  # Use the time value as the label
            axs[i].set_ylabel(f"{sample_type} protein concentration (g/L)")
            axs[i].set_xlabel("Absorbance")
            axs[i].set_title(f"{sample_type} protein")
            axs[i].legend()  # Add the legend

        plt.tight_layout()
        plt.show()
    # Plot the corrected total concentration for different sample types
    elif plot_type == "total_concentration":
        plt.style.use('ggplot')
        fig, axs = plt.subplots(1, len(Sample_types), figsize=(6 * len(Sample_types),8))

        for i, sample_type in enumerate(Sample_types):
            for j, abs_avg in enumerate(abs_avg_dict[sample_type]):
                axs[i].errorbar(abs_avg_dict[sample_type][j][:-1], total_conc_avg_dict[sample_type][j], fmt='-o', label=time_list[j])  # Use the time value as the label
            axs[i].set_xlabel(f"Absorbance")
            axs[i].set_ylabel("Corrected total Concentration (g/L)")
            axs[i].set_title(f"{sample_type} protein")
            axs[i].legend()  # Add the legend

        plt.tight_layout()
        plt.show()
    else:
        print("Invalid plot type")

plot_data("total_concentration") # You can choose for BSA, protein or total_concentration