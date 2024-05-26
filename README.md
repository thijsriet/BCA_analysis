# BCA_analysis
This script performs data analysis for BCA assay data performed in a 96-well plate format. 
It reads an Excel file containing the raw absorbance values, calculates the protein concentration based on a BSA calibration curve, and plots the results for different sample types. 
Depending on the export format you might have to adapt the slicing of the data. 
In this code it is based on the export formatting of the BioTek synergy 2 reader where the data looks like this:

Time	    T° 562	A1	  A2	  A3	  A4	  A5	  A6	  A7    ...
00:25:00	36,9	  0,513	1,461	2,42	1,764	0,65	1,031	0,322 ...
00:30:00	37	    0,547	1,713	2,574	2,028	0,709	1,117	0,345 ...

## Prerequisites

- Python 3.x
- pandas
- matplotlib
- numpy
- scikit-learn

## Installation

1. Clone the repository or download the code files.
2. Install the required dependencies using pip:
    ```
    pip install pandas matplotlib numpy scikit-learn
    ```

## Usage

1. Update the file path to the Excel file containing the BCA measurement data:
    ```
    file_path = "H:\\file\\project\\Enzymes\\ADH\\20240521_BCA_ADH.xlsx"
    ```

2. Modify the parameters according to your specific experiment:
    - `direction`: Specify whether the samples are in columns or rows.
    - `wavelength`: Set the wavelength used for the BCA assay.
    - `Multiplicity`: Specify the number of samples per time point.
    - `Sample_types`: Provide a list of sample types.
    - `BSA_start_conc`: Set the concentration of BSA calibration used for the BCA assay in g/L.
    - `BSA_dilution`: Specify the dilution factor of the BSA calibration used for the BCA assay.
    - `Dilution_factor`: Set the dilution factor of the samples.

3. Run the script:
    ```
    python BCA_analysis.py
    ```

4. Choose the plot type to display:
    - "BSA": Plot the BSA calibration curve.
    - "protein": Plot the protein concentration for different sample types.
    - "total_concentration": Plot the corrected total concentration for different sample types.

## License

This project is licensed under the [MIT License](LICENSE).

## Author

- Thijs van 't Riet

## Date

- 11-5-2024
