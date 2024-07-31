## Code to add concatenate MSDial software output, sort according to magnitude of QC average and Sample average intensities (CPS), output results csv, and output unknowns data to txt files for SIRUSID input
## This code is highly personalized. Please email brightpe@oregonstate.edu for help conforming to own data if needed. 
## Created 07/30/2024 By Peter Bright
## For Tap Drinking Water NTA Project, Garcia-Jaramillo Lab at OSU 
################Input Start##################
#### Notate your data file names####

hilic_neg_file = 'HILIC_Neg_short_example.csv'
hilic_pos_file = 'HILIC_Pos_short_example.csv'
rpcart_pos_file = 'RP_Cart_Pos_short_example.csv'
rpcart_neg_file = 'RP_Cart_Neg_short_example.csv'
rpplate_pos_file = 'RP_Plate_Pos_short_example.csv'
rpplate_neg_file = 'RP_Plate_Neg_short_example.csv'

#directory= 'path_to_your_target_directory_for_.txt_file_output/enter_here' #The output directory where you would like up to 50 txt files saved 
# .txt files will be for annotations with the name -Unknown- only, and will save with the adduct type, mass, and LC_type in the file name
# .txt file content is the m/z values with corresponding ms/ms data in the first column, and the intensity of these peaks in the second column.

threshold = 100 ###This is the QC average intensity threshold (CPS) that features must be over in order to make the final output. Likely not a sensitive parameter, as only the top features will be extracted to output
desired_number = 10 ###This is the desired number of features included in the final output (unknowns and w/o MS2: data)
ppm_threshold = 20 ###If comptox Monoisotopic MZ and measured MZ are available, this is the ppm threshold compounds must be under in order to make the final output.
################Input End##################
print('=========Start of Script=========\n')
#Import Required Libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from io import StringIO ###Package to work with Text files
current_directory = os.getcwd()
print("Current Directory:", current_directory)


hilic_neg = pd.read_csv(hilic_neg_file)
hilic_pos = pd.read_csv('{}'.format(hilic_pos_file))
rpcart_pos = pd.read_csv('{}'.format(rpcart_pos_file))
rpcart_neg = pd.read_csv('{}'.format(rpcart_neg_file))
rpplate_pos = pd.read_csv('{}'.format(rpplate_pos_file))
rpplate_neg = pd.read_csv('{}'.format(rpplate_neg_file))


##Add a sample preparation identifier
hilic_neg['LC_Type'] = 'HILIC Neg'
hilic_pos['LC_Type'] = 'HILIC Pos'
rpcart_pos['LC_Type'] = 'RP Cart Pos'
rpcart_neg['LC_Type'] = 'RP Cart Neg'
rpplate_pos['LC_Type'] = 'RP Plate Pos'
rpplate_neg['LC_Type'] = 'RP Plate Neg'


###Merge the raw MSDIAL output together from the HILIC, RP Plate, and RP Cartridge datasets (note, this filters according to the stipulations below: QC > 100, Reference M/Z = null, MS/MS assigned = True)
qc_columns = ['QC 1', 'QC 2', 'QC 3', 'QC 4']

def filter_dataframes(dataframes):
    filtered_dataframes = []
    for df in dataframes:

        df['QC average'] = df[qc_columns].mean(axis=1)
        # Filter rows where 'QC 3' column values are greater than threshold

        df = df[(df['QC average'] > threshold)]
        
        # Filter rows where 'Reference m/z' column values are not equal to 'null'
        df = df[df['Reference m/z'] != 'null']
        
        # Filter out rows where 'Ms/MS assigned' column values indicate FALSE
        df = df[df['MS/MS assigned'] == True]
        
        # Append the filtered DataFrame to the list
        filtered_dataframes.append(df)
    
    # Concatenate the filtered DataFrames into a single DataFrame
    final_df = pd.concat(filtered_dataframes, ignore_index=True)
    
    return final_df

# Apply filter
filtered_df = filter_dataframes([rpcart_pos, hilic_neg, hilic_pos, rpcart_neg, rpplate_pos, rpplate_neg])



# Filter out rows containing 'Blank' in 'Spectrum reference file name'
filtered_df = filtered_df[(3*filtered_df['average blank']) < (filtered_df['average sample'])]
filtered_df = filtered_df[(3*filtered_df['average ext. blank']) < (filtered_df['average sample'])]
filtered_df = filtered_df[~filtered_df['Spectrum reference file name'].str.contains('Blank|QC|ISTD', case=False, na=False)]
filtered_df.sort_values(by=['QC average','average sample'], ascending=False, inplace=True)

print("\nTotal number of unique MS-Dial Annotations (including w/o MS2:): {}".format(len(filtered_df['Metabolite name'].unique())), "\nTotal number of unique m/z values in MS-Dial output: {}\n".format(len(filtered_df['Average Mz'].unique())))


unannotated_mz_vals = filtered_df

def remove_annotated_duplicates(df, column_name):
    keep = df[column_name]=='Unknown'
    delete = df.duplicated(subset=[column_name], keep='first')
    function = ~keep & delete
    df = df[~function]
    return df

nonduplicates_and_mz = remove_annotated_duplicates(unannotated_mz_vals, 'Metabolite name')

# Sort the DataFrame by 'QC average' and 'average sample'
nonduplicates_and_mz.sort_values(by=['QC average', 'average sample'], ascending=False, inplace=True)

print("\nTotal number of unique m/z values in MS-Dial output excluding annotation duplicates: {}".format(len(nonduplicates_and_mz['Average Mz'].unique())))

# Determine the number of rows needed to capture 50 unique 'Metabolite name' values
unique_masses = set()  #Store only unique values
rows_required = 0 #start the count at zero

for index, row in nonduplicates_and_mz.iterrows():
    unique_masses.add(row['Average Mz'])
    rows_required += 1
    if len(unique_masses) >= desired_number:
        break

print(f"Number of rows needed to capture at least 50 unique 'Average Mz' values: {rows_required}")

# Select the top rows_required rows
highest_rows = nonduplicates_and_mz.head(rows_required)

# Find the most frequent 'Metabolite name' values within these top rows
nonduplicates_and_mz = highest_rows['Average Mz'].value_counts().head(rows_required).index.tolist()

# Filter the original DataFrame for the top 50 annotated metabolites
filtered_mzs = filtered_df[filtered_df['Average Mz'].isin(nonduplicates_and_mz)]



annotations = filtered_mzs[filtered_mzs['Metabolite name'] != 'Unknown']
annotations.sort_values(by='QC average', ascending=False).to_csv('{} annotations from the most abundant MS-Dial features EXAMPLE OUTPUT.csv'.format(len(annotations)))
unknowns = filtered_mzs[filtered_mzs['Metabolite name'] == 'Unknown']
unknowns.sort_values(by='QC average', ascending=False).to_csv('{} unknowns from the most abundant MS-Dial features EXAMPLE OUTPUT.csv'.format(len(unknowns)))


unique_hits_for_sirus = filtered_mzs[['Average Rt(min)', 'Average Mz', 'MS/MS spectrum','Adduct type', 'LC_Type', 'Metabolite name']]
unique_hits_for_sirus = unique_hits_for_sirus[unique_hits_for_sirus['Metabolite name']== 'Unknown']
unique_hits_for_sirus['name'] = unique_hits_for_sirus['Average Rt(min)'].astype(str) + '_' + unique_hits_for_sirus['Average Mz'].astype(str)
unique_hits_for_sirus.drop(['Average Rt(min)','Average Mz'], axis=1, inplace=True)
print("\nTop {} annotated compounds extracted for import into SIRUS ID".format(len(unique_hits_for_sirus)))

test = unique_hits_for_sirus
for index, row in test.iterrows():    #Initiate a loop with iterrows (yields index and row data for each row in df) itterates over each row in the dataframe.
    values = row.values[0].split(' ') #split the value in row [0] at each space '  ', store in new df called values
    values = [value.split(':') for value in values] # itterate through the values in your new df and further split them at the colon ';'
    
    # Write values to a StringIO object
    output = StringIO() #create storage location for values
    for value in values:
        output.write(' '.join(value) + '\n') #join each pair of split values (these are the values from the loop above) with a space ' ', and then joing each value with the new line character '\n' to create new columnns
    output.seek(0) #moves the cursor back to the first line,
    

    filename = '{}'.format(current_directory) + '/{}.txt'.format(test.loc[index, 'name'] +'_'+test.loc[index, 'LC_Type'] +'_'+test.loc[index, 'Adduct type']) # This line constructs the filename for the output text file using the 'name', 'LC-Type', and 'Adduct Type' value from the current row of test into the filename string.
    with open(filename, 'w') as f:  #This line opens the file specified by filename in write mode ('w'). The file is opened within a with statement, which ensures that the file is properly closed after writing.
        f.write(output.read()) # Inside the with block, you're reading the contents of the StringIO object output and writing them to the file opened earlier.

print('\n=========End of Script=========')