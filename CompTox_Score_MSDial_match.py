## Code to add a toxicity score to Comptox Hazard Dashboard Output, then merge scores back to MS-Dial concentration profiles
## CompTox Hazard Dashboard Output Score Criteria: VH = 4. H = 3, M = 2. L = 1, I = 0 (Change in lines 27-40 if necessary)
## Created 07/30/2024 By Peter Bright
## For Tap Drinking Water NTA Project, Garcia-Jaramillo Lab at OSU 
################Input Start##################
#### Notate your data file names####
file = 'Human_Toxicity_example.csv'
conc_file = 'MSDIAL_Conc_example.csv'

threshold = 100 ###This is the QC average intensity threshold (CPS) that features must be over in order to make the final output.
ppm_threshold = 20 ###If comptox Monoisotopic MZ and measured MZ are available, this is the ppm threshold compounds must be under in order to make the final output.
################Input End##################
print('=========Start of Script=========\n')
#Import Required Libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
current_directory = os.getcwd()
print("Current Directory:", current_directory)

###Import data File
toxdata = pd.read_csv('{}'.format(file))

# Define function to map letters to scores and store scores in columns VH, H, M, L, and I
def map_score(letter):
    if letter == 'VH':
        return 4
    elif letter == 'H':
        return 3
    elif letter == 'M':
        return 2
    elif letter == 'L':
        return 1
    elif letter == 'I':
        return 0
    else:
        return 0  # Handle NaN values
toxdata['Total_CompTox_Score'] = toxdata.apply(lambda row: sum(map_score(letter) for letter in row), axis=1)

def VH_score(VH):
    if VH == 'VH':
        return 4
    else:
        return 0
toxdata['VH'] = toxdata.apply(lambda row: sum(VH_score(VH) for VH in row), axis=1)

def H_score(H):
    if H == 'H':
        return 3
    else:
        return 0
toxdata['H'] = toxdata.apply(lambda row: sum(H_score(H) for H in row), axis=1)

def M_score(M):
    if M == 'M':
        return 2
    else:
        return 0
toxdata['M'] = toxdata.apply(lambda row: sum(M_score(M) for M in row), axis=1)

def L_score(L):
    if L == 'L':
        return 1
    else:
        return 0
toxdata['L'] = toxdata.apply(lambda row: sum(L_score(L) for L in row), axis=1)

### Imports the file with concentrations for individual compounds
msdial_data = pd.read_csv(r'{}'.format(conc_file))
msdial_data.rename(columns={'Salen':'Salem'},inplace=True)
#msdial_data = msdial_data1.drop_duplicates(subset=['Metabolite name','INCHIKEY'], keep='first')


#Query the dataset to see how many direct matches (based on metabolite name) there are between MS-Dial output and toxicity information
list = toxdata['INCHIKEY'].to_list()
total=len(list)
### match the data and store in 'hits'
hits = msdial_data[msdial_data['INCHIKEY'].isin(list)]
first_slice = len(hits['INCHIKEY'])
tox_duplicates = hits[hits['INCHIKEY'].duplicated()]
print('\n{}'.format(first_slice),'of the total {} suspect compounds matched'.format(total),'including {} duplicates from the tox data\n'.format(len(tox_duplicates)))




#msdial_data=msdial_data.rename(columns={'CASRN':'CAS'})
Tox_and_Conc_data = pd.merge(msdial_data, toxdata[['Total_CompTox_Score','CAS','Name','INCHIKEY']], on='INCHIKEY', how='left')
print('\nResulting dataframe is {}'.format(len(Tox_and_Conc_data)),'rows after merging {} rows from MS-Dial output\n'.format(len(msdial_data)))
print('\n--Matched {} Human Tox Scores to the Concentration Profile Data\n'.format(first_slice))



lost_vals=toxdata[~toxdata['INCHIKEY'].isin(Tox_and_Conc_data['INCHIKEY'])]
lost_vals[['INCHIKEY','PREFERRED_NAME','Total_CompTox_Score']]
names = lost_vals['PREFERRED_NAME'].to_list()
print('\nThese compounds were not matched to concentration profile data:\n{}\n\n'.format(names))



# Define function to map adduct masses for mass analysis. Values sourced from waters.com : https://support.waters.com/KB_Chem/Other/WKB67428_What_are_common_adducts_in_ESI_Mass_Spectrometry
# Note: any [2M + x] adducts will not be calculated as the M mass will be needed.
def adduct_mass(letter):
    if letter == '[M+H]+':
        return 1.0078
    elif letter == '[M+NH4]+':
        return 18.0344
    elif letter == '[M-H]-':
        return -1.0078
    elif letter == '[M+Na]+':
        return 22.9898
    elif letter == '[M+Cl]-':
        return 34.9689
    elif letter == '[M-H2O-H]-':
        return 27.0027
    elif letter == '[M+FA-H]-':
        return 44.9971
    elif letter == '[M-H2O+H]+':
        return -17.0027
    elif letter == '[M+ACN+H]+':
        return 42.0344
    elif letter == '[M+H-H2O]+':
        return 27.0027
    elif letter == '[M+Na-2H]-':
        return 20.9736
    else:
        return 0  # Handle NaN values
    
Tox_and_Conc_data['Adduct Mass'] = Tox_and_Conc_data.apply(lambda row: sum(adduct_mass(letter) for letter in row), axis=1)
Tox_and_Conc_data.replace('', np.nan, inplace=True)
Tox_and_Conc_data.replace(' ', np.nan, inplace=True)
Tox_and_Conc_data['Absolute Mass Difference'] = abs((Tox_and_Conc_data['Reference m/z'].astype(float) - Tox_and_Conc_data['Adduct Mass'].astype(float))-Tox_and_Conc_data['MONOISOTOPIC_MASS'].astype(float))
Tox_and_Conc_data['PPM Error From CompTox'] = abs((Tox_and_Conc_data['Absolute Mass Difference']/Tox_and_Conc_data['MONOISOTOPIC_MASS'].astype(float))*1e6)
questionable = ['[2M+H]+', '[2M+Na]+']
Tox_and_Conc_data.loc[Tox_and_Conc_data['Adduct type'].isin(questionable), ['Absolute Mass Difference','PPM Error From CompTox']] = 'No Calc'
Tox_and_Conc_data.fillna('',inplace=True)

print('\nUnfiltered Toxicity and Cocentration Profile dataset length: {}\n'.format(len(Tox_and_Conc_data)))
Tox_and_Conc_data['QC average'] = Tox_and_Conc_data[['QC 1', 'QC 2', 'QC 3', 'QC 4']].mean(axis=1)

Tox_and_Conc_data_prioritized = Tox_and_Conc_data[(Tox_and_Conc_data['QC average'] > threshold)]

Tox_and_Conc_data_prioritized['PPM Error From CompTox'] = pd.to_numeric(Tox_and_Conc_data_prioritized['PPM Error From CompTox'], errors='coerce')

# Filter rows where 'PPM Error From CompTox' is less than the threshold or is NaN (i.e., blank in original data)
Tox_and_Conc_data_prioritized = Tox_and_Conc_data_prioritized[(Tox_and_Conc_data_prioritized['PPM Error From CompTox'].isna()) | (Tox_and_Conc_data_prioritized['PPM Error From CompTox'] < ppm_threshold)]

# Filter out rows containing 'Blank' in 'Spectrum reference file name'
# Tox_and_Conc_data_prioritized = Tox_and_Conc_data_prioritized[(3*Tox_and_Conc_data_prioritized['average blank']) < (Tox_and_Conc_data_prioritized['average sample'])]
# Tox_and_Conc_data_prioritized = Tox_and_Conc_data_prioritized[(3*Tox_and_Conc_data_prioritized['average ext. blank']) < (Tox_and_Conc_data_prioritized['average sample'])]
Tox_and_Conc_data_prioritized = Tox_and_Conc_data_prioritized[~Tox_and_Conc_data_prioritized['Spectrum reference file name'].str.contains('Blank|QC|ISTD', case=False, na=False)]
Tox_and_Conc_data_prioritized.sort_values(by=['Confidence level', 'QC average','average sample'], ascending=[True, True, False], inplace=True)
print('\nFiltered Toxicity and Concentration profile dataset length: {} \n'.format(len(Tox_and_Conc_data_prioritized)))

Tox_and_Conc_data_prioritized.to_csv(r'Priority Compounds sorted from High Confidence List EXAMPLE.csv',index=False)

print('Successful\nSee ouput csv: Priority Compounds sorted from High Confidence List for results.\n=========End of Script=========')