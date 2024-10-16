#!/usr/bin/python3

#FastQC results aggregation
#---
#This script automates the aggregation of specific results from a fastQC analysis of multiple FASTQ files.


#Table of contents:
#---
#1. Unzip all fastQC data
#2. Per base sequence quality
#3. Quality score distribution
#4. Per sequence GC content
#5. Per base N content
#6. Sequence duplication levels
#7. Adapter content


#define global variables
DIR="/path/to/directory" #working directory

#import libraries
import os
import zipfile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


#Functions
#---
#define function to unzip all fastQC data
def unzip_data(dir):
    for file in os.listdir(dir):
        if file.endswith(".zip"):
            file_path = os.path.join(dir, file)
            
            #unzip the .zip file
            with zipfile.ZipFile(file_path, 'r') as zip_file:
                zip_file.extractall(dir)
                print(f"Extracted: {file}")

#define function to load in fastQC data
def load_fastQC_data(dir):
    #change directory
    os.chdir(dir)
    
    #define data variables
    data_block = []
    data_blocks = {}
    section = 0
    capture = False
    #open the fastQC data file
    with open("fastqc_data.txt", 'r') as file:
        for line in file:
            #tokenize each data line
            line = line.strip()
            
            #start capturing the current data block after section header
            if line.startswith("#"):
                capture = True
                header = line.lstrip("#").split("\t")  # Use as column headers
                data_block.append(header)
                continue
            
            #stop capturing for current data block
            if line.startswith(">>"):
                capture = False
                
                #convert data into data frame
                df = pd.DataFrame(data_block)
                
                #store the results
                data_blocks[section] = df
                
                #reset the data block
                data_block = []
                section += 1
                continue
            
            #capture data from current block
            if capture:
                data_block.append(line.split("\t"))
    
    #remove empty data frames
    data = {}
    section = 0
    for i in range(len(data_blocks)):
        df = data_blocks[i]
        if not df.empty:
            data[section] = df
            section += 1
    
    return data


#1. Unzip all fastQC data
#---
#unzip all fastQC data
for file in os.listdir(DIR):
    if file.endswith(".zip"):
        unzip_data(DIR)

#remove all .zip files
for file in os.listdir(DIR):
    if file.endswith(".zip"):
        os.remove(os.path.join(DIR, file))


#2. Per base sequence quality
#---
#list all file names
sample_names = []
for file in os.listdir(DIR):
    if os.path.isdir(os.path.join(DIR, file)):
        sample_names.append(file)

#create the plotting data
QC_data = pd.DataFrame(columns = ['Base','Mean','Sample'])
for sample_name in sample_names:
    #load the fastQC data
    data = load_fastQC_data(os.path.join(DIR, sample_name))
    
    #isolate the specific data
    df = data[2]
    df.columns = df.iloc[0]
    df = df[1:]
    df = df[['Base','Mean']]
    df['Sample'] = sample_name.replace('_fastqc', '')
    
    #store the results
    QC_data = pd.concat([QC_data, df])
    QC_data['Mean'] = pd.to_numeric(QC_data['Mean'])

#create the lineplot
plt.figure(figsize = (10, 6))
sns.lineplot(data = QC_data, x = 'Base', y = 'Mean', hue = 'Sample')
plt.ylim(0, 37)
plt.axhspan(0, 20, color = 'red', alpha = 0.2)
plt.axhspan(20, 28, color = 'yellow', alpha = 0.2)
plt.axhspan(28, 37, color = 'green', alpha = 0.2)
plt.title('per base sequence quality')
plt.xlabel('position in read (bp)')
plt.ylabel('mean quality score (Phred)')
plt.xticks(rotation = 45)
plt.savefig(os.path.join(DIR, "mean_read_quality_scores.png"), dpi = 300)


#3. Quality score distribution
#---
#create the plotting data
QC_data = pd.DataFrame(columns = ['Quality','Count','Sample'])
for sample_name in sample_names:
    #load the fastQC data
    data = load_fastQC_data(os.path.join(DIR, sample_name))
    
    #isolate the specific data
    df = data[3]
    df.columns = df.iloc[0]
    df = df[1:]
    df['Sample'] = sample_name.replace('_fastqc', '')
    
    #store the results
    QC_data = pd.concat([QC_data, df])
    QC_data['Count'] = pd.to_numeric(QC_data['Count'])
    QC_data['Quality'] = pd.to_numeric(QC_data['Quality'])

#create the lineplot
plt.figure(figsize = (10, 6))
sns.lineplot(data = QC_data, x = 'Quality', y = 'Count', hue = 'Sample')
plt.title('per sequence quality')
plt.xlabel('mean quality score (Phred)')
plt.ylabel('counts')
plt.savefig(os.path.join(DIR, "read_quality_score_distributions.png"), dpi = 300)


#4. Per sequence GC content
#---
#create the plotting data
QC_data = pd.DataFrame(columns = ['GC Content','Count','Sample'])
for sample_name in sample_names:
    #load the fastQC data
    data = load_fastQC_data(os.path.join(DIR, sample_name))
    
    #isolate the specific data
    df = data[5]
    df.columns = df.iloc[0]
    df = df[1:]
    df['Sample'] = sample_name.replace('_fastqc', '')
    
    #store the results
    QC_data = pd.concat([QC_data, df])
    QC_data['Count'] = pd.to_numeric(QC_data['Count'])
    QC_data['GC Content'] = pd.to_numeric(QC_data['GC Content'])

#create the lineplot
plt.figure(figsize = (10, 6))
sns.lineplot(data = QC_data, x = 'GC Content', y = 'Count', hue = 'Sample')
plt.title('per sequence GC content')
plt.xlabel('mean GC content (%)')
plt.ylabel('counts')
plt.savefig(os.path.join(DIR, "GC_content_distributions.png"), dpi = 300)


#5. Per base N content
#---
#create the plotting data
QC_data = pd.DataFrame(columns = ['Base','N-Count','Sample'])
for sample_name in sample_names:
    #load the fastQC data
    data = load_fastQC_data(os.path.join(DIR, sample_name))
    
    #isolate the specific data
    df = data[6]
    df.columns = df.iloc[0]
    df = df[1:]
    df['Sample'] = sample_name.replace('_fastqc', '')
    
    #store the results
    QC_data = pd.concat([QC_data, df])
    QC_data['N-Count'] = pd.to_numeric(QC_data['N-Count'])

#create the lineplot
plt.figure(figsize = (10, 6))
sns.lineplot(data = QC_data, x = 'Base', y = 'N-Count', hue = 'Sample')
plt.ylim(0, 100)
plt.title('per base N-content')
plt.xlabel('position in read (bp)')
plt.ylabel('percentange (%)')
plt.xticks(rotation = 45)
plt.savefig(os.path.join(DIR, "N_content_levels.png"), dpi = 300)


#6. Sequence duplication levels
#---
#create the plotting data
QC_data = pd.DataFrame(columns = ['Duplication Level','Percentage of total','Sample'])
for sample_name in sample_names:
    #load the fastQC data
    data = load_fastQC_data(os.path.join(DIR, sample_name))
    
    #isolate the specific data
    df = data[8]
    df.columns = df.iloc[1]
    df = df.iloc[2:]
    df['Sample'] = sample_name.replace('_fastqc', '')
    
    #store the results
    QC_data = pd.concat([QC_data, df])
    QC_data['Percentage of total'] = pd.to_numeric(QC_data['Percentage of total'])

#create the lineplot
plt.figure(figsize = (10, 6))
sns.lineplot(data = QC_data, x = 'Duplication Level', y = 'Percentage of total', hue = 'Sample')
plt.ylim(0, 100)
plt.title('sequence duplication levels')
plt.xlabel('sequence duplication level')
plt.ylabel('percentange (%)')
plt.xticks(rotation = 45)
plt.savefig(os.path.join(DIR, "sequence_duplication_levels.png"), dpi = 300)


#7. Adapter content
#---
#create the plotting data
QC_data = pd.DataFrame(columns = ['Position','Mean Adapter Count','Sample'])
for sample_name in sample_names:
    #load the fastQC data
    data = load_fastQC_data(os.path.join(DIR, sample_name))
    
    #isolate the specific data
    df = data[10]
    df.columns = df.iloc[0]
    df = df[1:]
    tmp = df.drop(columns = df.columns[0], axis = 1)
    tmp = tmp.apply(pd.to_numeric, axis = 1)
    mean_counts = tmp.mean(axis = 1)
    
    df = pd.concat([df['Position'], mean_counts], axis = 1)
    df['Sample'] = sample_name.replace('_fastqc', '')
    df.columns = ['Position','Mean Adapter Count','Sample']
    
    #store the results
    QC_data = pd.concat([QC_data, df])

#create the lineplot
plt.figure(figsize = (10, 6))
sns.lineplot(data = QC_data, x = 'Position', y = 'Mean Adapter Count', hue = 'Sample')
plt.ylim(0, 100)
plt.title('mean adapter content')
plt.xlabel('position in read (bp)')
plt.ylabel('mean percentange (%)')
plt.xticks(rotation = 45)
ticks = plt.gca().get_xticks()
plt.gca().set_xticks(ticks[::2])
plt.savefig(os.path.join(DIR, "mean_adapter_content.png"), dpi = 300)
