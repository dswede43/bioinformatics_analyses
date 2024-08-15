#Helper functions - Alliance Genome data mining app
#---
#Helper functions to query keywords from the Alliance Genome API and return lists of genes associated to those keywords.

#Functions
#---
#1. api_request: function to make API request to Alliance Genome
#2. query_limit: function to return the limit and offsets for keyword query
#3. query_keyword: function to return a list of genes associated to keywords
#4. convert_df: function to convert data frame encoding for saving as CSV
#5. download_data_button: function to download the data frame from the app
#6. convert_hgnc_symbols: function to return the HGNC id's for a list of HGNC symbols
#7. query_gene_functions: function to return the functions for a list of genes


#import packages
import requests
import json
import pandas as pd
import streamlit as st

#define global variables
INCREMENT = 300


#1. Function to make API request to Alliance Genome
#---
def api_request(api_url, headers = None):
    try:
        #send API request
        response = requests.get(api_url, headers = headers)
        
        if response.status_code != 200:
            #raise http error with custom message
            response.raise_for_status()
        
        else:
            #format the data from the API request
            data = response.json()
        
    except requests.exceptions.HTTPError as http_err:
        #display specific HTTP error
        print(f"HTTP error occurred: {http_err}")
    
    #return the API request data
    return data


#2. Function to return the limit and offsets for keyword query
#---
def query_limit(keyword):
    #define the API url
    api_url = 'https://www.alliancegenome.org/api/search?category=gene&limit=1&q=' + keyword
    
    #obtain the API request data
    data = api_request(api_url)
    
    #define the limit for the current keyword
    limit = data['total']
    
    #create an array of offsets
    start = 0
    increment = INCREMENT
    offsets = [start]
    while start + increment < limit:
        start += increment
        offsets.append(start)
    offsets.append(limit - 1)
    
    #return the number of genes from each query
    return limit, offsets


#3. Function to return a list of genes associated to keywords
#---
def query_keyword(keyword, limit, offset, species = 'Homo sapiens'):
    #define a list of column names
    columns = ['symbol','name','id','soTermName']
    
    #define the API url
    api_url = 'https://www.alliancegenome.org/api/search?category=gene&limit=' + str(INCREMENT) + '&offset=' + str(offset) + '&q=' + keyword
    
    #obtain the API request data
    data = api_request(api_url)
    
    #return API request results
    results = data['results']
    
    #format API request data as data frame
    df = pd.json_normalize(results)
    
    #filter for species
    df = df[df['species'] == species]
    
    #subset the data frame
    df = df[columns]
    
    #create new column names
    df.columns = ['hgnc','name','id','biotype']
    
    #add the current keyword
    df['keyword'] = keyword
    
    #return the gene data
    return df


#4. Function to convert data frame encoding for saving as CSV
#---
def convert_df(df):
    #convert data frame encoding
    return df.to_csv().encode("utf-8")


#5. Function to download the data frame from the app
#---
@st.experimental_fragment
def show_download_button(csv):
    #create download button
    st.download_button(
    label = "Download data as CSV",
    data = csv,
    file_name = "keyword_gene_associations.csv",
    mime = "text/csv")


#6. Function to return the HGNC id's for a list of HGNC symbols
#---
def convert_hgnc_symbols(hgnc_symbols):
    #define an empty list and dictionary
    invalid_symbols = []
    hgnc_ids = {}
    for hgnc_symbol in hgnc_symbols:
        #define the API url
        api_url = 'https://rest.genenames.org/fetch/symbol/' + hgnc_symbol
        
        #define the API request headers
        headers = {'Accept': 'application/json'}
        
        #obtain the API request data
        data = api_request(api_url, headers)
        
        #if the API request has returned a result
        if data['response']['numFound'] > 0:
            #obtain the HGNC id and append
            hgnc_id = data['response']['docs'][0]['hgnc_id']
            hgnc_ids[hgnc_symbol] = hgnc_id
        else:
            #mark the HGNC symbol as invalid
            invalid_symbols.append(hgnc_symbol)
    
    #return the list of HGNC id's and list of invalid HGNC symbols
    return hgnc_ids, invalid_symbols


#7. Function to return the functions for a list of genes
#---
def query_gene_functions(hgnc_symbols):
    #convert HGNC symbols to HGNC id's
    hgnc_ids, invalid_symbols = convert_hgnc_symbols(hgnc_symbols)
    
    #define an empty dictionary
    gene_functions = {}
    for hgnc_symbol in hgnc_ids:
        #define the HGNC id
        hgnc_id = hgnc_ids[hgnc_symbol]
        
        #define the API url
        api_url = 'https://www.alliancegenome.org/api/gene/' + hgnc_id
        
        #obtain the API request data
        data = api_request(api_url)
        
        #append the genes function to the dictionary
        gene_functions[hgnc_symbol] = data['geneSynopsis']
    
    #return the gene functions
    return gene_functions

