#Helper functions - Alliance Genome data mining app
#---
#Helper functions to query keywords from the Alliance Genome API and return lists of genes associated to those keywords.

#Functions
#---
#1. api_request: function to make API request to Alliance Genome
#2. query_keywords: function to return a list of genes associated to keywords
#3. convert_df: function to convert data frame encoding for saving as CSV

#import packages
import requests
import json
import pandas as pd


#1. Function to make API request to Alliance Genome
#---
def api_request(api_url):
    try:
        #send API request
        response = requests.get(api_url)
        
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


def query_limit(keywords):
    #define an empty list
    limits = {}
    for keyword in keywords:
        #define the API url
        api_url = 'https://www.alliancegenome.org/api/search?category=gene&limit=1&q=' + keyword
        
        #obtain the API request data
        data = api_request(api_url)
        
        #define the limit for the current keyword
        limits[keyword] = data['total']
    
    #return the number of genes from each query
    return limits


#2. Function to return a list of genes associated to keywords
#---
def query_keywords(keywords, limits, species = 'Homo sapiens'):
    #define a list of column names
    columns = ['symbol','name','id','soTermName']
    
    #define an empty data frame
    gene_data = pd.DataFrame(columns = ['hgnc','name','id','biotype','keyword'])
    
    for keyword in keywords:
        #define the limit for the current keyword
        limit = limits[keyword]
        
        #create an array of offsets
        start = 0
        increment = 300
        offsets = [start]
        while start + increment < limit:
            start += increment
            offsets.append(start)
        
        offsets.append(limit - 1)
        
        #API query for each page
        for offset in offsets:
            #re-define the API url
            api_url = 'https://www.alliancegenome.org/api/search?category=gene&limit=' + str(increment) + '&offset=' + str(offset) + '&q=' + keyword
            
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
            
            #store the results
            gene_data = pd.concat([gene_data, df], ignore_index = True)
    
    #return the gene data
    return gene_data


#3. Function to convert data frame encoding for saving as CSV
#---
def convert_df(df):
    #convert data frame encoding
    return df.to_csv().encode("utf-8")
