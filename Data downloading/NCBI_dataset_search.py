#!/usr/bin/python3

#NCBI dataset searching
#---
#Script to automate NCBI searches to find GEO or SRA omics datasets.


#Table of contents:
#---
#1. Search GEO
#2. Search SRA


#import libraries
from Bio import Entrez
import re

#define global variables
SEARCH = "environmental DNA" #search terms


#Functions
#---
def search_ncbi(database, term):
    """Searches a specific NCBI database with the provided search term."""
    handle = Entrez.esearch(db = database, term = term, retmax = 10)
    results = Entrez.read(handle)
    handle.close()
    return results["IdList"]

def fetch_summary(database, id_list):
    """Fetches a summary of the records for the provided IDs."""
    handle = Entrez.esummary(db=database, id=",".join(id_list))
    summaries = Entrez.read(handle)
    handle.close()
    return summaries


#1. Search GEO
#---
#find the GEO id's
print("Searching for GEO datasets...")
geo_ids = search_ncbi("gds", SEARCH)

if len(geo_ids) > 0:
    #fetch summaries of datasets
    geo_summaries = fetch_summary("gds", geo_ids)
    print("GEO results:")
    
    #print the summaries
    for summary in geo_summaries:
        print(f"Title: {summary['title']}")
        print(f"Accession: {summary['Accession']}")
else:
    print("No datasets found!")


#2. Search SRA
#---
#find the GEO id's
print("Searching for SRA datasets...")
sra_ids = search_ncbi("sra", SEARCH)

if len(sra_ids) > 0:
    #fetch summaries of datasets
    sra_summaries = fetch_summary("sra", sra_ids)
    print("SRA results:")
    
    #print the summaries
    for summary in sra_summaries:
        title = re.search(r"<Title>(.*?)</Title>", summary['ExpXml']).group(1)
        print(title)
        acc = re.search(r'<Submitter acc="(.*?)" center_name', summary['ExpXml']).group(1)
        print(acc)
else:
    print("No datasets found!")
