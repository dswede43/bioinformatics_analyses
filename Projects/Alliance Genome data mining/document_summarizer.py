#Document summarizer functions - Alliance Genome data mining app
#---
#Functions for machine learning methods to summarize gene functions.

#Functions
#---

#import packages
from transformers import pipeline
import pandas as pd
import re
import helpers

#define global variables
HGNC_SYMBOLS = ['IFNG','IFNA1','IFNB1']

#Obtain gene functions
#---
hgnc_ids, invalid_symbols = helpers.convert_hgnc_symbols(HGNC_SYMBOLS)

gene_functions = helpers.query_gene_functions(HGNC_SYMBOLS)

gene_texts = {}
for gene_name in gene_functions:
    gene_function = gene_functions[gene_name]
    gene_function = re.sub(r'\[.*?\]', '', gene_function)
    gene_function = gene_function.split()
    if len(gene_function) > 1:
        gene_function = gene_function[1:]
    else:
        gene_function = []
    gene_function = " ".join(gene_function)
    gene_texts[gene_name] = gene_name + " " + str(gene_function) + " "

gene_texts = pd.DataFrame(gene_texts.items())
gene_texts.columns = ['gene', 'function']

#define the summarization model
summarizer = pipeline("summarization", device = 0)

def summarize_text(text, summarizer, min_length = 5, max_length = 50):
    #apply the summarization
    gene_summary = summarizer(text, #text to summarize
                              min_length = min_length, #min length of summary
                              max_length = max_length, #max length of summary
                              do_sample = False) #use greedy-decoder
    return gene_summary[0]['summary_text']

gene_texts['summary'] = gene_texts['function'].apply(lambda x: summarize_text(x, summarizer))
