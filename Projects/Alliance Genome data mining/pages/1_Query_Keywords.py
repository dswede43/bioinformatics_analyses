#Query keywords page UI - Alliance Genome data mining app
#---
#Script to run the query keywords page via Streamlit.

#import packages
from PIL import Image
import streamlit as st
import json
import pandas as pd
import helpers

#define global variables
SPECIES = ['Homo sapiens',
           'Caenorhabditis elegans',
           'Danio rerio',
           'Drosophila melanogaster',
           'Mus musculus',
           'Rattus norvegicus',
           'Saccharomyces cerevisiae',
           'Xenopus laevis',
           'Xenopus tropicalis']


#Query keywords page UI
#---
def main():
    #set page customizations
    im = Image.open("static/favicon-16x16.png")
    st.set_page_config(
        page_title = "Query Keywords",
        page_icon = im)
    
    #set the web app title
    st.header('Query Keywords', divider = 'blue')
    
    #add text to explain the app
    st.write("This page is used to query the [Alliance Genome of Resources consortium](https://www.alliancegenome.org/) for biologically relevant keywords and returns lists of genes with functions associated to those keywords.")
    st.write("Please input a comma-separated list of keywords below!")
    st.write("**NOTE: the more specific the keyword is, the less genes that are associated and results in shorter query times.**")
    
    #create a form for the query input
    with st.form(key = 'my_form', clear_on_submit = False):
        species = st.selectbox('Species', SPECIES)
        keywords = st.text_input('Keywords', placeholder = 'Example: virus')
        submit_button = st.form_submit_button(label = 'Search')
    
    if submit_button or keywords:
        #format the inputted keywords
        keywords = [keyword.strip() for keyword in keywords.split(',')]
        
        #define an empty data frame
        gene_data = pd.DataFrame(columns = ['hgnc','name','id','biotype','keyword'])
        
        for keyword in keywords:
            #obtain the query limits and offsets for each keyword
            limit, offsets = helpers.query_limit(keyword)
            
            #create progress bar header
            st.write('Keyword: ' + keyword)
            
            #intialize the progress bar
            progress_bar = st.progress(0)
            total_progress = len(offsets)
            i = 0
            for offset in offsets:
                #query Alliance Genome for keywords
                df = helpers.query_keyword(keyword, limit, offset, species = species)
                
                #store the results
                gene_data = pd.concat([gene_data, df], ignore_index = True)
                
                #update progress bar
                i += 1
                progress_bar.progress(i / total_progress)
        
        #convert data for export
        csv = helpers.convert_df(gene_data)
        
        #display the download button
        helpers.show_download_button(csv)
        
        #display data
        st.dataframe(gene_data)

if __name__ == "__main__":
    main()
