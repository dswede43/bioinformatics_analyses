#App UI - Alliance Genome data mining app
#---
#Script to run the app UI via Streamlit.

#import packages
import streamlit as st
import json
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


#App UI
#---
def main():
    #set the web app title
    st.header('Alliance Genome data mining', divider = 'blue')
    
    #add text to explain the app
    st.write("This web application is used to query the [Alliance Genome of Resources consortium](https://www.alliancegenome.org/) for biologically relevant keywords and returns lists of genes with functions associated to those keywords. Data is obtained through API queries with [Alliance Genome's API](https://www.alliancegenome.org/swagger-ui/). Source code can be found on [Github](https://github.com/dswede43).")
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
        
        #obtain the query limits for each keyword
        limits = helpers.query_limit(keywords)
                
        #query Alliance Genome for keywords
        gene_data = helpers.query_keywords(keywords, limits, species = species)
        
        #add download data button        
        csv = helpers.convert_df(gene_data)
        st.download_button(
            label = "Download data as CSV",
            data = csv,
            file_name = "keyword_gene_associations.csv",
            mime = "text/csv")
        
        #display data
        st.dataframe(gene_data)

if __name__ == "__main__":
    main()
