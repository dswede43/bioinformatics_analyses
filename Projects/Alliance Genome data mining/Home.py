#Homepage UI - Alliance Genome data mining app
#---
#Script to run the app homepage via Streamlit.

#import packages
from PIL import Image
import streamlit as st

#define global variables
FUNCTIONS = {'**Query Keywords**': 'Input a list of biologically relevant keywords and return lists of genes with functions associated to those keywords',
             '**Summarize Genes**': 'Input a list of genes and use Natural Language Processing (NLP) to summarize the biological functions of those genes'}


#Homepage UI
#---
def main():
    #set page customizations
    im = Image.open("static/favicon-16x16.png")
    st.set_page_config(
        page_title = "Genome Data Mining",
        page_icon = im)
    
    #set the web app title
    st.header('Alliance Genome Data Mining', divider = 'blue')
    
    #create a sidebar
    st.sidebar.success('')
    
    #add text to explain the app
    st.write("This web application is used to obtain genome-related data through API queries with [Alliance Genome's API](https://www.alliancegenome.org/swagger-ui/.")
    
    #text to summarize the app functions
    st.write('**This application is capable of returning the following genome-related data:**')
    functions_list = ''
    for function in FUNCTIONS:
        functions_list += "- " + function + ": " + FUNCTIONS[function] + "\n"
    st.markdown(functions_list)
    st.write("Source code can be found on [Github](https://github.com/dswede43/bioinformatics_analyses/tree/main/Alliance%20Genome%20data%20mining).")

if __name__ == "__main__":
    main()
