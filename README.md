The entire analysis pipeline for this project is written in the combined_motif_analysis.py, which inputs a list of proteins and UniProt IDs to run BLAST, multiple sequence alignment, and MEME Suite FIMO to generate candidate motif matches. These candidate motifs are then used to generate PyMOL scripts to visualize each motif for putative interaction with the scorpion toxin LqhIII. 

To run this code, first clone the repository into your desired directory using the following command in the terminal:
```
git clone https://github.com/annalundeen/CompBioProject
```
Then navigate into the project directory: 
```
cd CompBioProject/
```
This code requires a virtual environment. Documentation on how to set up a virtual environment in VS Code can be found here: https://code.visualstudio.com/docs/python/environments. Once the virtual environment is set up in your project folder, make sure that it is activated using the command:
```
source .venv/bin/activate
```
When the virtual environment is active, the command line should show: 
```
(.venv) user/path_to_file
```
Then, you will need to install the packages necessary to run the code within the virtual environment:
```
.venv/bin/python -m pip install numpy pandas requests biopython pymemesuite
```
Once these packages are installed, you can run the comprehensive python script: 
```
.venv/bin/python combined_motif_analysis.py
```
