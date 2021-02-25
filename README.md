# SPEDE-SAMPLER

![alt text](https://github.com/CJMvS/spede-sampler-py/blob/main/spede_sampler.png?raw=true)

SPEDE-SAMPLER is a GUI program written in Python and R that offers the following functionality:

1. Reads in an aligned Fasta file, and randomly select *n* sequences *p* times. These files are saved in an output folder.
For example, a sequence alignment of ten COI sequences might be uploaded. The user may wish to randomly select 50% of this dataset without replacement, and repeat this process 15 times. The program will thus write out 15 Fasta files to a folder, where each file contains a random selection of five sequences.

2. Loops through the output folder to produce a Maximum Likelihood (ML) phylogeny for each resampled Fasta file. There are two ML programs available: FastTree and RAxML.

3. Loops through each ML phylogeny and run a GMYC species delimitation analysis on it in order to estimate the number of clusters and entities present.

Functions (1) and (2) are carried out in the main SPEDE-SAMPLER interface. 
Function (3) has been written as an R Shiny application, and can be run from R using only the following few lines of code:

`install.packages("shiny") # install the R Shiny package`

`library(shiny) # access the library`

`shiny::runGitHub("spede-sampler", "CJMvS", ref="main") # access the GUI app housed on GitHub`

### :white_check_mark: **Summary overview of the analysis pipeline:**

(1) Resample .fas files --> (2) Create a Maximum Likelihood phylogeny for each --> (3) Run a GMYC analysis on each tree --> (4) Compare species delimitation across trees

### DOWNLOAD
Download the content of this repository, unzip the folder, and open the **spede_sampler.exe** file to run the program.

:warning: Temporarily disable your antivirus before running this executable file, as it may incorrectly identify it as malware.
