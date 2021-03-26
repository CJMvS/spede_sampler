![alt text](https://github.com/CJMvS/spede-sampler-py/blob/main/spede_sampler_logo.png?raw=true)

*Clarke van Steenderen* <br />
*Department of Zoology and Entomology* <br />
[*The Centre for Biological Control*](https://www.ru.ac.za/centreforbiologicalcontrol/) <br />
*Rhodes University, Makhanda (Grahamstown), Eastern Cape, South Africa* <br />
*2021* <br />
*e-mail:* vsteenderen@gmail.com

---
### :pencil2: FUNCTIONALITY
---

SPEDE-SAMPLER is a GUI program written in Python and R for assessing sampling effects on species delimitation using the GMCY method. The program offers the following:

1. Reads in an aligned Fasta file, and randomly select *n* sequences *p* times. These files are saved in an output folder.
For example, a sequence alignment of ten COI sequences might be uploaded. The user may wish to randomly select 50% of this dataset without replacement, and repeat this process 15 times. The program will thus write out 15 Fasta files to a folder, where each file contains a random selection of five sequences.

2. Loops through the output folder to produce a Maximum Likelihood (ML) phylogeny for each resampled Fasta file. The two ML programs available are FastTree and RAxML.

3. Loops through each ML phylogeny and runs a GMYC species delimitation analysis on it in order to estimate the number of clusters and entities present.

Functions (1) and (2) are carried out in the main SPEDE-SAMPLER interface, and function (3) has been written as an R Shiny application.

### :white_check_mark: **Summary overview of the analysis pipeline:**

(1) Resample .fas files --> (2) Create a Maximum Likelihood phylogeny for each --> (3) Run a GMYC analysis on each tree --> (4) Compare species delimitation across trees

---
### :floppy_disk: DOWNLOAD
---

Download the content of this repository, unzip the folder, and open the **spede_sampler.py** file in Python to run the program.

To run the executable file version of the program, download the **spede_sampler.exe** file [here](https://drive.google.com/drive/folders/1wrN7-eYNBcLC5ANuu5DmVPbYSElppJtc?usp=sharing) and move it to the unzipped folder (or download the entire Google Drive folder and just open the spede_sampler.exe file from there). Open this executable file to run the program. 

:heavy_exclamation_mark: Temporarily disable your antivirus before running this executable file, as it may incorrectly identify it as malware.

---
### (1) RESAMPLE FASTA FILES
---

Before starting, please ensure that your sequence names do not exceed 10 characters. RAxML requires files in Phylip format, and names will automatically be trimmed. This can result in identical sample names, which will cause errors.

Open an aligned FASTA file (.fas or .fasta extension), and run the analysis by clicking on the 'Run Resampling' option in the menu. Input the number of sequences that you would like to subset from your data, and the number of iterations that you would like to run. Please ensure that your sequence lengths are equal. 

For example, if your input .fas file contained the following three sequences:

>seq1<br>
AAAGGGTTTAA<br>
>seq2<br>
CCCGGGAAAGT<br>
>seq3<br>
GGTTTAAGGGT<br>

If you wish to randomly select two of these sequences, and iterate the process four times, the program will create four .fas files, which could be for example:

**iteration1.fas:**

>seq2<br>
CCCGGGAAAGT<br>
>sequence1<br>
AAAGGGTTTAA<br>

**iteration2.fas**

>seq2<br>
CCCGGGAAAGT<br>
>seq3<br>
GGTTTAAGGGT<br>

**iteration3.fas**

>seq3<br>
GGTTTAAGGGT<br>
>sequence1<br>
AAAGGGTTTAA<br>

**iteration4.fas**

>seq3<br>
GGTTTAAGGGT<br>
>seq2<br>
CCCGGGAAAGT<br>

**NOTE**: 
The program automatically creates a folder called 'Iterations' into which these iterations are written. If you run the program again without changing this folder name, it will be overwritten.

If you wish to append outgroups to your iteration files, upload an outgroup .fas file. In the 'Run Resampling' window, select the 'Append Outgroup Sequences?' checkbox before running the process.

The standalone code for a console-driven Python script can be found here:
https://gist.github.com/CJMvS/badf0c340c862f46db55a7b3c240a0f0

---
### (2) RUN RAxML
---

The RAxML program requires files in Phylip format (.phy). Convert these using an online tool, or use the option in SPEDE-SAMPLER to carry out this task. Once this has been done, select the folder where the resampled files (now .phy) are saved. Input the evolutionary model of your choice, and the number of bootstrap iterations you wish to run for each tree.

The standalone code for a console-driven Python script can be found here:
https://gist.github.com/CJMvS/de80fb2b4a93ff10365131f11e2b6122

---
### (2) RUN FASTTREE
---

The FastTree program requires files in Fasta format (.fas or .fasta). Select the folder where the resampled files are saved. Input the number of bootstrap iterations you wish to run for each tree.

The standalone code for a console-driven Python script can be found here:
https://gist.github.com/CJMvS/2790ceb4e3999ac6cf6eba931336901c

---
### (3) RUN GMYC ANALYSES ON RESAMPLED TREE FILES
---

GMYC analyses need to be run in R, and so a separate R Shiny application needs to be launched. To do this, open R, and run the following three lines of code:

``` r
install.packages("shiny") # install the R Shiny package
library(shiny) # access the library
shiny::runGitHub("spede-sampler", "CJMvS", ref="main") # access the GUI app housed on GitHub
```
See the GitHub repository for the application [here](https://github.com/CJMvS/spede-sampler), with a detailed README and worked example.

---
### 💡 QUICK TIP: DOWNLOAD MULTIPLE GENBANK SEQUENCES FROM R
---

As an example, let's say that you want to download [100 caddisfly COI sequences](https://bioone.org/journals/freshwater-science/volume-26/issue-4/06-089.1/Associating-larvae-and-adults-of-Chinese-Hydropsychidae-caddisflies-Insecta-Trichoptera/10.1899/06-089.1.full) with Genbank IDs EF513774 to EF513874.

In R, follow these steps:

``` r
library(ape)
setwd("path/of/your/choice")

start_caddis = 513774
end_caddis = 513874
caddis_ids = c()

k = 0
for (i in start_caddis:end_caddis){
  caddis_ids[k] = paste("EF", i, sep = "")
  k = k + 1
 } 
 
caddis_seqs = ape::read.GenBank(caddis_ids)
write.dna(caddis_seqs, file ="caddis.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = "", colw = 10)

```
---
### CONTACT AND CITATION
--- 

For suggestions, or bug-reporting, please contact me at vsteenderen@gmail.com.

Please cite as:

van Steenderen, CJM. SPEDE-SAMPLER: assess sampling effects on species delimitation, version 1.1, 2020.
