# microDecon
An R package for removing contamination from metabarcoding (e.g., microbiome) datasets post-sequencing

## Description
Contamination is a serious and widespread problem in metabarcoding studies (particularly bacterial microbiome studies). microDecon provides 
an easy method of identifying and removing contaminant reads post-sequencing. It uses the data in one or more blank samples that were
carried throughout the entire project (collection, extraction, amplification, and sequencing) to identify and remove the reads in each
actual sample that are from contamination. Because it removes reads rather than entire OTUs, it is a large improvement over existing
methods. It was designed specifically for microbiome sequencing studies (e.g., 16S bacterial studies and ITS fungal studies), but it
should be applicatble to metabarcoding studies more generally. A more detailed explanation of the algorithms it uses can be found in the
User's Guide.

## Installation
Simply copy, paste, and run the following code in R.
```
install.packages("devtools") #Install devtools (if not already installed)
devtools::install_github("donaldtmcknight/microDecon") #Install microDecon
```

## Running the package
Load microDecon
```
library(microDecon)
```

Load a data frame with a column of OTU IDs, at least one column from a blank sample, at least one column from actual samples, and (optional)
a column with taxonomic information. For example, the following OTU data frame has two blanks, three samples (from two populations), and a taxa column (the columns must be in that order: OTU IDs, blanks, samples [samples should be grouped by population, species, etc.], taxa).
```
example <- cbind.data.frame(c("OTU1","OTU2","OTU3","OTU4","OTU5","OTU6"),
                        c(0,200,1000,50,0,25),
                        c(0,220,800,30,0,10),
                        c(0,180,1300,70,0,30),
                        c(60,660,1440,70,2400,30),
                        c(64,520,1000,48,1900,20),
                        c(40,480,700,35,2100,15),
                        c("K_Bacteria; P_Actinobacteria","K_Bacteria; P_Proteobacteria","K_Bacteria; P_Proteobacteria","K_Bacteria; P_Bacteroidetes","K_Bacteria","K_Bacteria"))
colnames(example) <- c("OTU_ID","Blank1","Blank2","Blank3","Pop1_Sample1","Pop1_Sample2","Pop2_Sample3","Taxa")
```

The `decon()` function is the primary function of the package (it is a warpper for all other functions). To decontaminate these samples on default settings, use the following code 
(`data` is your data frame, `numb.blanks` specifies the number of blanks in your data frame, `numb.ind` specifies the number of individuals per group (e.g. population), `taxa` specifies whether your data frame
includes a taxa column).
```
decontaminated <- decon(data = example,numb.blanks=3,numb.ind=c(2,1),taxa=T)
```
For details on additional options and when you should use settings other than the defaults, see the User's Guide and help files in R.
