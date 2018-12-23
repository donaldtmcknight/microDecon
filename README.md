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

library(devtools) #Load devtools

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
example <- cbind.data.frame(c("OTU1","OTU2","OTU3","OTU4","OTU5"),
                        c(0,200,1000,50,0),
                        c(100,0,1300,10,0),
                        c(100,50,1000,2,500),
                        c(120,60,1200,4,400),
                        c(60,20,1400,3,600),
                        c("k_Fungi","k__Fungi","k__Fungi; p__Ascomycota","p__Basidiomycota","k__Fungi"))
colnames(example) <- c("OTU_ID","Blank1","Blank2","Pop1_Sample1","Pop1_Sample2","Pop2_Sample3","Taxa")
```

The `decon()` function is the primary function of the package. To decontaminate these samples on default settings, use the following code 
(`data` is your data frame, `numb.blanks` specifies the number of blanks in your data frame, `taxa` specifies whether your data frame
includes a taxa column).
```
decontaminated <- decon(data = example,numb.blanks=2,taxa=T)
```

In cases with many samples and OTUs that were entirely from contamination in all samples, you may find that microDecon correctly eliminates
the contamination in most samples, but in a few samples, a handful of reads remain. This can be corrected by applying filtering thresholds
using the `remove.thesh()` function. This takes the output of `decon()`. It is best to apply the filtering within species, population, or
some other sensible a priori grouping criteria (as with `decon()` only OTUs that amplified in the blank are affected). To run it on the 
default settings use the following code (`data` is your output from `decon()`, `taxa` specifies whether your data frame includes a taxa 
column, `numb.ind` gives the number of individuals in each group [in the order that they occur in your data frame]; blanks will have been
condensed to a single mean column by `decon()` so you do not need to specify them here).
```
decontaminated.final <- remove.thresh(data = decontaminated, taxa=T, numb.ind = c(2,1))
```

For details on additional options and when you should use settings other than the defaults, see the User's Guide and help files in R.
