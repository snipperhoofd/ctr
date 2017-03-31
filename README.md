# Congruent Transcrtiptional Response (CTR)

![alt text](https://github.com/Oyserman/ctr/blob/master/CTR_workflow.png "Workflow") 
 

# Introduction:

Micrbiobial Ecologists often use genomic content to infer the biogeochemical processes and interactions within a microbial community. When organisms overlap in their genomic content, they are often inferred to possess similar traits. Conversely, divergent genomic content is indicative of niche differentiation.

Transcriptomics data may aid in making these ecological inferences by providing information on the regulation of this genetic content. Congruent transcriptional responses (CTRs) of particular metabolic modules in disparate organisms within a community may be indicative of common metabolic features shared by many organisms. Conversely, when divergent transcriptional responses are observed across numerous lineages, this is indicative of the disparate niches in the community.

In this package we present a statistical framework to identify CTRs of (pre-defined) metabolic modules in microbial communties. In this manner, genomic bins may be clustered into sub-networks that share responses for particular modules, thereby facilitating inferences on functional redundancy (e.g. intra-cluster), niche partitioning (inter-cluster), and the emergence of complex traits (integrating infromation about multiple modules).

To demonstrate the utility of this approach, we identify CTRs in polymer storage (polyphosphate, glycogen, and polyhydroxyalkanoates) and amino acid biosynthesis modules of 38 genome bins recovered from a bioreactor operated under conditions that select for organisms capable of storing polymers.

# Installation:

The software is provided as an R package. At this moment the package is still under development and not hosted on CRAN, however a preview can already be installed using the following commands in R:

```r
library(devtools)
devtools::install_github("oyserman/ctr")
```
## Dependencies 

The following packages are required and will be automatically installed when using devtools.

- plyr
- reshape
- igraph
- knitr
- hexbin
- RColorBrewer
