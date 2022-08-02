# RStatistics_ThesisProject

## GO Analysis on Duplicated Genes 

### grab_geneids_fromduplicationslist.R

#### Description:

This may require some edits depending on input files. I used this script with OrthoFinder's Duplications.tsv files to find duplications
on my node's of interest. For certain groups, I parsed out only the species ID's that I wanted.<br>

The duplications gotten from this list were used for the duplication_groups_analysis scripts below.<br>

### BPduplication_groups_analysis.R
### CCduplication_groups_analysis.R
### MFduplications_groups_analysis

#### Description:
These scripts all run this same but they consider different GO categories (for my own organizational purposes). <br>
They utilize the topGO package for a functional enrichment analysis. From there, I used a Fisher's Exact test to consider significant counts of GO terms
in my list of duplicated genes against a list of all background genes.<br>

I used a compareCluster using Enricher from the YuLab (https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html) to visualize
the most enriched terms between all my groups of interest<br>

#### REQUIRED FILES:
<ul>
  <li>GO Mapping file for all genes.</li>
  <ul>
    <li>More detail can be found in the ReadMappings() instrictuions in the topGO manual. (https://bioconductor.org/packages/devel/bioc/vignettes/topGO/inst/doc/topGO.pdf) </li>
    <li>I retrieved a mapping file using HMMER2GO (Staton 2018). File has default extension 'orfs_Pfam-A_GO_GOterm_mapping.tsv' </li>
    </ul>
  <li>Full List of Genes (Background). List seperated by new-line characters. </li>
  <li>List of Genes of Interest</li>
  <ul>
    <li> I had multiple interesting lists to compare and denoted them with group numbers.</li> 
    </ul>
  <li>Mapping between Gene ID TERMS and GO terms</li>
  <ul>
    <li> I used columns from HMMER2GO files with default extension 'orfs_Pfam-A_GO.tsv'. The first column is the PFAM term and the second column is the GO id. </li> 
    </ul>
</ul>


## GO Analysis for Evolution Models

### go_correlations_with_omega.R

#### Description:

This script utilizes the topGO package in R to do a functional enrichment analysis of dN/dS (omega) values from branch-model results returned by PAML's codeML feature. <br>
It starts by calculating with omega values fit significantly better under the branch-model than the null. <br>

Next, you can explore the distribution of omega values to find a cut-off values. Or, you can utilize the quantile feature in R to define a cut-off value. The cut-off 
values if provided to a topDiff genes function (more instructions can be found in the topGO manual linked above).

After running the the topGO object function, I used the KSelim statistical test to find enrichment based on scores (the omega values). I visualized the results
using enricher's BarPlot function. 

#### REQUIRED FILES:
<ul>
  <li> Tab-delimited output from codeML's branch model. Need: Alternative and Null LNL values. Alternative and Null dN/dS values. </li>
  <li> A GO mapping between genes (or, in my case, gene familes) and GO terms </li>
    <ul>
      <li> I created a GO map using HMMER2GO. Then, I manually parsed that files so the ID names would match my branch-model output. </li>
    </ul>
</ul>

### branchsite_go_enrichment.R

#### Description: 

This script function similarly to the above the script. It first finds which results from PAML's codeML branch-site model fit significantly better under a branch-site 
model M2a than the null models. Then, it can determine what cut-off value for positive site results. You can also just evaluate all significant genes with positively
selected sites.<br>

This script uses the topGO object and then applies the KS elim statistical test to find enrichment based on how many positively selection sites were found. I visualized the results
using enricher's BarPlot function.

#### REQUIRED FILES:

<ul>
  <li> Tab-delimited output from codeML's branch-site model. Need: Alternative and Null LNL values. </li>
  <li> List of gene families and the  number of sites under positive selection (parsed from branch-site model output) </li>
  <li> A GO mapping between genes (or, in my case, gene familes) and GO terms </li>
    <ul>
      <li> I created a GO map using HMMER2GO. Then, I manually parsed that files so the ID names would match my branch-site model output. </li>
    </ul>
</ul>

## Evolution Model Statistics

### single_copy_dnds_dist.R

#### Description: 

This script can make visualization for the distribution of omega values for background lineages and foreground lineages returned by PAML's codeML branch-model.
First, it will calculate with genes fit significantly better under the alternative model than the null. <br>

From there, it parse and plot distriubtions for the background and foreground values under the alternative model. The Komolgorov-Smirnov test is run to determine
if the distributions are significantly different.

#### REQUIRED FILES:

Tab-delmited output file from codeML's branch-model.


### Notes:

<ul>
  <li> Scripts for creating output files and parsing results from PAML and Orthofinder can be found in my EvolutionModels repository. (https://github.com/leanngward/EvolutionModels) </li>
  <li> Required installations for these scripts are usually commented out. </li>




