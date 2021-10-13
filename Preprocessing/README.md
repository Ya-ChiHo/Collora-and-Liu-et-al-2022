20210521scRFEobjectsUnstim.R - generating scRFE objects for figure 7. 

20210602GSEAHIVinfectedcells_ranking.R - ranking gene expression in HIV-1 RNA+ cells in various condtions for figure 4 and figure 7. 

20210605scRFEprocessing.R - wrapping the replicates of scRFE into an R parsable summary. 

20210727Clonalcorrelation.R - Generating the clonecorrelation and GSEA ranking for figure 2 and 3. 

20210727elasticnetcsvoutput.R - exporting the matricies for elastic net regression notebook (figure 2). 

20210730ModuleID.R - script is for identifying modules, including several R functions written by Sam Kazer (and appear in his 2019 Nature Med paper), and wrappers for some of those functions written by Jack Collora. 

20210803GenerateCorMatrixForMods.R - generates pearson correlation matrix for module plotting. 

2021Mastergeneration.R - script generates the "master" seurat objects with binarized ADT data, bulk clone quantification, and whatnot. 

210126stimintegrationstep1.R - preps the libraries from stim data for integration. 

Firstintegrationscript_unstim.R - integrates the libraries from unstim data. 

JackUtilities.R - series of functions Jack wrote to do various plottings, GSEA steps, etc. Each function is commented. 

Secondintegrationscript_unstim.R - reintigrates unstimed data after contaminating (non-CD4+ T cells and non-T cells) are removed. 

stimintegrationstep2.R - integrates libraries from the stim data. 

