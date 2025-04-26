# ONT Pipeline
**Broeckel Lab Bioinformatics (2024)**

Our data processing and analysis pipeline for sample data pulled from the PromethION:

1. Prepare data for transfer                                                                               
      1.1. Concatinate FASTQ files                                            
      1.2. (optional) QC the FASTQs via pycoQC                                                            
2. Read Alignment via Minimap2 (FASTQs w/ GRCh38)                                        
       2.1. (optional) Alignment QC via pycoQC                                         
3. SNP/INDEL via PMDV                                                           
4. Large SV via Spectre                                                          
5. SNP Annotation via snpEff                                      
       5.1. Filter snpEff VCF (NOT AUTOMATIC. SHOULD BE DONE AFTERWARDS)                                                 
6. SV-Calling via SNIFFLES2
   

## Development Notes 📝
  - Cannot run more than two instances of clair3_run.sh - too resource intensive on MAPLE.

### To-Do 
  - Pipeline errors out if there aren't .txt files in the prom dir. This can happen if a run errors out.
  - Option to skip steps?
  - Should make a custom conda environment for this pipeline
       -> This should check for that
  - Check to configure sshd_config to not timeout during transfers (or run on mobaxterm)
  - Investigate Sniffles2 subprocesses not terminating
  - Add snpEff outside of home directory for all user access.
  

### Completed additions ✅
  - Option to skip user confirmations?
  - Check whether input ends in '/' and add if not.
  - Sequencing summary file naming? Not *_combined.txt
  - Filter_Snp_Eff.py should be taken out and just run seperately. Doesn't work in pipeline rn. 
  - Print version of each tool to a .log file
