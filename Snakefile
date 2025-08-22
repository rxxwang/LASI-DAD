dir = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/"

rule all:
  input: expand(dir + "meth_pheno_data.{chunk}.RDS", chunk = range(1, 101))

biomarker_data_file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/Alzheimers_Disease/Updated LASI-DAD AD biomarker data/lasidad_w12adbio_final_methylID.dta"

manifest_file = '/net/orion/skardia_lab/clubhouse/research/projects/LASI/Methylation_Nov14_2024/Scott_QC_test/Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024.rda'

gap_cpg_file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250528_EWAS1/HRS Gap Probes v1.xlsx"

methylation_data_file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250528_EWAS1/Methylation.Rda"
phenotype_data_file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/Methylation_Nov14_2024/Scott_QC_test/LASI_phenos_2340.RData"
blood_cell_data_file = '/net/orion/skardia_lab/clubhouse/research/projects/LASI/Methylation_Nov14_2024/Scott_QC_test/LASI_WBC.csv'

# Data prepocessing
rule process_data:
  input: biomarker = biomarker_data_file, 
         manifest = manifest_file, 
         gap = gap_cpg_file, 
         meth = methylation_data_file, 
         pheno = phenotype_data_file, 
         wbc = blood_cell_data_file
  output: merged_data = expand(dir + "meth_pheno_data.{chunk}.RDS", chunk = range(1, 101)), 
          cpg_list = dir + "cpg_list.RDS", 
          manifest = dir + "manifest.RDS" 
  script: "1-Data_preprocessing.R"
