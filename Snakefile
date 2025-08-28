dir = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/"
nchunks = 500
cor_threshold = 0.1
BIOMARKERS = ["w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final", "abeta_ratio"]
rule all:
    input:
        expand(
            dir + "results_model{model}/results_model{model}_{biomarker}_{chunk}.RDS",
            model=range(1, 7),
            biomarker=BIOMARKERS,
            chunk=range(1, nchunks+1)
        )
        
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
  output: merged_data = expand(dir + "data/meth_pheno_data.{chunk}.RDS", chunk = range(1, nchunks+1)), 
          cpg_list = dir + "cpg_list.RDS", 
          manifest = dir + "manifest.RDS" 
  params: nchunks = nchunks
  script: "1-Data_preprocessing.R"

# Run null model
rule null_model:
  input: data = dir + "data/meth_pheno_data.{chunk}.RDS"
  output: out = dir + "null_model/null_residuals.{chunk}.RDS"
  script: "2.0-null_model.R"
  
# Combine residuals
rule combine_res:
  input: res_chunk = expand(dir + "null_model/null_residuals.{chunk}.RDS", chunk = range(1, nchunks+1)),
         manifest = dir + "manifest.RDS"
  output: residuals = expand(dir + "null_model/null_residuals.chr{chr}.RDS", chr = range(1,23))
  params: nchunks = nchunks
  script: "2.1-combine_res.R"
  
# Filter correlation residuals
rule cor_filter:
  input: residuals = dir + "null_model/null_residuals.chr{chr}.RDS"
  output: residuals_filter = dir + "residuals_chr/keep_cpg.{chr}.RDS"
  params: cor_threshold = cor_threshold
  script: "3.0-cor_filter.R"

# PCA + plotting + table
rule PCA:
  input: residuals_filter = expand(dir + "residuals_chr/keep_cpg.{chr}.RDS", chr = range(1, 23)),
         phenos = dir + "data/meth_pheno_data.1.RDS",
         manifest = dir + "manifest.RDS"
  output: pca_result = dir + "pca_result.RDS",
          scree_plot = dir + "plots/scree_plot.png",
          cor_table = dir +"tables/cor_table_PC.txt",
          pca_loading = dir + "plots/pca_loading.pdf"
  script: "3.5-PCA.R"

rule run_model:
    input:
        data = dir + "data/meth_pheno_data.{chunk}.RDS",
        manifest = dir + "manifest.RDS",
        pca_result = dir + "pca_result.RDS"
    output:
        results = dir + "results_model{model}/results_model{model}_{biomarker}_{chunk}.RDS"
    params:
        biomarker = lambda wc: wc.biomarker,
        model = lambda wc: wc.model
    script:
        lambda wildcards: f"4.{wildcards.model}-EWAS_model{wildcards.model}.R"
