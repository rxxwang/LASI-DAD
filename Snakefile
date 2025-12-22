dir = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/"
nchunks = 500
cor_threshold = 0.1
BIOMARKERS = ["w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final", "abeta_ratio"]
#BIOMARKERS = ["w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final", "abeta_ratio", "w1abeta42_final", "w1abeta40_final"]
variables = ["biomarker", "age", "gender", "smoke"]

rule all:
  input:
    dir + "gfadata.Rda"
    

biomarker_data_file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/Alzheimers_Disease/Updated LASI-DAD AD biomarker data/lasidad_w12adbio_final_methylID.dta"
manifest_file = '/net/orion/skardia_lab/clubhouse/research/projects/LASI/Methylation_Nov14_2024/Scott_QC_test/Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024.rda'
gap_cpg_file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250528_EWAS1/HRS Gap Probes v1.xlsx"
meth_data_file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/Methylation_Nov14_2024/Scott_QC_test/Meth_signal_noob_928074x2290_beadlt4NA.RData"
unmeth_data_file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/Methylation_Nov14_2024/Scott_QC_test/Unmeth_signal_noob_928074x2290_beadlt4NA.RData"
methylation_data_file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/Methylation.Rda"
phenotype_data_file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/Methylation_Nov14_2024/Scott_QC_test/LASI_phenos_2340.RData"
blood_cell_data_file = '/net/orion/skardia_lab/clubhouse/research/projects/LASI/Methylation_Nov14_2024/Scott_QC_test/LASI_WBC.csv'

# Data prepocessing
# rule process_data:
#  input: 
#    biomarker = biomarker_data_file, 
#    manifest = manifest_file, 
#    gap = gap_cpg_file, 
#    # meth = methylation_data_file, 
#    meth = meth_data_file,
#    unmeth = unmeth_data_file,
#    pheno = phenotype_data_file, 
#    wbc = blood_cell_data_file
#  output: 
#    merged_data = expand(dir + "data4/meth_pheno_data.{chunk}.RDS", chunk = range(1, nchunks+1)), 
#    cpg_list = dir + "cpg_list4.RDS", 
#    manifest = dir + "manifest4.RDS" 
#  params: 
#    nchunks = nchunks
#  resources: 
#    mem_mb = 100000,
#    runtime = 1200
#  script: 
#    "1-Data_preprocessing.R"

# Run null model
#rule null_model:
#  input: 
#    data = dir + "data4/meth_pheno_data.{chunk}.RDS"
#  output: 
#    out = dir + "null_model4/null_residuals.{chunk}.RDS"
#  resources: 
#    mem_mb = 10000,
#    runtime = 2400 
#  script: 
#    "2.0-null_model.R"
  
# Combine residuals
#rule combine_res:
#  input: 
#    res_chunk = expand(dir + "null_model4/null_residuals.{chunk}.RDS", chunk = range(1, nchunks+1)),
#    manifest = dir + "manifest4.RDS"
#  output: 
#    residuals = expand(dir + "null_model4/null_residuals.chr{chr}.RDS", chr = range(1,23))
#  params: 
#    nchunks = nchunks
#  resources: 
#    mem_mb = 50000,
#    runtime = 300
#  script: 
#    "2.1-combine_res.R"
  
# Filter correlation residuals
#rule cor_filter:
#  input: residuals = dir + "null_model4/null_residuals.chr{chr}.RDS"
#  output: residuals_filter = dir + "residuals_chr4/keep_cpg.{chr}.RDS"
#  params: cor_threshold = cor_threshold
#  resources: 
#    mem_mb = 7000,
#    runtime = 300
#  script: "3.0-cor_filter.R"

# PCA + plotting + table
#rule PCA:
#  input: residuals_filter = expand(dir + "residuals_chr4/keep_cpg.{chr}.RDS", chr = range(1, 23)),
#         phenos = dir + "data4/meth_pheno_data.1.RDS",
#         manifest = dir + "manifest4.RDS"
#  output: 
#    pca_result = dir + "pca_result4.RDS",
#    scree_plot = dir + "plots/scree_plot4.png",
#    cor_table = dir +"tables/cor_table_PC4.txt",
#    pca_loading = dir + "plots/pca_loading4.pdf"
#  resources: 
#    mem_mb = 7000,
#    runtime = 300
#  script: "3.5-PCA.R"

# rule run_model:
#   input:
#     data = dir + "data3/meth_pheno_data.{chunk}.RDS",
#     manifest = dir + "manifest3.RDS",
#     pca_result = dir + "pca_result3.RDS"
#   output:
#     results = dir + "results3_model{model}/results.model{model}.{biomarker}.{chunk}.RDS"
#   params:
#     biomarker = lambda wc: wc.biomarker,
#     model = lambda wc: wc.model
#   resources: 
#     mem_mb = 10000,   # 10 GB
#     runtime = 2400    # 40 hours
#   script:
#     "4-EWAS_model.R"
#   
# rule qqplot:
#   input: 
#     data = expand(dir + "results3_model{{model}}/results.model{{model}}.{{biomarker}}.{chunk}.RDS", chunk = range(1, nchunks+1))  
#   output: 
#     combined_result = dir + "results3/model{model}.{biomarker}.RDS",
#     qqplot = expand(dir + "qqplot3/qqplot_model{{model}}.{{biomarker}}.{variable}.png", variable = variables)
#   params:
#     nchunks = nchunks,
#     biomarker = lambda wc: wc.biomarker,
#     model = lambda wc: wc.model
#   resources: 
#     mem_mb = 20000,   # 20 GB
#     runtime = 300    # 5 hours
#   script:
#     "5.0-qqplot.R"
# 
# rule combine_qqplot:
#   input: 
#     qqplot = expand(dir + "qqplot3/qqplot_model{model}.{{biomarker}}.biomarker.png", model = range(1, 7))
#   output: 
#     combined_qqplot = dir + "combined_qqplot3/qqplot.{biomarker}.biomarker.png"
#   params:
#     biomarker = lambda wc: wc.biomarker,
#   resources: 
#     mem_mb = 5000,   # 20 GB
#     runtime = 60    # 5 hours
#   script:
#     "5.2-qqplot_combine.R"

rule run_model1_nr:
  input:
    data = dir + "data4/meth_pheno_data.{chunk}.RDS",
    manifest = dir + "manifest4.RDS",
    pca_result = dir + "pca_result4.RDS"
  output:
    results = dir + "results4_model{model}/results.model{model}.{biomarker}.{chunk}.RDS"
  params:
    biomarker = lambda wc: wc.biomarker,
    model = lambda wc: wc.model
  resources: 
    mem_mb = 10000,   # 10 GB
    runtime = 2400    # 40 hours
  script:
    "4.0-EWAS_model_nr.R"
    
rule combine_results_nr:
  input: 
    data = expand(dir + "results4_model{{model}}/results.model{{model}}.{{biomarker}}.{chunk}.RDS", chunk = range(1, nchunks+1))  
  output: 
    combined_result = dir + "results4/model{model}.{biomarker}.RDS"
  params:
    nchunks = nchunks,
    biomarker = lambda wc: wc.biomarker,
    model = lambda wc: wc.model
  resources: 
    mem_mb = 20000,   # 20 GB
    runtime = 300    # 5 hours
  script:
    "6.0-combine_results_nr.R"
    
rule fillin_nan:
  input: 
    model10_data = dir + "results4/model{model}0.{biomarker}.RDS",
    model11_data = dir + "results4/model{model}1.{biomarker}.RDS"
  output: 
    fillin_nan_result = dir + "results4_fillin_nan/model{model}.fillin.nan.{biomarker}.RDS",
    qqplot = dir + "qqplot4/qqplot_model{model}.fillin.nan.{biomarker}.png"
  params:
    biomarker = lambda wc: wc.biomarker
  resources: 
    mem_mb = 20000,   # 20 GB
    runtime = 300    # 5 hours
  script:
    "6.1-fillin_nan.R"

rule manhattan:
  input: 
    data = dir + "results4_fillin_nan/model{model}.fillin.nan.{biomarker}.RDS"
  output: 
    manhattan = dir + "plots/manhattan.model{model}.{biomarker}.png"
  params:
    biomarker = lambda wc: wc.biomarker
  resources: 
    mem_mb = 20000,   # 20 GB
    runtime = 300    # 5 hours
  script:
    "6.2-manhattan.R"

rule combine_qqplot:
  input:
    qqplot = expand( dir + "qqplot4/qqplot_model{{model}}.fillin.nan.{biomarker}.png", biomarker = BIOMARKERS)
  output:
    combined_qqplot = dir + "qqplot4/qqplot_model{model}.fillin.nan.png"
  resources:
    mem_mb = 20000,   # 20 GB
    runtime = 300    # 5 hours
  script:
    "6.3-qqplot_combine.R"

rule heatmap_upset:
  input:
    data = expand(dir + "results4_fillin_nan/model{{model}}.fillin.nan.{biomarker}.RDS", biomarker = BIOMARKERS),
    manifest = dir + "manifest4.RDS"
  output:
    pval = dir + "pval.model{model}.RDS",
    heatmap = dir + "plots/heatmap_model{model}.pdf",
    upset_data = dir + "upset_data_model{model}.RDS"
  resources:
    mem_mb = 20000,   # 20 GB
    runtime = 300    # 5 hours
  script:
    "6.4-heatmap_upset.R"

rule cor_filter_rearrange:
  input:
    pval = dir + "pval.model2.RDS",
    residuals = dir + "null_model4/null_residuals.chr{chr}.RDS",
    manifest = dir + "manifest4.RDS"
  output:
    residuals_filter = dir + "residuals_chr4_rearrange/keep_cpg.{chr}.RDS"
  params:
    chr = lambda wc: wc.chr,
    cor_threshold = cor_threshold
  resources:
    mem_mb = 20000,   # 20 GB
    runtime = 300    # 5 hours
  script:
    "7.0-cor_filter_rearrange.R"

rule GFA:
  input:
    data = expand(dir + "results4_fillin_nan/model2.fillin.nan.{biomarker}.RDS", biomarker = BIOMARKERS),
    residuals_filter = expand(dir + "residuals_chr4_rearrange/keep_cpg.{chr}.RDS", chr = range(1, 23)),
    pheno = dir + "data4/meth_pheno_data.1.RDS"
  output: 
    gfadata = dir + "gfadata.Rda"
  resources:
    mem_mb = 20000,
    runtime = 20
  script:
    "7.1-GFA.R"