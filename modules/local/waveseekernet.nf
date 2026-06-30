process WAVESEEKERNET_ANALYZE {
  tag "$sample"
  label 'process_medium'
  maxForks 1

  container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? 'waveseekernet.sif' : 'waveseekernet:latest' }"

  input:
  tuple val(sample), path(consensus_fasta), path(cds_fasta)
  path config
  val weights_dir

  output:
  tuple val(sample), path("${prefix}.waveseekernet.tsv"),       emit: predictions
  tuple val(sample), path("${prefix}.waveseekernet.html"),      emit: html_report, optional: true
  tuple val(sample), path("**/*_shap_*.tsv"), optional: true,   emit: shap_tsv
  tuple val(sample), path("**/*_shap_*.png"), optional: true,   emit: shap_png
  path "versions.yml",                                          emit: versions

  script:
  def args = task.ext.args ?: ''
  prefix   = task.ext.prefix ?: "${sample}"
  """
  waveseekernet_analyze.py \\
    --consensus-fasta $consensus_fasta \\
    --cds-fasta $cds_fasta \\
    --config $config \\
    --weights-dir $weights_dir \\
    --output-predictions ${prefix}.waveseekernet.tsv \\
    --output-shap-prefix ${prefix} \\
    $args

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      waveseekernet: \$(waveseekernet_analyze.py --version 2>&1 | echo "1.0")
  END_VERSIONS
  """
}