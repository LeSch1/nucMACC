process danpos_mono{
  label 'big'
  container 'uschwartz/danpos'

  publishDir "${params.outDir}/RUN/05_MONO-NUCS_PROFILE", mode: 'copy', pattern: "*_monoNucs_profile.bw"

  input:
  tuple val(sampleID), file(bam)
  file(chrSizes)

  output:
  file("*_monoNucs_profile.bw")
  tuple val(sampleID), file("result/pooled/*.xls")


  script:
  """
  danpos.py dpos $bam -m 1 --extend 70 -c $params.genomeSize \
  -u 0 -z 20 -a 10 -e 1  > $sampleID"_DANPOS_stats.txt"
  wigToBigWig result/pooled/*.wig -clip $chrSizes $sampleID"_monoNucs_profile.bw"
  """
}

process danpos_sub{
  label 'big'
  container 'uschwartz/danpos'

  publishDir "${params.outDir}/RUN/06_SUB-NUCS_PROFILE", mode: 'copy', pattern: "*_subNucs_profile.bw"

  input:
  tuple val(sampleID), file(bam)
  file(chrSizes)

  output:
  file("*_subNucs_profile.bw")
  tuple val(sampleID), file("result/pooled/*.xls")

  script:
  """
  danpos.py dpos $bam -m 1 --extend 70 -c $params.genomeSize \
  -u 0 -z 70 -a 20 -e 1  > $sampleID"_DANPOS_stats.txt"
  wigToBigWig result/pooled/*.wig -clip $chrSizes $sampleID"_subNucs_profile.bw"
  """
}
