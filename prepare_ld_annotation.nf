nextflow.enable.dsl=2

process make_bed_file_from_vcf_add_cm_info {
    publishDir "${params.intermediates}/${id}", mode: 'rellink', overwrite: true
    input:
      tuple val(id), val(chr), path(vcfin)
    output:
      tuple val(id), path('extract_positions_from_vcf_and_create_index')
    script:
      """
      make_bed_file_from_vcf_add_cm_info.sh ${vcfin} ${id}
      """
}

workflow prepare_1000G {
    take: vcf_file
    main:
    foo(vcf_file)
    bar(foo.out)
}

workflow {
   // read in metafile
  Channel
   .fromPath("${params.input}")
   .splitCsv( header:false, sep:" " )
   .map { it1, it2 -> tuple(file(it1).getBaseName(),file(it1),file("${it1}.tbi"),file(it2)) }
   .set { vcf_filename_tracker_added }
  
  // Make channel
  channel.fromPath("${params.input}")
    .map { file -> tuple(file.baseName - ~/\.vcf.gz/, file) }
    .transpose()
    .set { vcf_filename_tracker_added }

  // Start processing
  extract_positions_from_vcf_and_create_index(vcf_filename_tracker_added)

}

