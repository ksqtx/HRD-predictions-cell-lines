#!/usr/bin/env nextflow
// Pre-process VCF files for CHORD

if( !params.sampleInfo ) error "Missing sampleInfo parameter"

Channel
    .fromPath(params.sampleInfo)
    .splitCsv(header: ['id', 'source', 'sage_vcf', 'sage_index', 'gripss_vcf', 'gripss_index'], skip: 1, strip: true)
    .into {sageInfo; grippsInfo; grippsPanel}

chFasta = Channel.value(file(params.fasta))
chFastaIdx = Channel.value(file(params.fastaIndex))


////////////////////////////////////////////////////
/* --          Filter SNVs and Indels          -- */
////////////////////////////////////////////////////

process splitDatabaseMultiallelics {
    label 'cpus_4'

    input:
    path 'gnomad.vcf.gz' from params.gnomad
    path 'gnomad.vcf.gz.tbi' from params.gnomadIndex
    path 'dbsnp.vcf.gz' from params.dbsnp
    path 'dbsnp.vcf.gz.tbi' from params.dbsnpIndex
    path 'known_indels.vcf.gz' from params.knownIndels
    path 'known_indels.vcf.gz.tbi' from params.knownIndelsIndex

    output:
    path 'gnomad.norm.vcf.gz' into splitGnomadVcf
    path 'gnomad.norm.vcf.gz.tbi' into splitGnomadIndex
    path 'dbsnp.norm.vcf.gz' into splitDbsnpVcf
    path 'dbsnp.norm.vcf.gz.tbi' into splitDbsnpIndex
    path 'known_indels.norm.vcf.gz' into splitKnownIndelsVcf
    path 'known_indels.norm.vcf.gz.tbi' into splitKnownIndelsIndex

    script:
    """
    bcftools norm \
      --multiallelics -any \
      --threads ${task.cpus} \
      gnomad.vcf.gz | \
    bgzip -c > gnomad.norm.vcf.gz

    tabix -p vcf gnomad.norm.vcf.gz

    bcftools norm \
      --multiallelics -any \
      --threads ${task.cpus} \
      dbsnp.vcf.gz | \
    bgzip -c > dbsnp.norm.vcf.gz

    tabix -p vcf dbsnp.norm.vcf.gz

    bcftools norm \
      --multiallelics -any \
      --threads ${task.cpus} \
      known_indels.vcf.gz | \
    bgzip -c > known_indels.norm.vcf.gz

    tabix -p vcf known_indels.norm.vcf.gz
    """
}

process normalizeSageVcfs {
    input:
    tuple val(id), val(source), file('sage.vcf.gz'), file('sage.vcf.gz.tbi') from sageInfo.map{ [it.id, it.source, file(it.sage_vcf), file(it.sage_index)] }
    file('reference.fasta') from chFasta
    file('reference.fasta.fai') from chFastaIdx

    output:
    tuple file('*.sage.pave.norm.vcf.gz'), file('*.sage.pave.norm.vcf.gz.tbi') into normalizedSageVcfsPanel
    tuple val(id), val(source), file('*.sage.pave.norm.vcf.gz'), file('*.sage.pave.norm.vcf.gz.tbi') into normalizedSageVcfsFilter
    path "*_path.txt" into normalizedSagePaths

    script:
    """
    printf "%s_%s" "${id}" "${source}" > sample_name.txt

    bcftools norm \
      --fasta-ref reference.fasta \
      --atomize \
      --atom-overlaps . \
      --multiallelics -any \
      --check-ref x \
      sage.vcf.gz | \
    bcftools reheader \
      --samples sample_name.txt | \
    bgzip -c > ${id}_${source}.sage.pave.norm.vcf.gz

    tabix -p vcf ${id}_${source}.sage.pave.norm.vcf.gz

    printf "%s_%s.sage.pave.norm.vcf.gz" "${id}" "${source}" > ${id}_${source}_sage_norm_path.txt
    """
}

process createSagePON {
    label 'cpus_4'
    
    publishDir "${params.outDir}/", pattern: '*_panel.{vcf.gz,vcf.gz.tbi}', mode: 'copy'

    input:
    file(vcfs) from normalizedSageVcfsPanel.collect()
    path 'vcf_paths.txt' from normalizedSagePaths.collectFile(name: 'norm_vcf_paths.txt', newLine: true)

    output:
    path 'sage_panel.vcf.gz' into sagePanelVcf
    path 'sage_panel.vcf.gz.tbi' into sagePanelIndex

    script:
    """
    bcftools merge \
      --file-list vcf_paths.txt \
      --merge indels  \
      --threads ${task.cpus} | \
    bcftools +setGT \
      -- --target-gt q --new-gt c:0/1 -i 'FMT/DP>0' | \
    bcftools +fill-tags \
      -- --tags NS | \
    bcftools view \
      -i 'INFO/NS > 2' | \
    bgzip -c > sage_panel.vcf.gz

    tabix -p vcf sage_panel.vcf.gz
    """
}

process removeCommonVariants {
    input:
    tuple val(id), val(source), file('norm.vcf.gz'), file('norm.vcf.gz.tbi') from normalizedSageVcfsFilter
    path 'gnomad.vcf.gz' from splitGnomadVcf
    path 'gnomad.vcf.gz.tbi' from splitGnomadIndex
    path 'dbsnp.vcf.gz' from splitDbsnpVcf
    path 'dbsnp.vcf.gz.tbi' from splitDbsnpIndex
    path 'known_indels.vcf.gz' from splitKnownIndelsVcf
    path 'known_indels.vcf.gz.tbi' from splitKnownIndelsIndex
    path 'sage_panel.vcf.gz' from sagePanelVcf.first()
    path 'sage_panel.vcf.gz.tbi' from sagePanelIndex.first()

    output:
    tuple val(id), val(source), file('*.sage.pave.norm.no_common.vcf.gz'), file('*.sage.pave.norm.no_common.vcf.gz.tbi') into somaticSageVcfs

    script:
    """
    bcftools isec \
      norm.vcf.gz \
      gnomad.vcf.gz \
      dbsnp.vcf.gz \
      known_indels.vcf.gz \
      sage_panel.vcf.gz \
      --complement \
      --collapse indels \
      -w1 | \
    bgzip -c > ${id}.sage.pave.norm.no_common.vcf.gz

    tabix -p vcf ${id}.sage.pave.norm.no_common.vcf.gz
    """
}

process filterQualitySage {
    publishDir "${params.outDir}/${id}_${source}/", pattern: '*.{vcf.gz,vcf.gz.tbi}', mode: 'copy'

    input:
    tuple val(id), val(source), file('sage.vcf.gz'), file('sage.vcf.gz.tbi') from somaticSageVcfs

    output:
    tuple val(id), val(source), file('*.sage.pave.norm.no_common.quality.vcf.gz'), file('*.sage.pave.norm.no_common.quality.vcf.gz.tbi') into processedSageVcfs

    script:
    """
    bcftools view \
      -i 'FMT/AF[0:*] >= 0.1' \
      sage.vcf.gz | \
    bcftools view \
      -e 'INFO/TIER = "LOW_CONFIDENCE" && QUAL < 240' | \
    bcftools view \
      -e 'INFO/TIER = "HIGH_CONFIDENCE" && QUAL < 170' | \
    bgzip -c > ${id}_${source}.sage.pave.norm.no_common.quality.vcf.gz

    tabix -p vcf ${id}_${source}.sage.pave.norm.no_common.quality.vcf.gz
    """
}


////////////////////////////////////////////////////
/* --        Filter Structural Variants        -- */
////////////////////////////////////////////////////

process listGridssVcfs {
    input:
    tuple val(id), val(source), file('gripss.vcf.gz'), file('gripss.vcf.gz.tbi') from grippsPanel.map{ [it.id, it.source, file(it.gripss_vcf), file(it.gripss_index)] }

    output:
    tuple file('*.gripss.vcf.gz'), file('*.gripss.vcf.gz.tbi') into gripssVcfsForPanel
    path "*_path.txt" into gripssPaths

    script:
    """
    printf "%s_%s" "${id}" "${source}" > sample_name.txt

    bcftools reheader \
      --samples sample_name.txt \
      --output ${id}_${source}.gripss.vcf.gz \
      gripss.vcf.gz

    tabix -p vcf ${id}_${source}.gripss.vcf.gz

    printf "%s_%s.gripss.vcf.gz" "${id}" "${source}" > ${id}_${source}_gridss_path.txt
    """
}

process createGridssPON {
    publishDir "${params.outDir}/", pattern: '*_panel.{vcf.gz,vcf.gz.tbi}', mode: 'copy'

    input:
    file(vcfs) from gripssVcfsForPanel.collect()
    path 'vcf_paths.txt' from gripssPaths.collectFile(name: 'gripss_vcf_paths.txt', newLine: true)

    output:
    path 'gridss_panel.vcf.gz' into gridssPanelVcf
    path 'gridss_panel.vcf.gz.tbi' into gridssPanelIndex

    script:
    """
    bcftools merge \
      --apply-filters 'PASS' \
      --file-list vcf_paths.txt \
      --merge all | \
    bcftools +setGT \
      -- --target-gt q --new-gt c:0/1 -i 'FMT/AF>0' | \
    bcftools +fill-tags \
      -- --tags NS | \
    bcftools view \
      -i 'INFO/NS > 2' | \
    bgzip -c > gridss_panel.vcf.gz

    tabix -p vcf gridss_panel.vcf.gz
    """
}

process filterGridssVcfs {
    publishDir "${params.outDir}/${id}_${source}/", pattern: '*.{vcf.gz,vcf.gz.tbi}', mode: 'copy'

    input:
    tuple val(id), val(source), file('gripss.vcf.gz'), file('gripss.vcf.gz.tbi') from grippsInfo.map{ [it.id, it.source, file(it.gripss_vcf), file(it.gripss_index)] }
    path 'gridss_panel.vcf.gz' from gridssPanelVcf.first()
    path 'gridss_panel.vcf.gz.tbi' from gridssPanelIndex.first()

    output:
    tuple val(id), val(source), file('*.gripss.processed.vcf.gz'), file('*.gripss.processed.vcf.gz.tbi') into processedGridssVcfs

    script:
    """
    bcftools isec \
      gripss.vcf.gz \
      gridss_panel.vcf.gz \
      --complement \
      --collapse all \
      -w1 | \
    bcftools view \
      -f 'PASS' \
      -i 'FMT/AF[0:*] >= 0.1 && QUAL >= 1000' | \
    bcftools view \
      -i 'INFO/AS > 0 && INFO/RAS > 0' | \
    bgzip -c > ${id}_${source}.gripss.processed.vcf.gz

    tabix -p vcf ${id}_${source}.gripss.processed.vcf.gz
    """
}

process categorizeGridssVcfs {
    publishDir "${params.outDir}/${id}_${source}/", pattern: '*chord.txt', mode: 'copy'

    input:
    tuple val(id), val(source), file('gripss.processed.vcf.gz'), file('gripss.processed.vcf.gz.tbi') from processedGridssVcfs

    output:
    tuple val(id), val(source), file('*chord.txt') into chordSvDfs

    script:
    rscript = "$baseDir/bin/gridss_simple_events.R"
    """
    Rscript ${rscript} gripss.processed.vcf.gz ${id}.gripss.chord.txt
    """
}


////////////////////////////////////////////////////
/* --                Run CHORD                 -- */
////////////////////////////////////////////////////

process chord {
  input:
  tuple val(id), val(source), file('sage.vcf.gz'), file('sage.vcf.gz.tbi'), file('svs.txt') from processedSageVcfs.combine(chordSvDfs, by:[0, 1])

  output:
  path '*_contexts.txt' into chordContexts
  path '*_chord.txt' into chordPredictions

  script:
  rscript = "$baseDir/bin/chord.R"
  """
  Rscript ${rscript} ${id} ${source} sage.vcf.gz svs.txt
  """
}

chordContexts.collectFile(name: "${params.outDir}/chord_mutation_contexts.txt", sort: { file -> file.text }, keepHeader: true, skip: 1)
chordPredictions.collectFile(name: "${params.outDir}/chord_predictions.txt", sort: { file -> file.text }, keepHeader: true, skip: 1)


