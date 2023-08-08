while getopts s:b:x:y:f:o:r: flag
do
  case "${flag}" in
    s) sample=${OPTARG};;
    x) fq1=${OPTARG};;
    y) fq2=${OPTARG};;
    o) outputDir=${OPTARG};;
    r) ref=${OPTARG};;
esac
done

flowcell=$(zcat ${fq1} | head -n 1 | perl -ne 'print "$1" while /@[A-Z0-9]+:[A-Z0-9]+:([A-Z0-9]+):/gi')

idpu="ID:${sample}:${flowcell}"

RG="@RG\\tID:${idpu}\\tPL:ILLUMINA\\tPU:${idpu}\\tLB:not_specified\\tSM:${sample}"

bwa mem \
  -k 19 \
  -T 30 \
  -R ${RG} \
  -t 12 \
  ${ref} \
  ${fq1} \
  ${fq2} > ${outputDir}/${sample}.sam

