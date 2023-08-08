while getopts i:d:o: flag
do
  case "${flag}" in
    i) input=${OPTARG};;
    d) denylist=${OPTARG};;
    o) output=${OPTARG};;
esac
done

bedtools intersect -v -abam ${input} -b ${denylist} > ${output}

