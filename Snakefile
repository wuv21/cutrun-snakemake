import os
import sys
from collections import defaultdict
import itertools

snakefileDir = os.path.dirname(workflow.snakefile)
scriptDir = snakefileDir + "/scripts"

###############################################################################
# input paths dict
###############################################################################
inputPaths = defaultdict(str)
for inKey in config["paths"]:
  inPath = os.path.normpath(config["paths"][inKey])
  inputPaths[inKey] = inPath

###############################################################################
# make all out directories if needed
###############################################################################
outDirs = defaultdict(str)
for outDirKey in config["out_dirs"]:
  outDir = os.path.normpath(config["out_dirs"][outDirKey])

  if not os.path.exists(outDir):
    os.mkdir(outDir)

  outDirs[outDirKey] = outDir

###############################################################################
# resource allocation
###############################################################################
availThreads = config["general"]["threads"]
availMem = config["general"]["memory"]

###############################################################################
# determine what are the final outputs and check if req'd files exist
###############################################################################
def generateCombo(prev, step, stepN, conditions):
  op = step["type"]
  
  if op == "append" and len(prev) == 0:
    return(conditions[step["field"]])

  elif op == "append":
    newlist = []
    pass

  elif op == "product":
    selectCondList = [conditions[x] for x in step["field"]]
    
    return(list(itertools.product(prev, *selectCondList)))

  elif op == "groupedProduct" and prev != None:
    groupBy = conditions[step["groupBy"]]
    newlist = []

    for i in range(0, len(prev), len(groupBy)):
      section = prev[i:i + len(groupBy)]
      
      for suffix in conditions[step["field"]]:
        newsection = tuple(map(lambda x: x + (suffix,),  section))
        newlist.extend(newsection)

    return(newlist)


def final_output(modes):
  # file check
  if modes["withEColiNormalization"] and not os.path.exists(inputPaths["ecoli_ref_dir"]):
    raise Exception("E.coli normalization mode is turned but no reference directory specified")

  if not os.path.exists(inputPaths["fastq_dir"]):
    raise Exception("Fastq directory does not exist.")

  if not os.path.exists(inputPaths["ref_dir"]):
    raise Exception("No genome reference exists.")

  # generate auto sample naming
  if modes["autoloadSamplesRecursively"]:
    conditions = config["samples_auto"]["conditions"]
    bclSampleIds = config["samples_auto"]["bcl2fastqSampleID"]
    sampleIndexingStart = config["samples_auto"]["sampleIndexingStart"]

    if bclSampleIds != None:
      combs = []
      for i, step in enumerate(bclSampleIds):
        combs = generateCombo(combs, step, i, conditions)
        print(combs)
      
      
      combs2 = []
      for i, x in enumerate(combs):
        if type(x) is not tuple and type(x) is not list:
          x = (x,)
        
        combs2.append(x + ("S" + str(i + sampleIndexingStart),))
      
    finalCombs = [config["samples_auto"]["conditionDelimiter"].join(map(lambda a: str(a), x)) for x in combs2]

  else:
    raise Exception("I'm sorry :( This functionality isn't ready yet")
  
  print(finalCombs)
  outs = []

  # test only
#  for file in os.scandir("test/fastq"):
#    os.remove(file.path)
#  
#  for x in finalCombs:
#    for i in range(1, 3):
#      with open(f"test/fastq/{x}_R{i}_001.fastq.gz", "w"):
#        pass

 
  outs.append(
    expand(
      f"outs/{{ref}}_stats/{{smpl}}.stats",
      smpl = finalCombs,
      ref = ["ecoli", "actual"]))
 
  outs.append(
    expand(
      f"outs/bw/{{smpl}}.bw",
      smpl = finalCombs))
 

  return(outs)

###############################################################################
# rules
###############################################################################

rule all: 
  input:
    *final_output(config["general"]["run_modes"])


initialFastqFiles = lambda wildcards: expand(inputPaths["fastq_dir"] + "/" + wildcards.smpl + "_R{read}_" + config["samples_auto"]["fileSuffix"], read = [1,2])

rule align_actual:
  input:
    initialFastqFiles
  output:
    "outs/actual_sam/{smpl}.sam"
  params:
    s = "{smpl}",
    o = "outs/actual_sam",
    r = inputPaths["ref_dir"]
  threads: availThreads["align"]
  resources:
    mem_mb = 120000
  shell:
    """
      sh {scriptDir}/bwa_fx.sh \
       -s {params.s} \
       -x {input[0]} \
       -y {input[1]} \
       -o {params.o} \
       -r {params.r}
    """

rule align_ecoli:
  input:
    initialFastqFiles
  output:
    "outs/ecoli_sam/{smpl}.sam"
  params:
    s = "{smpl}",
    o = "outs/ecoli_sam",
    r = inputPaths["ecoli_ref_dir"]
  threads: availThreads["align"]
  shell:
    """
      sh {scriptDir}/bwa_fx.sh \
       -s {params.s} \
       -x {input[0]} \
       -y {input[1]} \
       -o {params.o} \
       -r {params.r}
    """


rule sam_to_sorted_bam:
  input:
    "outs/{ref}_sam/{smpl}.sam"
  output:
    "outs/{ref}_bam/{smpl}.bam"
  params:
    bamDir = "outs/{ref}_bam"
  threads: availThreads["sortBam"]
  shell:
    """
    samtools sort -@ {threads} \
      -T {params.bamDir}/{wildcards.smpl}.tmp \
      -O bam \
      -o {output} \
      {input}
    """

rule dedupe:
  input:
    rules.sam_to_sorted_bam.output
  output:
    "outs/{ref}_bam/{smpl}.dedup.bam"
  params:
    picard = inputPaths["picardExe"],
    metricsFile = "outs/{ref}_bam/{smpl}.metrics.txt"
  shell:
    """
      java -jar {params.picard} MarkDuplicates \
        --INPUT {input} \
        --OUTPUT {output} \
        --METRICS_FILE {params.metricsFile}
    """
    
rule index_bam:
  input:
    rules.sam_to_sorted_bam.output
  output:
    touch(".tmp_snek/{ref}_{smpl}.index.done")
  threads: availThreads['index']
  shell:
    """
    samtools index {input} -@ 12
    """

rule generate_stats:
  input:
    rules.dedupe.output
  output:
    "outs/{ref}_stats/{smpl}.stats"
  shell:
    """
    samtools flagstats {input} > {output}
    """

rule dedeny:
  input:
    tmp = rules.index_bam.output,
    real = rules.dedupe.output
  output:
    "outs/{ref}_bam/{smpl}.dedup.dedeny.bam"
  params:
    denylist = inputPaths["denylist"]
  wildcard_constraints:
    ref = "actual"
  shell:
    """
    sh {scriptDir}/dedeny_fx.sh \
      -i {input.real} \
      -d {params.denylist} \
      -o {output}
    """

rule index_dedeny_bam:
  input:
    rules.dedeny.output
  output:
    touch(".tmp_snek/{ref}_{smpl}.dedenyIndex.done")
  threads: availThreads['index']
  wildcard_constraints:
    ref = "actual"
  shell:
    """
    samtools index {input} -@ 12
    """


rule compile_stats:
  input:
    tmp = expand(".tmp_snek/{ref}_{{smpl}}.dedenyIndex.done", ref = ["actual"]),
    stats = expand("outs/{ref}_stats/{{smpl}}.stats", ref = ["ecoli", "actual"]),
    actualBam = expand("outs/{ref}_bam/{{smpl}}.dedup.dedeny.bam", ref = ["actual"])
  output:
    "outs/bw/{smpl}.bw"
  threads: availThreads["compileStats"]
  shell:
    """
      actualUniq=$(cat {input.stats[1]} | grep -oP "\d+ (?=\+ \d+ primary mapped)")
      ecoliUniq=$(cat {input.stats[0]} | grep -oP "\d+ (?=\+ \d+ primary mapped)")
      smplNormFactor=$(echo "scale = 4; 1 / ($ecoliUniq / $actualUniq * 100)" | bc -l)

      echo "Uniq actual reads: $actualUniq"
      echo "Uniq E.coli reads: $ecoliUniq"
      echo "Normalization factor: $smplNormFactor"

      bamCoverage \
        -b {input.actualBam[0]} \
        -o {output} \
        --MNase \
        --scaleFactor=$smplNormFactor \
        --normalizeUsing="CPM" \
        --samFlagExclude=1024 \
        --numberOfProcessors=8
    """

