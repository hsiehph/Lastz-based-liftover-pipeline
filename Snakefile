shell.prefix("source config.sh;")

configfile: "config.json"

if not os.path.exists("log"):
	os.makedirs("log")

import os
import glob
import csv

TN=config["target_names_file"]
REGIONS=[line.rstrip('\n') for line in open(config["target_names_file"])]


TMPDIR = config["tmp_dir"]
POST   = config["lastz_settings"]

TARGET  = config["target_path"] + config["target_name"]
TARGETI = TARGET + ".fai"
QUERY   = config["query_path"] + config["query_name"]
QUERYI  = QUERY + ".fai"

LQ=config["local_path"] + config["query_name"]

RS_LT= "mkdir -p %s ; rsync  --ignore-existing --bwlimit=50000 -v %s %s" % (config["local_path"], TARGET, config["local_path"])
RS_LTI="mkdir -p %s ; rsync  --ignore-existing --bwlimit=50000 -v %s %s" % (config["local_path"], TARGETI, config["local_path"])
RS_LQ= "mkdir -p %s ; rsync  --ignore-existing --bwlimit=50000 -v %s %s" % (config["local_path"], QUERY, config["local_path"])
RS_LQI="mkdir -p %s ; rsync  --ignore-existing --bwlimit=50000 -v %s %s" % (config["local_path"], QUERYI, config["local_path"])

FAIDX="samtools faidx %s%s" % (config["local_path"], config["target_name"])

rule all: 
    input: "net/all-no-class.net", "filt_net/all-no-class.filt.net"

rule filtNet:
     message : "INFO: filtering net file"
     input:  "net/all-no-class.net"
     output: "filt_net/all-no-class.filt.net"
     params: sge_opts=config["cluster_settings"]["lite"]
     shell:  "netFilter {input} > {output} && netToBed {output} filt_net/all-no-class-filter.bed"

rule net:
     message : "INFO: netting merged chain"
     input:  CHAIN="merged_chains/all.chain.filter.prenet.chain", TS="stats/target.sizes", QS="stats/query.sizes"
     output: "net/all-no-class.net"
     params: sge_opts=config["cluster_settings"]["lite"]
     shell:  "chainNet {input.CHAIN}  -minSpace=1 {input.TS} {input.QS}  stdout /dev/null  |  netSyntenic stdin {output}"

rule filter:
     message : "INFO: filtering merged chain"
     input:  CHAINS="chained_psl/chained.chain", TS="stats/target.sizes", QS="stats/query.sizes"
     params: sge_opts=config["cluster_settings"]["lite"]
     output: "merged_chains/all.chain.filter.prenet.chain"
     shell:  "chainPreNet {input.CHAINS} {input.TS} {input.QS} {output}"

rule chain   :
     message : "INFO: chaining alignment"
     input   :  PSL="merged_psl/merged.psl", T2BIT="2bits/target.2bit", Q2BIT="2bits/query.2bit"
     output  : temp("chained_psl/chained.chain")
     params  : sge_opts=config["cluster_settings"]["heaviest"]
     shell   : "axtChain -linearGap=medium -psl {input.PSL} {input.T2BIT} {input.Q2BIT} {output}"

rule mergePsl:
     message : "INFO: merging psl"
     input   : expand("raw_psl/{contig}.psl", contig=REGIONS)
     output  : "merged_psl/merged.psl"
     params  : sge_opts=config["cluster_settings"]["lite"]
     shell   : "pslCat {input} > {output}"

rule lav_to_psl:
     message   : "INFO: coverting lav to psl"
     input     : "raw_lav/{contig}.lav"
     output    : "raw_psl/{contig}.psl"
     params    :  sge_opts=config["cluster_settings"]["heavy"]
     shell     :  "lavToPsl {input} {output}"

rule runLastZ:
     message : "INFO: running lastz, please wait, it could take a while"
     input   :  T={TARGET}, Q={QUERY}
     output  : "raw_lav/{contig}.lav"
     params  : sge_opts=config["cluster_settings"]["heaviest"]
     shell   :  "{RS_LT} ; {RS_LQ} ; {RS_LTI} ; {RS_LQI} ; {FAIDX} {wildcards.contig} > {TMPDIR}/{wildcards.contig}.fasta ; lastz  --ambiguous=iupac {TMPDIR}/{wildcards.contig}.fasta {LQ} {POST} > {output}"

rule stats:
     message : "INFO: getting target and query stats"
     input:  T={TARGET}, Q={QUERY}, TI={TARGETI}, QI={QUERYI}
     output: "stats/target.sizes", "stats/query.sizes"
     params: sge_opts=config["cluster_settings"]["lite"]
     shell:  "faSize {input.T} -detailed > {output[0]} && faSize {input.Q} -detailed > {output[1]}"

rule twobit:
     message : "INFO: creating 2bit files"
     input:  T={TARGET}, Q={QUERY}
     output: "2bits/target.2bit", "2bits/query.2bit"
     params: sge_opts=config["cluster_settings"]["heavy"]
     shell:  "faToTwoBit {input.T} {output[0]} && faToTwoBit {input.Q} {output[1]}"

rule faIndexT:
     message : "INFO: indexing target"
     input:  T={TARGET}
     output: {TARGETI}
     params: sge_opts=config["cluster_settings"]["heavy"]
     shell:  "samtools faidx {input.T}"

rule faIndexQ:
     message : "INFO: indexing query"
     input:  Q={QUERY}
     output: {QUERYI}
     params: sge_opts=config["cluster_settings"]["heavy"]
     shell:  "samtools faidx {input.Q}"
