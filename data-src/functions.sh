# ensure all primary (chr1~22, X, Y, M) chromosome features float to the top of the file
# for backwards compatiblity with older versions of GenomeKit
wget_ordered_chroms() {
  finalfile=$1
  url=$2
  tmpfile=$(mktemp)
  regex=$3
  wget -O ${tmpfile} ${url}
  (zgrep --perl-regexp ${regex} ${tmpfile}; zgrep --perl-regexp --invert-match ${regex} ${tmpfile}) | gzip > ${finalfile}
  rm ${tmpfile}
}

function wget_gencode() {
  wget_ordered_chroms $1 $2 "^chr(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|M)\t"
}


function wget_ncbi_hg19() {
  # source: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chromAlias.txt
  #chr1: NC_000001.10
  #chr10: NC_000010.10
  #chr11: NC_000011.9
  #chr12: NC_000012.11
  #chr13: NC_000013.10
  #chr14: NC_000014.8
  #chr15: NC_000015.9
  #chr16: NC_000016.9
  #chr17: NC_000017.10
  #chr18: NC_000018.9
  #chr19: NC_000019.9
  #chr2: NC_000002.11
  #chr20: NC_000020.10
  #chr21: NC_000021.8
  #chr22: NC_000022.10
  #chr3: NC_000003.11
  #chr4: NC_000004.11
  #chr5: NC_000005.9
  #chr6: NC_000006.11
  #chr7: NC_000007.13
  #chr8: NC_000008.10
  #chr9: NC_000009.11
  #chrM: NC_001807.4
  #chrX: NC_000023.10
  #chrY: NC_000024.9

  wget_ordered_chroms $1 $2 "^(NC_0000(01\.10|10\.10|11\.9|12\.11|13\.10|14\.8|15\.9|16\.9|17\.10|18\.9|19\.9|02\.11|20\.10|21\.8|22\.10|03\.11|04\.11|05\.9|06\.11|07\.13|08\.10|09\.11|23\.10|24\.9))|(NC_001807\.4)\t"
}

function wget_ncbi_hg38() {
  # source: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.chromAlias.txt
  # chr1: NC_000001.11
  # chr10: NC_000010.11
  # chr11: NC_000011.10
  # chr12: NC_000012.12
  # chr13: NC_000013.11
  # chr14: NC_000014.9
  # chr15: NC_000015.10
  # chr16: NC_000016.10
  # chr17: NC_000017.11
  # chr18: NC_000018.10
  # chr19: NC_000019.10
  # chr2: NC_000002.12
  # chr20: NC_000020.11
  # chr21: NC_000021.9
  # chr22: NC_000022.11
  # chr3: NC_000003.12
  # chr4: NC_000004.12
  # chr5: NC_000005.10
  # chr6: NC_000006.12
  # chr7: NC_000007.14
  # chr8: NC_000008.11
  # chr9: NC_000009.12
  # chrM: NC_012920.1
  # chrX: NC_000023.11
  # chrY: NC_000024.10
  wget_ordered_chroms $1 $2 "^(NC_0000(01\.11|10\.11|11\.10|12\.12|13\.11|14\.9|15\.10|16\.10|17\.11|18\.10|19\.10|02\.12|20\.11|21\.9|22\.11|03\.12|04\.12|05\.10|06\.12|07\.14|08\.11|09\.12|23\.11|24\.10))|(NC_012920\.1)\t"
}

function wget_ncbi_macfas5() {
  # source: https://hgdownload.soe.ucsc.edu/goldenPath/macFas5/bigZips/macFas5.chromAlias.txt
  # chr1: NC_022272.1
  # chr10: NC_022281.1
  # chr11: NC_022282.1
  # chr12: NC_022283.1
  # chr13: NC_022284.1
  # chr14: NC_022285.1
  # chr15: NC_022286.1
  # chr16: NC_022287.1
  # chr17: NC_022288.1
  # chr18: NC_022289.1
  # chr19: NC_022290.1
  # chr2: NC_022273.1
  # chr20: NC_022291.1
  # chr3: NC_022274.1
  # chr4: NC_022275.1
  # chr5: NC_022276.1
  # chr6: NC_022277.1
  # chr7: NC_022278.1
  # chr8: NC_022279.1
  # chr9: NC_022280.1
  # chrM: NC_012670.1
  # chrX: NC_022292.1
  wget_ordered_chroms $1 $2 "^(NC_0222(72\.1|81\.1|82\.1|83\.1|84\.1|85\.1|86\.1|87\.1|88\.1|89\.1|90\.1|73\.1|91\.1|74\.1|75\.1|76\.1|77\.1|78\.1|79\.1|80\.1|92\.1))|(NC_012670\.1)\t"
}


function wget_ncbi_mm10() {
  # source: https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/p6/mm10.p6.chromAlias.txt
  # chr1: NC_000067.6
  # chr10: NC_000076.6
  # chr11: NC_000077.6
  # chr12: NC_000078.6
  # chr13: NC_000079.6
  # chr14: NC_000080.6
  # chr15: NC_000081.6
  # chr16: NC_000082.6
  # chr17: NC_000083.6
  # chr18: NC_000084.6
  # chr19: NC_000085.6
  # chr2: NC_000068.7
  # chr3: NC_000069.6
  # chr4: NC_000070.6
  # chr5: NC_000071.6
  # chr6: NC_000072.6
  # chr7: NC_000073.6
  # chr8: NC_000074.6
  # chr9: NC_000075.6
  # chrM: NC_005089.1
  # chrX: NC_000086.7
  # chrY: NC_000087.7
  wget_ordered_chroms $1 $2 "^(NC_0000(67\.6|76\.6|77\.6|78\.6|79\.6|80\.6|81\.6|82\.6|83\.6|84\.6|85\.6|68\.7|69\.6|70\.6|71\.6|72\.6|73\.6|74\.6|75\.6|86\.7|87\.7))|(NC_005089\.1)\t"
}

function wget_ncbi_mm39() {
  # source: https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chromAlias.txt
  # chr1: NC_000067.7
  # chr10: NC_000076.7
  # chr11: NC_000077.7
  # chr12: NC_000078.7
  # chr13: NC_000079.7
  # chr14: NC_000080.7
  # chr15: NC_000081.7
  # chr16: NC_000082.7
  # chr17: NC_000083.7
  # chr18: NC_000084.7
  # chr19: NC_000085.7
  # chr2: NC_000068.8
  # chr3: NC_000069.7
  # chr4: NC_000070.7
  # chr5: NC_000071.7
  # chr6: NC_000072.7
  # chr7: NC_000073.7
  # chr8: NC_000074.7
  # chr9: NC_000075.7
  # chrM: NC_005089.1
  # chrX: NC_000086.8
  # chrY: NC_000087.8
  wget_ordered_chroms $1 $2 "^(NC_0000(67\.7|76\.7|77\.7|78\.7|79\.7|80\.7|81\.7|82\.7|83\.7|84\.7|85\.7|68\.8|69\.7|70\.7|71\.7|72\.7|73\.7|74\.7|75\.7|86\.8|87\.8))|(NC_005089\.1)\t"
}

function wget_ensembl() {
  wget_ordered_chroms $1 $2 "^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|MT)\t"
}
