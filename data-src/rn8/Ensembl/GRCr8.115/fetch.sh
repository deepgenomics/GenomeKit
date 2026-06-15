#!/usr/bin/env bash

wget -O GRCr8.115.gff3.gz ftp://ftp.ensembl.org/pub/release-115/gff3/rattus_norvegicus/Rattus_norvegicus.GRCr8.115.gff3.gz

# MT appears in the gff3 but not in the GRCr8 assembly files.
# Remove all MT features
# Using a literal tab character to support different versions of zgrep
tmpfile=$(mktemp)
zgrep -v "^MT	" GRCr8.115.gff3.gz | gzip > $tmpfile
mv $tmpfile GRCr8.115.gff3.gz
