#!/usr/bin/env bash

wget -O Macfas_6.0.115.gff3.gz ftp://ftp.ensembl.org/pub/release-115/gff3/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_6.0.115.gff3.gz

# MT appears in the gff3 but not in the macFas6 assembly files.
# Remove all MT features
# Using a literal tab character to support different versions of zgrep
tmpfile=$(mktemp)
zgrep -v "^MT	" Macfas_6.0.115.gff3.gz | gzip > $tmpfile
mv $tmpfile Macfas_6.0.115.gff3.gz
