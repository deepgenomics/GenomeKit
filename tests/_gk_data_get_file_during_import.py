# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from genome_kit import gk_data

gk_data.get_file('test_file.txt')

print("_gk_data_get_file_during_import: this line should be unreachable!")
