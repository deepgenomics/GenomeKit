from importlib.metadata import entry_points

from genome_kit import Interval, Variant
from genome_kit.genome_annotation import Gene, Transcript, Exon, Intron, Cds, Utr
from genome_kit.genome import Genome
from genome_kit import marshal

# dynamically inject functions to the marshal module for backwards compatibility
eps = list(entry_points(group="genomekit.plugins.compat", name="module"))
if len(eps) > 0:
    objects = list(eps)[0].load()(Gene, Transcript, Exon, Intron, Cds, Interval, Variant, Genome)
    for obj in objects:
        setattr(marshal, obj.__name__, obj)

# monkey patch serialization functions in GK classes
eps = list(entry_points(group="genomekit.plugins.compat", name="monkeypatch"))
if len(eps) > 0:
    list(eps)[0].load()(Gene, Transcript, Exon, Intron, Cds, Utr, Interval, Variant, Genome)
