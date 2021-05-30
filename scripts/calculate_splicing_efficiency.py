from statistics import mean


class GenicRegion():

    def __init__(self, chromo, start, end, name, strand, score):
        self.chromo
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.name = name
    
    @property
    def reads_per_pase(self):
        return len(self) / self.score
    
    def __len__(self):
        return self.end - self.start
    
    
class Intron(GenicRegion):

    def __init__(upstream_exon, downstream_exon, *args, **kwargs):
        self.upstream_exon = upstream_exon
        self.downstream_exon = downstream_exon
        super().__init__(*args, **kwargs)
    
    @property
    def splicing_efficiency(self):
        return 1 - (self.reads_per_base / mean(
                (
                self.upstream_exon.reads_per_base,
                self.downstream_exon.reads_per_base
                )
            )
        )
    
    def to_bed_line(self):
        return [self.chromo, 
                self.start, 
                self.end, 
                self.name, 
                self.splicing_efficiency
        ]


    