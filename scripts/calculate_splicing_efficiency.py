from statistics import mean
import csv

class GenicRegion():

    def __init__(self, chromo, start, end, name, strand, score):
        self.chromo
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.name = name

        assert start < end
        assert strand in ['+', '-']
    
    @property
    def reads_per_pase(self):
        return len(self) / self.score
    
    def __len__(self):
        return self.end - self.start
    
    def __hash__(self):
        hash(''.join(self.__dict__.values()))
    
    
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
    

# def read_next_intersection_intron(csv_reader):
#     '''Take a csv Reader object of the input intersection of intron and exon
#     regions and return an Intron object by reading the data from the next two
#     lines. This will not work if there are not two exons for every intron
#     which I think is an incorrect assumption so need to correct that actually. 

#     Args:
#         csv_reader (Reader): csv.Reader object pointed to input file.

#     Returns:
#         Intron: Intron object instance.
#     '''
#     exon_1_inter = next(csv_reader)
#     exon_2_inter = next(csv_reader)

#     # make sure intron is the same
#     for i in range(0, 4):  # first 4 columns describe intron
#         assert exon_1_inter[i] == exon_2_inter[i]
    
#     exon_1 = GenicRegion(*exon_1_inter[4:])
#     exon_2 = GenicRegion(*exon_2_inter[4:])

#     intron = Intron(exon_1, exon_2, *exon_1_inter[:4])

#     return intron


def match_introns_to_exons(input_rows):
    intron_exon_dict = {}
    for each_row in input_rows:
        intron = tuple(each_row[:4])
        exon = tuple(each_row[5:])
        if intron in intron_exon_dict:
            intron_exon_dict[intron].append(exon)
        else:
            intron_exon_dict[intron] = [exon]
        
    return intron_exon_dict


def remove_unflanked_introns(intron_exon_dict):
    for intron, exons in intron_exon_dict.items():
        if len(exons) < 2:
            intron_exon_dict.pop(intron)
        elif len(exons) >= 3:
            raise Exception('Found intron with > 2 exons!')
        
    return intron_exon_dict


def make_intron_objects_from_dict(intron_exon_dict):
    for intron, exons in intron_exon_dict.items():
        exon_1, exon_2 = GenicRegion(*exons[0]), GenicRegion(*exons[1])
        intron = Intron(exon_1, exon_2, *intron)
        yield intron


def write_intron_SE_to_output(introns, output_file):
    with open(output_file, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        for each_intron in introns:
            write.writerow(
                intron.to_bed_line()
            )
    return output_file



def main():
    input_file = str(snakemake.input)
    output_file = str(snakemake.output)

    input_rows = [row for row in csv.reader(open(input_file), delimiter='\t')]
    intron_exon_dict = match_introns_to_exons(input_rows)
    intron_exon_dict = remove_unflanked_introns(intron_exon_dict)
    introns = make_intron_objects_from_dict(intron_exon_dict)
    write_intron_SE_to_output(introns, output_file)

    