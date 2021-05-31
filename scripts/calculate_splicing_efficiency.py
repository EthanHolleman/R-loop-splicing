from statistics import mean
import csv

class GenicRegion():

    def __init__(self, chromo, start, end, name, strand, score):
        self.chromo = chromo
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.name = name

        assert start < end
        assert strand in ['+', '-'], f'Tried to set strand to {strand}'
    
    @property
    def score(self):
        return self._score
    
    @score.setter
    def score(self, new_score):
        try:
            self._score = float(new_score)
        except Exception:
            self._score = 0
    
    @property
    def reads_per_base(self):
        return  self.score / len(self)
    
    def __len__(self):
        return self.end - self.start
    
    def __hash__(self):
        hash(''.join(self.__dict__.values()))
    
    
class Intron(GenicRegion):

    def __init__(self, upstream_exon, downstream_exon, *args, **kwargs):
        self.upstream_exon = upstream_exon
        self.downstream_exon = downstream_exon
        super().__init__(*args, **kwargs)
    
    @property
    def splicing_efficiency(self):
        try:
            return 1 - (self.reads_per_base / mean(
                    (
                    self.upstream_exon.reads_per_base,
                    self.downstream_exon.reads_per_base
                    )
                )
            )
        except ZeroDivisionError:
            # this would occur is there is no coverage
            # over the exonic regions
            return 'NA'
    
    def to_bed_line(self):
        return [self.chromo, 
                self.start, 
                self.end, 
                self.name, 
                self.splicing_efficiency,
                self.strand
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
        intron = tuple(each_row[:7])
        exon = tuple(each_row[7:])
        intron_gene = intron[3].split('.')[0]
        exon_gene = exon[3].split('.')[0]
        if intron_gene == exon_gene:
            if intron in intron_exon_dict:
                intron_exon_dict[intron].append(exon)
            else:
                intron_exon_dict[intron] = [exon]
        
    return intron_exon_dict


def remove_unflanked_introns(intron_exon_dict):
    trimmed_introns = {}
    for intron, exons in intron_exon_dict.items():
        if len(exons) == 2:
            trimmed_introns.update({intron: exons})
    for intron, exons in trimmed_introns.items():
        assert len(exons) == 2
    return trimmed_introns


def make_intron_objects_from_dict(intron_exon_dict):
    for intron, exons in intron_exon_dict.items():
        # have to slice out the EXON keyword from bed line
        exon_1 = GenicRegion(*exons[0][:5] + exons[0][6:])
        exon_2 = GenicRegion(*exons[1][:5] + exons[1][6:])
        
        assert exon_1.start != exon_2.start
        assert exon_1.end != exon_2.end
        assert exon_1.chromo == exon_2.chromo
        assert exon_1.name == exon_2.name

        intron = Intron(exon_1, exon_2, *intron[:5] + intron[6:])
        yield intron


def write_intron_SE_to_output(introns, output_file):
    with open(output_file, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        for each_intron in introns:
            writer.writerow(
                each_intron.to_bed_line()
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

if __name__ == '__main__':
    main()

    