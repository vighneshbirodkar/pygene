import commands
import numpy as np
from skimage.transform import resize


ROMAN_MAP = {'II': 2, 'VI': 6, 'XI': 11, 'IX': 9, 'I': 1, 'XII': 12, 'VII': 7,
             'IV': 4, 'XVI': 16, 'Mito': 0, 'XIII': 13, 'VIII': 8, 'XIV': 14,
             'V': 5, 'X': 10, 'III': 3, 'XV': 15}


def resize_signal(a, size):
    a = a.reshape(-1, 1)
    b = resize(a, (size, 1), preserve_range=True)
    return b[:, 0]


def samtools_reads_iter(text, gene_start, gene_end):
    "Produces iterator of (start, end) indices from samtools output"

    c = 0
    for line in text.split('\n'):
        c += 1
        words = line.split('\t')
        flag = int(words[1])
        if len(words) <= 1:
            continue

        alignment_start = int(words[3])
        alignment_length = len(words[9])

        yield alignment_start, alignment_start + alignment_length, flag


def get_samtools_view_command(filename, chromosome, start, end):
    return 'samtools view %s %s:%d-%d' % (filename, chromosome, start, end)


class GeneTable(object):

    def __init__(self, filename, interested_type='CDS', verbose=False,
                 convert_roman=True, pos_neg=True):
        self.verbose = verbose
        self.total_seqs = 0
        self.error_seqs = 0
        self.convert_roman = convert_roman
        self.pos_neg = pos_neg

        # Maps chromosome to gene table
        self.chromo_map = {}
        # Gene table maps gene name to start and end indices

        if filename:
            with open(filename) as gtf_file:
                for line in gtf_file:

                    if line.startswith('#'):
                        continue

                    words = line.split()
                    chromosome = words[0]
                    gene_type = words[2]

                    if gene_type != interested_type:
                        continue

                    gene_name = words[9].strip(';').strip('"')

                    if chromosome in self.chromo_map:
                        gene_map = self.chromo_map[chromosome]
                    else:
                        gene_map = {}
                        self.chromo_map[chromosome] = gene_map

                        gene_map = self.chromo_map[chromosome]
                        start = int(words[3])
                        stop = int(words[4])
                        gene_map[gene_name] = (start, stop)

        if verbose:
            total_genes = sum([len(table) for
                               table in self.chromo_map.values()])
            print('Parsed %d Chromosomes with %d total genes ' %
                  (len(self.chromo_map), total_genes))

    def get_chromosomes(self):
        return self.chromo_map.keys()

    def get_counts(self, chromosome, gene_name, bamfile, nbins=100,
                   offset=1000):

        counts = [0 for i in range(nbins)]

        gene_start, gene_end = self.chromo_map[chromosome][gene_name]
        gene_start -= 1000
        gene_end += 1000
        gene_length = gene_end - gene_start

        cmd = get_samtools_view_command(bamfile, chromosome,
                                        gene_start, gene_end)
        status, output = commands.getstatusoutput(cmd)

        if len(output) == 0:
            print('No reads in Chromosome %s Gene %s File %s' %
                  (chromosome, gene_name, bamfile))
            return counts

        if status != 0:
            print('Error running samtools')

        for (alignment_start,
             alignment_end) in samtools_reads_iter(output, gene_start,
                                                   gene_end):

            mid = (alignment_start + alignment_end)/2
            idx = (float(mid - gene_start)/gene_length)*nbins
            idx = int(idx)
            idx -= 1

            self.total_seqs += 1
            if idx < 0 or idx >= nbins:
                self.error_seqs += 1
                continue

           # print(alignment_start, alignment_length, gene_start, gene_end)
            counts[idx] += 1

        return counts

    def _get_choromsomes_from_spec(self, spec):
        "Get all chromosomes which match the specification"

        if spec == 'all':
            chromosomes = self.chromo_map.keys()
        elif isinstance(spec, str):
            chromosomes = [spec]
        chromosomes = spec
        return chromosomes

    def get_all_counts(self, bamfile, chromosomes='all', nbins=100):
        count_list = []

        chromosomes = self._get_choromsomes_from_spec(chromosomes)

        for chromo in chromosomes:
            if self.verbose:
                print('get_all_counts: Processing Chromosome %s' % chromo)
            for gene in self.chromo_map[chromo]:
                counts = self.get_counts(chromo, gene, bamfile, nbins)
                count_list.append(counts)

        return np.array(count_list)

    def get_error(self):
        return float(self.error_seqs)/self.total_seqs

    def get_genes_in_chromo(self, chromosome):
        return self.chromo_map[chromosome].keys()

    def get_utr_3_5_all(self, bamfile, chromosomes, offset=0):

        chromosomes = self._get_choromsomes_from_spec(chromosomes)

        start = np.zeros(2*offset)
        end = np.zeros(2*offset)
        count = 0

        for chromo in chromosomes:
            for gene in self.chromo_map[chromo]:
                s, e = self.get_utr_3_5_gene(bamfile, chromo, gene, offset)
                start += s
                end += e
                count += 1

        return (start/count, end/count)

    def get_pre_mid_post_all(self, bamfile, chromosome, offset=0,
                             normal_gene_length=3000):

        chromosomes = self._get_choromsomes_from_spec(chromosome)

        pre = np.zeros(offset)
        mid = np.zeros(normal_gene_length)
        post = np.zeros(offset)

        count = 0

        for chromo in chromosomes:
            for gene in self.chromo_map[chromo]:
                s, m, e = self.get_pre_mid_post(bamfile, chromo, gene, offset,
                                                normal_gene_length)
                pre += s
                mid += m
                post += e
                count += 1

        return (pre/count, mid/count, post/count)

    def get_pre_mid_post(self, bamfile, chromosome, gene_name, offset=0,
                         normal_gene_length=3000):

        start, end = self.chromo_map[chromosome][gene_name]
        pre_region = self.get_counts_by_location(bamfile, chromosome,
                                                 start - offset, start)
        mid_region = self.get_counts_by_location(bamfile, chromosome, start,
                                                 end)
        post_region = self.get_counts_by_location(bamfile, chromosome, end,
                                                  end + offset)
        mid_region = resize_signal(mid_region, normal_gene_length)

        return (pre_region, mid_region, post_region)

    def get_gene_count(self, bamfile, chromosome, gene_name):

        start, end = self.chromo_map[chromosome][gene_name]
        return self.get_counts_by_location(bamfile, chromosome,
                                           start, end)

    def get_counts_by_location(self, bamfile, chromosome, start, end):

        if self.convert_roman:
            chromosome = ROMAN_MAP[chromosome]

        cmd = get_samtools_view_command(bamfile, chromosome,
                                        start, end)
        status, output = commands.getstatusoutput(cmd)
        if status != 0:
            print('Error with samtools : %s' % output)
        array = np.zeros(end - start)
        for (s, e, flag) in samtools_reads_iter(output, start, end):

                sign = 1
                if self.pos_neg:
                    if (flag & 0x10):
                        sign = -1
                    else:
                        sign = 1

                s -= start
                e -= start

                s = np.clip(s, 0, array.shape[0] - 1)
                e = np.clip(e, 0, array.shape[0] - 1)

                array[s:e] += sign
        return array

    def get_utr_3_5_gene(self, bamfile, chromosome, gene_name, offset=500):
        "Return the counts array of reads that map to start and end of a gene"

        start_array = np.zeros(2*offset)
        end_array = np.zeros(2*offset)

        gene_start, gene_end = self.chromo_map[chromosome][gene_name]

        arrays = [start_array, end_array]
        locations = [gene_start, gene_end]
        for (array, gene_location) in zip(arrays, locations):

            start = gene_location - offset
            end = gene_location + offset
            cmd = get_samtools_view_command(bamfile, chromosome,
                                            start, end)
            status, output = commands.getstatusoutput(cmd)

            for (alignment_start,
                 alignment_end) in samtools_reads_iter(output, start, end):
                alignment_start -= start
                alignment_end -= start

                if alignment_start < 0 or alignment_end >= array.shape[0]:
                    continue

                array[alignment_start:alignment_end] += 1

        return start_array, end_array

    def get_all_gene_names(self, chromosome):
        d = self.chromo_map[chromosome]
        return d.keys()

    def get_transition_significance(self, bamfile, chromosome, start, end):
        values = self.get_counts_by_location(bamfile, chromosome, start, end)

        mid = len(values)/2
        left = values[:mid]
        right = values[mid:]

        up_left = np.sum(left[left > 0])
        down_right = np.sum(-right[right < 0])
        down_left = np.sum(-left[left < 0])
        up_right = np.sum(right[right > 0])

        assert up_left >= 0
        assert down_right >= 0
        assert down_left >= 0
        assert up_right >= 0

        numerator = up_left + down_right
        denominator = down_left + up_right

        if numerator == 0:
            return -np.inf
        if denominator == 0:
            return np.inf
        else:
            significance = np.log2(float(numerator)/float(denominator))
            return significance
