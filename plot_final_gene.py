from gtf import GeneTable
import numpy as np
from matplotlib import pyplot as plt
import csv


#bamfile_template = ('/home/vighnesh/data/yeast/bamfiles/accepted_hits_%d.bam')


bamfile = '/home/vighnesh/data/yili/R1.bam'
gt = GeneTable(None, verbose=True, convert_roman=True)

total_points = 0
positive_points = 0
with open('/home/vighnesh/data/yili/pombe_origin.csv') as csv_file:
    csv_reader = csv.reader(csv_file)

    for line in csv_reader:
        chromosome, start, end = map(int, line)
        if chromosome == 1:
            score = gt.get_transition_significance(bamfile, 'I', start, end)
            if score > 0:
                positive_points += 1
            total_points += 1
            print(positive_points, total_points)
