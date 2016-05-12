from gtf import GeneTable
import numpy as np
from matplotlib import pyplot as plt


#bamfile_template = ('/home/vighnesh/data/yeast/bamfiles/accepted_hits_%d.bam')
gtf_file = ('Saccharomyces_cerevisiae.R64-1-1.84.gtf')
bamfile = 'R1.bam'
gt = GeneTable(gtf_file, verbose=True, num_flag=True)


a = gt.get_pre_mid_post(bamfile, 'I', 'YAR053W')

b = gt.get_counts_by_location(bamfile, 'I', 110578, 114578)

plt.plot(np.hstack(b))

plt.show()
