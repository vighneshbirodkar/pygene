from gtf import GeneTable
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys


#start = 518619
#end = 519119

start = int(sys.argv[1])
end = int(sys.argv[2])

bamfile = '/home/vighnesh/data/yili/R1.bam'
gt = GeneTable(None, verbose=True, convert_roman=False)
values = gt.get_counts_by_location(bamfile, 3, start, end)

fig = plt.figure(figsize=(16, 9))


plt.plot(range(start, end), values, linewidth=2)
plt.grid(True)
plt.ylabel('Counts')
plt.xlabel('Genome Index')
plt.axhline(c='r', linestyle='--')
plt.xticks([start, end])
plt.gca().set_xticklabels([start, end])
plt.xlim(start, end)

plt.savefig('figure.png')
plt.show()
