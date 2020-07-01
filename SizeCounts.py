# author:siayouyang
from optparse import OptionParser
import collections
import os

##OptionParser
parser = OptionParser()
parser.add_option("-s", "--sam", dest="sam_name", type=str, help="SAM file")
parser.add_option("-o", "--output", dest="output_name", type=str, help="output prefix")
parser.add_option("--min", dest="min_bp", type=int, default=1,help="minimum fragment length, default = 1")
parser.add_option("--max", dest="max_bp", type=int, default=500,help="maximum fragment length, default = 500")
parser.add_option("-r", "--region", dest="region", type=str, help="region to be analyzed(eg.-r chrIII:12000:13000)")
parser.add_option("-n", "--normalized", dest="normalized", action="store_true", help="set to get normalized counts")
(options, args) = parser.parse_args()

sam_name = options.sam_name
output_name = options.output_name
min_bp = options.min_bp
max_bp = options.max_bp
region = options.region
normalized = options.normalized

##Filter(length and position)
if region is None:
    sam = open(sam_name, 'r')
    sam_readlines = sam.readlines()
    size = open(output_name + "_size_temp", 'w')
    for r in range(0, len(sam_readlines)):
        if "@" in sam_readlines[r].split()[0]:
            pass
        else:
            if min_bp <= int(sam_readlines[r].split()[8]) < max_bp:
                size.write(f'{int(sam_readlines[r].split()[8])}\n')
    sam.close()
    size.close()
else:
    chr_name = region.split(":")[0]
    start = int(region.split(":")[1])
    end = int(region.split(":")[2])
    sam = open(sam_name, 'r')
    sam_readlines = sam.readlines()
    size = open(output_name + "_size_temp", 'w')
    for r in range(0, len(sam_readlines)):
        if "@" in sam_readlines[r].split()[0]:
            pass
        else:
            if sam_readlines[r].split()[2] == chr_name:
                if start <= (int(sam_readlines[r].split()[3]) + int(int(sam_readlines[r].split()[8]) / 2)) < end:
                    if min_bp <= int(sam_readlines[r].split()[8]) < max_bp:
                        size.write(f'{int(sam_readlines[r].split()[8])}\n')
    sam.close()
    size.close()

##Counting
count = open(output_name + "_count.txt", 'w')
c = collections.Counter()
with open(output_name + "_size_temp") as st:
    for text in st:
        c.update([(text.strip())])

c_sorted = sorted(c.most_common())

count.write(f'length\tcounts\n')
for key, val in c_sorted:
    count.write(f'{key}\t{val}\n')

count.close()
st.close()

##Normalized counts(area = 1 * bp range)
if normalized:
    count = open(output_name + "_count.txt", 'r')
    normalized_count = open(output_name + "_norm_count.txt", 'w')
    normalized_area = int(max_bp - min_bp)
    count_readlines = count.readlines()
    sum = 0
    for r in range(1, len(count_readlines)):
        sum += int(count_readlines[r].split()[1])
    normalization_factor = normalized_area / sum      #get normalization factor

    normalized_count.write(f'length\tcounts\n')
    for r in range(1, len(count_readlines)):
        normalized_count.write(f'{count_readlines[r].split()[0]}\t{float(count_readlines[r].split()[1])*float(normalization_factor)}\n')
    count.close()
    normalized_count.close()
    os.remove(output_name + "_count.txt")

##End
os.remove(output_name + "_size_temp")
print("Done!")