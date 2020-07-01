# author:siayouyang
from optparse import OptionParser

##OptionParser
parser = OptionParser()
parser.add_option("-d", "--data", dest="data_name", type=str, help="length-counts data")
(options, args) = parser.parse_args()

data_name = options.data_name

data = open(data_name, 'r')
data_readlines = data.readlines()
sum = 0
for r in range(1, len(data_readlines)):
    sum += float(data_readlines[r].split()[0]) * float(data_readlines[r].split()[1])
print(sum / float(len(data_readlines) - 1))
data.close()