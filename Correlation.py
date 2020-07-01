#author:siayouyang
from optparse import OptionParser
import os

##OptionParser
parser = OptionParser()
parser.add_option("--data1", dest="data1", help="data1")
parser.add_option("--data2", dest="data2", help="data2")
parser.add_option("-o", "--output", dest="output_name", help="output name")
(options, args) = parser.parse_args()

data1 = options.data1
data2 = options.data2
output_name = options.output_name

#sort
data1 = open(data1, 'r')
output_temp = open(output_name + "_cor_temp", 'w')
data1_readlines = data1.readlines()
output_temp.write(data1_readlines[0])
for n in range(1, 1000):
    int(n)
    for r in range(1, len(data1_readlines)):
        if int(data1_readlines[r].split()[0]) == n:
            output_temp.write(data1_readlines[r])
output_temp.close()

#add data2 counts to last column
data2 = open(data2, 'r')
data2_readlines = data2.readlines()
output_temp = open(output_name + "_cor_temp", 'r')
output_temp_readlines = output_temp.readlines()
output = open(output_name + "_cor.txt", 'w')

x1 = output_temp_readlines[0].strip("\n")
output.write(f'{x1}\t{output_name}\n')
for r1 in range(1, len(output_temp_readlines)):
    for r2 in range(1, len(data2_readlines)):
        if int(data2_readlines[r2].split()[0]) == int(output_temp_readlines[r1].split()[0]):
            x2 = output_temp_readlines[r1].strip("\n")
            output.write(f'{x2}\t{data2_readlines[r2].split()[1]}\n')
output.close()
data2.close()
output_temp.close()

#End
os.remove(output_name + "_cor_temp")
print("Done！")