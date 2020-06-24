#author: siayouyang
from optparse import OptionParser

##OptionParser
parser = OptionParser()
parser.add_option("--sam1", dest="sam1", help="input data in sam format")
parser.add_option("--sam2", dest="sam2", help="input data in sam format")
parser.add_option("--sam3", dest="sam3", help="input data in sam format")
parser.add_option("-r", "--region", dest="region", type=str, help="region to be analyzed(eg.-r chrIII:12000:13000)")
parser.add_option("-o", "--output", dest="output_name", help="output prefix")
(options, args) = parser.parse_args()

region = options.region
sam1 = options.sam1
sam2 = options.sam2
sam3 = options.sam3
output_name = options.output_name

chr_name = region.split(":")[0]
start = int(region.split(":")[1])
end = int(region.split(":")[2])

sam1 = open(sam1, 'r')
sam1_readlines = sam1.readlines()
size_position = open(output_name + "_sizeposition.txt", 'w')
size_position.write(f'position\tsize\tsample\n')
for r in range(0, len(sam1_readlines)):
    if "@" in sam1_readlines[r].split()[0]:
        pass
    else:
        if sam1_readlines[r].split()[2] == chr_name:
            if start <= (int(sam1_readlines[r].split()[3]) + int(int(sam1_readlines[r].split()[8]) / 2)) < end:
                if 70 <= int(sam1_readlines[r].split()[8]) < 200:
                    size_position.write(f'{int(sam1_readlines[r].split()[3]) + int(int(sam1_readlines[r].split()[8]) / 2)}\t{sam1_readlines[r].split()[8]}\tS1\n')
sam1.close()
print(f"sam1 for {output_name} was completed!")

if sam2 is None:
    pass
else:
    sam2 = open(sam2, 'r')
    sam2_readlines = sam2.readlines()
    for r in range(0, len(sam2_readlines)):
        if "@" in sam2_readlines[r].split()[0]:
            pass
        else:
            if sam2_readlines[r].split()[2] == chr_name:
                if start <= (int(sam2_readlines[r].split()[3]) + int(int(sam2_readlines[r].split()[8]) / 2)) < end:
                    if 70 <= int(sam2_readlines[r].split()[8]) < 200:
                        size_position.write(f'{int(sam2_readlines[r].split()[3]) + int(int(sam2_readlines[r].split()[8]) / 2)}\t{sam2_readlines[r].split()[8]}\tS2\n')
    sam2.close()
    print(f"sam2 for {output_name} was completed!")

if sam3 is None:
    pass
else:
    sam3 = open(sam3, 'r')
    sam3_readlines = sam3.readlines()
    for r in range(0, len(sam3_readlines)):
        if "@" in sam3_readlines[r].split()[0]:
            pass
        else:
            if sam3_readlines[r].split()[2] == chr_name:
                if start <= (int(sam3_readlines[r].split()[3]) + int(int(sam3_readlines[r].split()[8]) / 2)) < end:
                    if 70 <= int(sam3_readlines[r].split()[8]) < 200:
                        size_position.write(f'{int(sam3_readlines[r].split()[3]) + int(int(sam3_readlines[r].split()[8]) / 2)}\t{sam3_readlines[r].split()[8]}\tS3\n')
    sam3.close()
    print(f"sam3 for {output_name} was completed!")

size_position.close()
