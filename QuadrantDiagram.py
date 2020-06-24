#author: siayouyang
import os
from optparse import OptionParser

##OptionParser
parser = OptionParser()
parser.add_option("-c", "--chr", dest="chr_name", type=str, help="chromosome to be analyzed(eg.-c chrIII)")
parser.add_option("-s", "--sub", dest="subnucleosome", help="subnucleosome data in wiggle format")
parser.add_option("-n", "--nuc", dest="nucleosome", help="nucleosome data in wiggle format")
parser.add_option("-p", "--pos", dest="position", help="position data in xls format from DANPOS")
parser.add_option("-o", "--out", dest="output_name", type=str, help="a name for output file(.xlsx)")
parser.add_option("--ResectionStart", dest="resection_start", type=int, default=0, help="resection region start position (defualt = 0)")
parser.add_option("--ResectionEnd", dest="resection_end", type=int, default=0, help="resection region end position (default = 0)")
parser.add_option("--ArbitraryStart", dest="arbitrary_start", type=int, default=0, help="arbitrary region start position (defualt = 0)")
parser.add_option("--ArbitraryEnd", dest="arbitrary_end", type=int, default=0, help="arbitrary region end position (default = 0)")
(options, args) = parser.parse_args()

chr_name = options.chr_name
subnucleosome = options.subnucleosome
nucleosome = options.nucleosome
position = options.position
output_name = options.output_name
resection_start = options.resection_start
resection_end = options.resection_end
arbitrary_start = options.arbitrary_start
arbitrary_end = options.arbitrary_end

#Extract data
nucleosome = open(nucleosome, 'r')
subnucleosome = open(subnucleosome, 'r')
nucleosome_readlines = nucleosome.readlines()
subnucleosome_readlines = subnucleosome.readlines()
nucleosome_occupancy = open(output_name + "_nucleosome_occ_temp", 'w')
subnucleosome_occupancy = open(output_name + "_subnucleosome_occ_temp", 'w')
position = open(position, 'r')
position_readlines = position.readlines()
position_summit = open(output_name + "_position_summit_temp", 'w')

for rp in range(0, len(position_readlines)):
    if ((position_readlines[rp]).split()[0]).upper() == chr_name.upper():
        rp_value = (position_readlines[rp]).split()[3]
        position_summit.write(f'{rp_value}\n')
position_summit.close()


for rn in range(0, len(nucleosome_readlines)):
    rn_value = (nucleosome_readlines[rn]).split()[3]
    nucleosome_occupancy.write(f'{rn_value}\n')
nucleosome_occupancy.close()

for rs in range(0, len(subnucleosome_readlines)):
    rs_value = (subnucleosome_readlines[rs]).split()[3]
    subnucleosome_occupancy.write(f'{rs_value}\n')
subnucleosome_occupancy.close()

#get nucleosome and subnucleosome occupancy mean change
nucleosome_occupancy = open(output_name + "_nucleosome_occ_temp", 'r')
nucleosome_occupancy_readlines = nucleosome_occupancy.readlines()
subnucleosome_occupancy = open(output_name + "_subnucleosome_occ_temp", 'r')
subnucleosome_occupancy_readlines = subnucleosome_occupancy.readlines()
position_summit = open(output_name + "_position_summit_temp", 'r')
position_summit_readlines = position_summit.readlines()
nucleosome_mean_change = open(output_name + "_nucleosome_mean_temp", 'w')
for pos_r in range(0, len(position_summit_readlines)):
    pos_n = int(position_summit_readlines[pos_r])
    nuc_sum = 0
    for nuc_r in range((pos_n - 1)-73, pos_n + 73):
        nuc_sum += float(nucleosome_occupancy_readlines[nuc_r])
    nuc_mean = nuc_sum/147
    nucleosome_mean_change.write(f'{nuc_mean}\n')
nucleosome_mean_change.close()
nucleosome_occupancy.close()

subnucleosome_mean_change = open(output_name + "_subnucleosome_mean_temp", 'w')
for pos_r in range(0, len(position_summit_readlines)):
    pos_n = int(position_summit_readlines[pos_r])
    sub_sum = 0
    for sub_r in range((pos_n - 1)-55, pos_n + 55):
        sub_sum += float(subnucleosome_occupancy_readlines[sub_r])
    sub_mean = sub_sum/111
    subnucleosome_mean_change.write(f'{sub_mean}\n')
subnucleosome_mean_change.close()
subnucleosome_occupancy.close()

#Combine lists
nucleosome_mean_change = open(output_name + "_nucleosome_mean_temp", 'r')
nucleosome_mean_change_readlines = nucleosome_mean_change.readlines()
subnucleosome_mean_change = open(output_name + "_subnucleosome_mean_temp", 'r')
subnucleosome_mean_change_readlines = subnucleosome_mean_change.readlines()
full_list = open(output_name + "_chr.txt", 'w')
full_list.write(f'position\tnucleosome\tsubnucleosome\n')
for pos_r in range(1, len(position_summit_readlines)):  #exclude the first nucleosome occupancy data
    full_list.write(f'{position_summit_readlines[pos_r].split()[0]}\t{nucleosome_mean_change_readlines[pos_r].split()[0]}\t{subnucleosome_mean_change_readlines[pos_r].split()[0]}\n')
full_list.close()
nucleosome_mean_change.close()
subnucleosome_mean_change.close()
position_summit.close()

#classification
full_list = open(output_name + "_chr.txt", 'r')
full_list_readlines = full_list.readlines()
full_list_class = open(output_name + "_chr_class.txt", 'w')
II_chr = 0
ID_chr = 0
DI_chr = 0
DD_chr = 0
for r in range(1, len(full_list_readlines)):                #exclude title row
    r_nuc = float(full_list_readlines[r].split()[1])
    r_sub = float(full_list_readlines[r].split()[2])
    if r_nuc >= 0:
        if r_sub >= 0:
            II_chr += 1
        elif r_sub < 0:
            ID_chr += 1
    elif r_nuc < 0:
        if r_sub >= 0:
            DI_chr += 1
        elif r_sub < 0:
            DD_chr += 1
IN_chr = ID_chr + II_chr
DN_chr = DD_chr + DI_chr
IS_chr = II_chr + DI_chr
DS_chr = ID_chr + DD_chr
total_chr = II_chr + ID_chr + DI_chr + DD_chr
II_percent_chr = (II_chr / total_chr) * 100
ID_percent_chr = (ID_chr / total_chr) * 100
DI_percent_chr = (DI_chr / total_chr) * 100
DD_percent_chr = (DD_chr / total_chr) * 100
IN_percent_chr = (IN_chr / total_chr) * 100
DN_percent_chr = (DN_chr / total_chr) * 100
IS_percent_chr = (IS_chr / total_chr) * 100
DS_percent_chr = (DS_chr / total_chr) * 100
full_list_class.write(f'Total counts: {total_chr}\n')
full_list_class.write(f'Type\tCounts\t% Percentage\n')
full_list_class.write(f'II\t{II_chr}\t{II_percent_chr}\n')
full_list_class.write(f'ID\t{ID_chr}\t{ID_percent_chr}\n')
full_list_class.write(f'DI\t{DI_chr}\t{DI_percent_chr}\n')
full_list_class.write(f'DD\t{DD_chr}\t{DD_percent_chr}\n')
full_list_class.write(f'IN\t{IN_chr}\t{IN_percent_chr}\n')
full_list_class.write(f'DN\t{DN_chr}\t{DN_percent_chr}\n')
full_list_class.write(f'IS\t{IS_chr}\t{IS_percent_chr}\n')
full_list_class.write(f'DS\t{DS_chr}\t{DS_percent_chr}\n')
full_list_class.write(f'abbreviations:\nII = increased nucleosome and increased subnucleosome;\nDI = decreased nucleosome and increased subnucleosome;\nID = increased nucleosome and decreased subnucleosome;\nDD = decreased nucleosome and decreased subnucleosome;\nIN = increased nucleosome;\nDN = decreased nucleosome;\nIS = increased subnucleosome;\nDS = decreased subnucleosome')
full_list_class.close()

#get resection region
if resection_start == 0 and resection_end == 0:
    print('no resection region was defined')
else:
    full_list = open(output_name + "_chr.txt", 'r')
    full_list_readlines = full_list.readlines()
    resection_list = open(output_name + "_resection.txt", 'w')
    resection_list.write(f'position\tnucleosome\tsubnucleosome\n')   #add title
    for r in range(1, len(full_list_readlines)):                #exclude title row
        r_position = int(full_list_readlines[r].split()[0])
        if resection_start <= r_position < resection_end:
            resection_list.write(full_list_readlines[r])
    resection_list.close()

    # classification(resection)
    resection_list = open(output_name + "_resection.txt", 'r')
    resection_list_readlines = resection_list.readlines()
    resection_list_class = open(output_name + "_resection_class.txt", 'w')
    II_chr = 0
    ID_chr = 0
    DI_chr = 0
    DD_chr = 0
    for r in range(1, len(resection_list_readlines)):  # exclude title row
        r_nuc = float(resection_list_readlines[r].split()[1])
        r_sub = float(resection_list_readlines[r].split()[2])
        if r_nuc >= 0:
            if r_sub >= 0:
                II_chr += 1
            elif r_sub < 0:
                ID_chr += 1
        elif r_nuc < 0:
            if r_sub >= 0:
                DI_chr += 1
            elif r_sub < 0:
                DD_chr += 1
    IN_chr = ID_chr + II_chr
    DN_chr = DD_chr + DI_chr
    IS_chr = II_chr + DI_chr
    DS_chr = ID_chr + DD_chr
    total_chr = II_chr + ID_chr + DI_chr + DD_chr
    II_percent_chr = (II_chr / total_chr) * 100
    ID_percent_chr = (ID_chr / total_chr) * 100
    DI_percent_chr = (DI_chr / total_chr) * 100
    DD_percent_chr = (DD_chr / total_chr) * 100
    IN_percent_chr = (IN_chr / total_chr) * 100
    DN_percent_chr = (DN_chr / total_chr) * 100
    IS_percent_chr = (IS_chr / total_chr) * 100
    DS_percent_chr = (DS_chr / total_chr) * 100
    resection_list_class.write(f'Total counts: {total_chr}\n')
    resection_list_class.write(f'Type\tCounts\tPercentage\n')
    resection_list_class.write(f'II\t{II_chr}\t{II_percent_chr}\n')
    resection_list_class.write(f'ID\t{ID_chr}\t{ID_percent_chr}\n')
    resection_list_class.write(f'DI\t{DI_chr}\t{DI_percent_chr}\n')
    resection_list_class.write(f'DD\t{DD_chr}\t{DD_percent_chr}\n')
    resection_list_class.write(f'IN\t{IN_chr}\t{IN_percent_chr}\n')
    resection_list_class.write(f'DN\t{DN_chr}\t{DN_percent_chr}\n')
    resection_list_class.write(f'IS\t{IS_chr}\t{IS_percent_chr}\n')
    resection_list_class.write(f'DS\t{DS_chr}\t{DS_percent_chr}\n')
    resection_list_class.write(
        f'abbreviations:\nII = increased nucleosome and increased subnucleosome;\nDI = decreased nucleosome and increased subnucleosome;\nID = increased nucleosome and decreased subnucleosome;\nDD = decreased nucleosome and decreased subnucleosome;\nIN = increased nucleosome;\nDN = decreased nucleosome;\nIS = increased subnucleosome;\nDS = decreased subnucleosome')
    resection_list_class.close()

#get arbitraty region
if arbitrary_start == 0 and arbitrary_end == 0:
    print('no arbitrary region was defined')
else:
    full_list = open(output_name + "_chr.txt", 'r')
    full_list_readlines = full_list.readlines()
    arbitrary_list = open(output_name + "_arbitrary.txt", 'w')
    arbitrary_list.write(f'position\tnucleosome\tsubnucleosome\n')   #add title
    for r in range(1, len(full_list_readlines)):                #exclude title row
        r_position = int(full_list_readlines[r].split()[0])
        if arbitrary_start <= r_position < arbitrary_end:
            arbitrary_list.write(full_list_readlines[r])
    arbitrary_list.close()

    # classification(arbitrary)
    arbitrary_list = open(output_name + "_arbitrary.txt", 'r')
    arbitrary_list_readlines = arbitrary_list.readlines()
    arbitrary_list_class = open(output_name + "_arbitrary_class.txt", 'w')
    II_chr = 0
    ID_chr = 0
    DI_chr = 0
    DD_chr = 0
    for r in range(1, len(arbitrary_list_readlines)):  # exclude title row
        r_nuc = float(arbitrary_list_readlines[r].split()[1])
        r_sub = float(arbitrary_list_readlines[r].split()[2])
        if r_nuc >= 0:
            if r_sub >= 0:
                II_chr += 1
            elif r_sub < 0:
                ID_chr += 1
        elif r_nuc < 0:
            if r_sub >= 0:
                DI_chr += 1
            elif r_sub < 0:
                DD_chr += 1
    IN_chr = ID_chr + II_chr
    DN_chr = DD_chr + DI_chr
    IS_chr = II_chr + DI_chr
    DS_chr = ID_chr + DD_chr
    total_chr = II_chr + ID_chr + DI_chr + DD_chr
    II_percent_chr = (II_chr / total_chr) * 100
    ID_percent_chr = (ID_chr / total_chr) * 100
    DI_percent_chr = (DI_chr / total_chr) * 100
    DD_percent_chr = (DD_chr / total_chr) * 100
    IN_percent_chr = (IN_chr / total_chr) * 100
    DN_percent_chr = (DN_chr / total_chr) * 100
    IS_percent_chr = (IS_chr / total_chr) * 100
    DS_percent_chr = (DS_chr / total_chr) * 100
    arbitrary_list_class.write(f'Total counts: {total_chr}\n')
    arbitrary_list_class.write(f'Type\tCounts\t% Percentage\n')
    arbitrary_list_class.write(f'II\t{II_chr}\t{II_percent_chr}\n')
    arbitrary_list_class.write(f'ID\t{ID_chr}\t{ID_percent_chr}\n')
    arbitrary_list_class.write(f'DI\t{DI_chr}\t{DI_percent_chr}\n')
    arbitrary_list_class.write(f'DD\t{DD_chr}\t{DD_percent_chr}\n')
    arbitrary_list_class.write(f'IN\t{IN_chr}\t{IN_percent_chr}\n')
    arbitrary_list_class.write(f'DN\t{DN_chr}\t{DN_percent_chr}\n')
    arbitrary_list_class.write(f'IS\t{IS_chr}\t{IS_percent_chr}\n')
    arbitrary_list_class.write(f'DS\t{DS_chr}\t{DS_percent_chr}\n')
    arbitrary_list_class.write(
        f'abbreviations:\nII = increased nucleosome and increased subnucleosome;\nDI = decreased nucleosome and increased subnucleosome;\nID = increased nucleosome and decreased subnucleosome;\nDD = decreased nucleosome and decreased subnucleosome;\nIN = increased nucleosome;\nDN = decreased nucleosome;\nIS = increased subnucleosome;\nDS = decreased subnucleosome')
    arbitrary_list_class.close()

#End
os.remove(output_name + "_nucleosome_occ_temp")
os.remove(output_name + "_subnucleosome_occ_temp")
os.remove(output_name + "_position_summit_temp")
os.remove(output_name + "_subnucleosome_mean_temp")
os.remove(output_name + "_nucleosome_mean_temp")

