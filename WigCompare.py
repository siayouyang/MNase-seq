#author: siayouyang
#substraction
import os
from optparse import OptionParser
from decimal import *
import math

#OptionParser
parser = OptionParser()
parser.add_option("-c", "--control", dest="control", help="control data in wiggle format")
parser.add_option("-t", "--treatment", dest="treatment", help="treatment data in wiggle format")
parser.add_option("-n", "--normalizeRegion", dest="normalize_region", help="region selected as a reference for normalization, <chromosome>:<start position>:<end position> (eg. chrIII:50:950)")
parser.add_option("--operation", dest="operation",default="substract",type=str, help="posible choices: substract, log2, ratio(default: \"substract\")")
parser.add_option("--chrI", dest="chrI", action="store_true", help="process chrI data")
parser.add_option("--chrII", dest="chrII", action="store_true", help="process chrII data")
parser.add_option("--chrIII", dest="chrIII", action="store_true", help="process chrIII data")
parser.add_option("--chrIV", dest="chrIV", action="store_true", help="process chrIV data")
parser.add_option("--chrV", dest="chrV", action="store_true", help="process chrV data")
parser.add_option("--chrVI", dest="chrVI", action="store_true", help="process chrVI data")
parser.add_option("--chrVII", dest="chrVII", action="store_true", help="process chrVII data")
parser.add_option("--chrVIII", dest="chrVIII", action="store_true", help="process chrVIII data")
parser.add_option("--chrIX", dest="chrIX", action="store_true", help="process chrIX data")
parser.add_option("--chrX", dest="chrX", action="store_true", help="process chrX data")
parser.add_option("--chrXI", dest="chrXI", action="store_true", help="process chrXI data")
parser.add_option("--chrXII", dest="chrXII", action="store_true", help="process chrXII data")
parser.add_option("--chrXIII", dest="chrXIII", action="store_true", help="process chrXIII data")
parser.add_option("--chrXIV", dest="chrXIV", action="store_true", help="process chrXIV data")
parser.add_option("--chrXV", dest="chrXV", action="store_true", help="process chrXV data")
parser.add_option("--chrXVI", dest="chrXVI", action="store_true", help="process chrXVI data")
parser.add_option("--chrM", dest="chrM", action="store_true", help="process chrM data")
(options, args) = parser.parse_args()

control = options.control
treatment = options.treatment
normalize_region = options.normalize_region
operation = options.operation
chrI = options.chrI
chrII = options.chrII
chrIII = options.chrIII
chrIV= options.chrIV
chrV = options.chrV
chrVI = options.chrVI
chrVII = options.chrVII
chrVIII = options.chrVIII
chrIX = options.chrIX
chrX = options.chrX
chrXI = options.chrXI
chrXII = options.chrXII
chrXIII = options.chrXIII
chrXIV= options.chrXIV
chrXV = options.chrXV
chrXVI = options.chrXVI
chrM = options.chrM

def WigNormalize(wig):
    #remove lines start with "#"
    wig_file_temp = open(wig, 'r')
    wig_file_temp_readlines = wig_file_temp.readlines()
    wig_file = open(wig + "_clean", 'w')

    for r in wig_file_temp_readlines:
        if r.startswith("#"):
            pass
        else:
            wig_file.write(r)
    wig_file_temp.close()
    wig_file.close()

    def split_chr_bp(_split_chr, chr):
        #saperate chromosome
        wig_file = open(wig + "_clean", 'r')
        wig_file_readlines = wig_file.readlines()
        wig_file_chr = open(wig + _split_chr, 'w')

        for a in range(0, len(wig_file_readlines)):
            if (wig_file_readlines[a]).split()[0] == chr:
                wig_file_chr.write(wig_file_readlines[a])
            else:
                pass

        wig_file.close()
        wig_file_chr.close()

        #split lines to 1bp
        wig_file = open(wig + _split_chr, 'r')
        wig_file_readlines = wig_file.readlines()
        wig_file = open(wig + _split_chr, 'w')
        for a in range(0, len(wig_file_readlines)):
            minus = int((wig_file_readlines[a]).split()[2]) - int((wig_file_readlines[a]).split()[1])
            b = int((wig_file_readlines[a]).split()[1])
            if minus > 1:
                for c in range(b, b + minus):
                    wig_file.write(f'{(wig_file_readlines[a]).split()[0]}\t{c}\t{c + 1}\t{(wig_file_readlines[a]).split()[3]}\n')
            elif minus == 1:
                wig_file.write(wig_file_readlines[a])
        wig_file.close()

    #Get normalization factor
    norm = normalize_region
    norm = (norm.upper().split(":"))
    norm_chr = norm[0]
    norm_chr_start = int(norm[1])
    norm_chr_end = int(norm[2])

    def get_normalization_factor(_split_chr, _split_chr_norm_factor_temp):
        wig_split_chr = open(wig + _split_chr, 'r')  ###
        wig_split_chr_readlines = wig_split_chr.readlines()
        wig_split_chr_norm_factor_temp = open(wig + _split_chr_norm_factor_temp, 'w')
        coverage_sum = 0
        for r in range(norm_chr_start, norm_chr_end):
            wig_split_chr_readlines_coverage = float((wig_split_chr_readlines[r]).split('\t')[3])
            coverage_sum += wig_split_chr_readlines_coverage
        coverage_mean = coverage_sum / (norm_chr_end - norm_chr_start)
        normalization_factor = 1 / coverage_mean
        wig_split_chr_norm_factor_temp.write(str(normalization_factor))
        wig_split_chr.close()

    #Normalization
    def normalization(_split_chr, _split_chr_norm, _split_chr_norm_factor_temp):
        wig_split_chr = open(wig + _split_chr, 'r')
        wig_split_chr_norm_factor_temp = open(wig + _split_chr_norm_factor_temp, 'r')
        wig_split_chr_readlines = wig_split_chr.readlines()
        wig_split_chr_norm_factor_temp_readlines = wig_split_chr_norm_factor_temp.readlines()
        normalization_factor = wig_split_chr_norm_factor_temp_readlines[0]
        wig_split_chr_norm = open(wig + _split_chr_norm, 'w')

        for r in range(0, len(wig_split_chr_readlines)):
            line_3_value = (wig_split_chr_readlines[r]).split()[3]
            line_3_value_normalized = Decimal(str(line_3_value)) * Decimal(str(normalization_factor))
            wig_split_chr_norm.write(
                f'{(wig_split_chr_readlines[r]).split()[0]}\t{(wig_split_chr_readlines[r]).split()[1]}\t{(wig_split_chr_readlines[r]).split()[2]}\t{line_3_value_normalized}\n')
        wig_split_chr_norm.close()

    # split and get normalization factor
    def split_and_getnorm(chr, _split_chr, chr_name, chr_name_UPPER, _split_chr_norm):
        if chr == True:
            split_chr_bp(_split_chr, chr_name)
            if norm_chr == chr_name_UPPER:
                get_normalization_factor(_split_chr, "_split_chr_norm_factor_temp")
                normalization(_split_chr, _split_chr_norm, "_split_chr_norm_factor_temp")
                os.remove(wig + _split_chr)
        elif norm_chr == chr_name_UPPER:
            split_chr_bp(_split_chr, chr_name)
            get_normalization_factor(_split_chr, "_split_chr_norm_factor_temp")
            os.remove(wig + _split_chr)

    split_and_getnorm(chrI, "_split_chrI", "chrI", "CHRI", "_split_chrI_norm")
    split_and_getnorm(chrII, "_split_chrII", "chrII", "CHRII", "_split_chrII_norm")
    split_and_getnorm(chrIII, "_split_chrIII", "chrIII", "CHRIII", "_split_chrIII_norm")
    split_and_getnorm(chrIV, "_split_chrIV", "chrIV", "CHRIV", "_split_chrIV_norm")
    split_and_getnorm(chrV, "_split_chrV", "chrV", "CHRV", "_split_chrV_norm")
    split_and_getnorm(chrVI, "_split_chrVI", "chrVI", "CHRVI", "_split_chrVI_norm")
    split_and_getnorm(chrVII, "_split_chrVII", "chrVII", "CHRVII", "_split_chrVII_norm")
    split_and_getnorm(chrVIII, "_split_chrVIII", "chrVIII", "CHRVIII", "_split_chrVIII_norm")
    split_and_getnorm(chrIX, "_split_chrIX", "chrIX", "CHRIX", "_split_chrIX_norm")
    split_and_getnorm(chrX, "_split_chrX", "chrX", "CHRX", "_split_chrX_norm")
    split_and_getnorm(chrXI, "_split_chrXI", "chrXI", "CHRXI", "_split_chrXI_norm")
    split_and_getnorm(chrXII, "_split_chrXII", "chrXII", "CHRXII", "_split_chrXII_norm")
    split_and_getnorm(chrXIII, "_split_chrXIII", "chrXIII", "CHRXIII", "_split_chrXIII_norm")
    split_and_getnorm(chrXIV, "_split_chrXIV", "chrXIV", "CHRXIV", "_split_chrXIV_norm")
    split_and_getnorm(chrXV, "_split_chrXV", "chrXV", "CHRXV", "_split_chrXV_norm")
    split_and_getnorm(chrXVI, "_split_chrXVI", "chrXVI", "CHRXVI", "_split_chrXVI_norm")
    split_and_getnorm(chrM, "_split_chrM", "chrM", "CHRM", "_split_chrM_norm")

    # normalize part
    def normalize_others(chr, _split_chr, chr_name_UPPER, _split_chr_norm):
        if norm_chr == chr_name_UPPER:
            pass
        elif chr == True:
            normalization(_split_chr, _split_chr_norm, "_split_chr_norm_factor_temp")
            os.remove(wig + _split_chr)

    normalize_others(chrI, "_split_chrI", "CHRI", "_split_chrI_norm")
    normalize_others(chrII, "_split_chrII", "CHRII", "_split_chrII_norm")
    normalize_others(chrIII, "_split_chrIII", "CHRIII", "_split_chrIII_norm")
    normalize_others(chrIV, "_split_chrIV", "CHRIV", "_split_chrIV_norm")
    normalize_others(chrV, "_split_chrV", "CHRV", "_split_chrV_norm")
    normalize_others(chrVI, "_split_chrVI", "CHRVI", "_split_chrVI_norm")
    normalize_others(chrVII, "_split_chrVII", "CHRVII", "_split_chrVII_norm")
    normalize_others(chrVIII, "_split_chrVIII", "CHRVIII", "_split_chrVIII_norm")
    normalize_others(chrIX, "_split_chrIX", "CHRIX", "_split_chrIX_norm")
    normalize_others(chrX, "_split_chrX", "CHRX", "_split_chrX_norm")
    normalize_others(chrXI, "_split_chrXI", "CHRXI", "_split_chrXI_norm")
    normalize_others(chrXII, "_split_chrXII", "CHRXII", "_split_chrXII_norm")
    normalize_others(chrXIII, "_split_chrXIII", "CHRXIII", "_split_chrXIII_norm")
    normalize_others(chrXIV, "_split_chrXIV", "CHRXIV", "_split_chrXIV_norm")
    normalize_others(chrXV, "_split_chrXV", "CHRXV", "_split_chrXV_norm")
    normalize_others(chrXVI, "_split_chrXVI", "CHRXVI", "_split_chrXVI_norm")
    normalize_others(chrM, "_split_chrM", "CHRM", "_split_chrM_norm")

    #End
    os.remove(wig + "_split_chr_norm_factor_temp")
    os.remove(wig + "_clean")

WigNormalize(control)
WigNormalize(treatment)

#Subtraction
def substraction(control_input, treatment_input, _split_chr_norm, _chr_norm_substract):
    control = open(control_input + _split_chr_norm, 'r')
    treatment = open(treatment_input + _split_chr_norm, 'r')
    control_readlines = control.readlines()
    treatment_readlines = treatment.readlines()
    change = open(treatment_input + _chr_norm_substract, 'w')

    for rc in range(0, len(control_readlines)):
        rc_value = (control_readlines[rc]).split()[3]
        rt_value = (treatment_readlines[rc]).split()[3]
        substract = Decimal(str(rt_value)) - Decimal(str(rc_value))   #substract
        change.write(f'{(control_readlines[rc]).split()[0]}\t{(control_readlines[rc]).split()[1]}\t{(control_readlines[rc]).split()[2]}\t{substract}\n')

    change.close()

#ratio
def ratio(control_input, treatment_input, _split_chr_norm, _chr_norm_ratio):
    control = open(control_input + _split_chr_norm, 'r')
    treatment = open(treatment_input + _split_chr_norm, 'r')
    control_readlines = control.readlines()
    treatment_readlines = treatment.readlines()
    change = open(treatment_input + _chr_norm_ratio, 'w')

    for rc in range(0, len(control_readlines)):
        rc_value = (control_readlines[rc]).split()[3]
        rt_value = (treatment_readlines[rc]).split()[3]
        if Decimal(str(rc_value)) == 0:
            change.write(f'{(control_readlines[rc]).split()[0]}\t{(control_readlines[rc]).split()[1]}\t{(control_readlines[rc]).split()[2]}\t{0}\n')
        elif Decimal(str(rt_value)) == 0:
            change.write(f'{(control_readlines[rc]).split()[0]}\t{(control_readlines[rc]).split()[1]}\t{(control_readlines[rc]).split()[2]}\t{0}\n')
        else:
            ratio = float(rt_value) / float(rc_value)  # log2
            change.write(
                f'{(control_readlines[rc]).split()[0]}\t{(control_readlines[rc]).split()[1]}\t{(control_readlines[rc]).split()[2]}\t{ratio}\n')

    change.close()

#log2
def log_2(control_input, treatment_input, _split_chr_norm, _chr_norm_log2):
    control = open(control_input + _split_chr_norm, 'r')
    treatment = open(treatment_input + _split_chr_norm, 'r')
    control_readlines = control.readlines()
    treatment_readlines = treatment.readlines()
    change = open(treatment_input + _chr_norm_log2, 'w')

    for rc in range(0, len(control_readlines)):
        rc_value = (control_readlines[rc]).split()[3]
        rt_value = (treatment_readlines[rc]).split()[3]
        if Decimal(str(rc_value)) == 0:
            change.write(f'{(control_readlines[rc]).split()[0]}\t{(control_readlines[rc]).split()[1]}\t{(control_readlines[rc]).split()[2]}\t{0}\n')
        elif Decimal(str(rt_value)) == 0:
            change.write(f'{(control_readlines[rc]).split()[0]}\t{(control_readlines[rc]).split()[1]}\t{(control_readlines[rc]).split()[2]}\t{0}\n')
        else:
            log_2 = math.log2(float(rt_value) / float(rc_value))   #log2
            change.write(f'{(control_readlines[rc]).split()[0]}\t{(control_readlines[rc]).split()[1]}\t{(control_readlines[rc]).split()[2]}\t{log_2}\n')

    change.close()


def choose_operation(chr, _split_chr_norm, _chr_norm_substract, _chr_norm_ratio, _chr_norm_log2):
    if chr == True:
        if operation == "substract":
            substraction(control, treatment, _split_chr_norm, _chr_norm_substract)
        if operation == "ratio":
            ratio(control, treatment, _split_chr_norm, _chr_norm_ratio)
        if operation == "log2":
            log_2(control, treatment, _split_chr_norm, _chr_norm_log2)


choose_operation(chrI, "_split_chrI_norm", "_chrI_norm_substract", "_chrI_norm_ratio", "_chrI_norm_log2")
choose_operation(chrII, "_split_chrII_norm", "_chrII_norm_substract", "_chrII_norm_ratio", "_chrII_norm_log2")
choose_operation(chrIII, "_split_chrIII_norm", "_chrIII_norm_substract", "_chrIII_norm_ratio", "_chrIII_norm_log2")
choose_operation(chrIV, "_split_chrIV_norm", "_chrIV_norm_substract", "_chrIV_norm_ratio", "_chrIV_norm_log2")
choose_operation(chrV, "_split_chrV_norm", "_chrV_norm_substract", "_chrV_norm_ratio", "_chrV_norm_log2")
choose_operation(chrVI, "_split_chrVI_norm", "_chrVI_norm_substract", "_chrVI_norm_ratio", "_chrVI_norm_log2")
choose_operation(chrVII, "_split_chrVII_norm", "_chrVII_norm_substract", "_chrVII_norm_ratio", "_chrVII_norm_log2")
choose_operation(chrVIII, "_split_chrVIII_norm", "_chrVIII_norm_substract", "_chrVIII_norm_ratio", "_chrVIII_norm_log2")
choose_operation(chrIX, "_split_chrIX_norm", "_chrIX_norm_substract", "_chrIX_norm_ratio", "_chrIX_norm_log2")
choose_operation(chrX, "_split_chrX_norm", "_chrX_norm_substract", "_chrX_norm_ratio", "_chrX_norm_log2")
choose_operation(chrXI, "_split_chrXI_norm", "_chrXI_norm_substract", "_chrXI_norm_ratio", "_chrXI_norm_log2")
choose_operation(chrXII, "_split_chrXII_norm", "_chrXII_norm_substract", "_chrXII_norm_ratio", "_chrXII_norm_log2")
choose_operation(chrXIII, "_split_chrXIII_norm", "_chrXIII_norm_substract", "_chrXIII_norm_ratio", "_chrXIII_norm_log2")
choose_operation(chrXIV, "_split_chrXIV_norm", "_chrXIV_norm_substract", "_chrXIV_norm_ratio", "_chrXIV_norm_log2")
choose_operation(chrXV, "_split_chrXV_norm", "_chrXV_norm_substract", "_chrXV_norm_ratio", "_chrXV_norm_log2")
choose_operation(chrXVI, "_split_chrXVI_norm", "_chrXVI_norm_substract", "_chrXVI_norm_ratio", "_chrXVI_norm_log2")
choose_operation(chrX, "_split_chrM_norm", "_chrM_norm_substract", "_chrM_norm_ratio", "_chrM_norm_log2")
print("Done!")




