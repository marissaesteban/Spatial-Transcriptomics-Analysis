f = open("marissa_analysis/SNP_analysis/CoMutPlotter/all_snps_tsv.tsv")
w = open("marissa_analysis/SNP_analysis/CoMutPlotter/all_snps_tsv2.tsv", "w")

for line in f:
    line = line.strip().split()

    for item in line:
        w.write(item.strip() + "\t")

    w.write("\n")

f.close()
w.close()
