FAIL_SV_PATHS = """
plots/svs/hifi.CHM13v2.200081.chr13_5643741_5644104_trsolve.png
plots/svs/hifi.CHM13v2.200081.chr20_9901_10876_TRF.png
plots/svs/hifi.CHM13v2.200082.chr15_20233035_20233346_trsolve.png
plots/svs/hifi.CHM13v2.200085.chr11_122042881_122048125_trsolve.png
plots/svs/hifi.CHM13v2.200086.chr9_149775419_149775650_trsolve.png
plots/svs/hifi.CHM13v2.200086.chrX_1315514_1316318_trsolve.png
plots/svs/hifi.CHM13v2.200087.chr11_62267107_62267219_trsolve.png
plots/svs/hifi.CHM13v2.200104.chr14_13859812_13859823_trsolve.png
plots/svs/hifi.CHM13v2.200106.chr12_132263404_132264200_trsolve.png
plots/svs/hifi.CHM13v2.200106.chr14_4163985_4172555_trsolve.png
plots/svs/hifi.CHM13v2.NA12879.chr6_61591864_61591891_TRF.png
plots/svs/hifi.CHM13v2.NA12879.chr12_10899728_10899759_trsolve.png
plots/svs/hifi.CHM13v2.NA12886.chr8_7527052_7527888_trsolve.png
plots/svs/hifi.CHM13v2.NA12886.chr16_16444137_16444777_trsolve.png
plots/svs/hifi.CHM13v2.NA12886.chr22_4496132_4498541_trsolve.png
plots/svs/hifi.CHM13v2.NA12886.chr22_16398030_16398052_trsolve.png
""".split()

FAIL_SVS = [t.split("/")[-1].split(".")[2] for t in FAIL_SV_PATHS]


FAIL_VNTR_PATHS = """
plots/vntrs/hifi.CHM13v2.200106.chr12_132263404_132264200_trsolve.png
plots/vntrs/hifi.CHM13v2.NA12879.chr6_61591864_61591891_TRF.png
plots/vntrs/hifi.CHM13v2.NA12884.chr20_51386936_51387345_trsolve.png
plots/vntrs/hifi.CHM13v2.NA12885.chr2_242399998_242404718_trsolve.png
""".split()

FAIL_VNTRS = [t.split("/")[-1].split(".")[2] for t in FAIL_VNTR_PATHS]

FAIL_RECURRENT_PATHS = """
fig/recurrent/CHM13v2.chr7_152705910_152706128_trsolve.png
fig/recurrent/CHM13v2.chr7_153042968_153043065_trsolve.png
fig/recurrent/CHM13v2.chr13_100612222_100612334_trsolve.png
fig/recurrent/CHM13v2.chr17_6906664_6906688_trsolve.png
fig/recurrent/CHM13v2.chr18_49776737_49777085_trsolve.png"""