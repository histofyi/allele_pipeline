


motif_alleles = ["hla_a_01_01", "hla_a_02_01", "hla_a_02_02", "hla_a_02_03", "hla_a_02_04", "hla_a_02_05", "hla_a_02_06", "hla_a_02_07", "hla_a_02_11", "hla_a_02_20", "hla_a_02_52", "hla_a_03_01", "hla_a_03_02", "hla_a_11_01", "hla_a_11_02", "hla_a_23_01", "hla_a_24_02", "hla_a_24_07", "hla_a_25_01", "hla_a_26_01", "hla_a_26_08", "hla_a_29_02", "hla_a_30_01", "hla_a_30_02", "hla_a_31_01", "hla_a_32_01", "hla_a_33_01", "hla_a_33_03", "hla_a_34_01", "hla_a_34_02", "hla_a_36_01", "hla_a_66_01", "hla_a_68_01", "hla_a_68_02", "hla_a_69_01", "hla_a_74_01", "hla_b_07_02", "hla_b_07_04", "hla_b_08_01", "hla_b_13_01", "hla_b_13_02", "hla_b_14_01", "hla_b_14_02", "hla_b_15_01", "hla_b_15_02", "hla_b_15_03", "hla_b_15_10", "hla_b_15_11", "hla_b_15_13", "hla_b_15_17", "hla_b_15_18", "hla_b_18_01", "hla_b_18_03", "hla_b_18_05", "hla_b_27_04", "hla_b_27_05", "hla_b_27_09", "hla_b_35_01", "hla_b_35_02", "hla_b_35_03", "hla_b_35_07", "hla_b_35_08", "hla_b_37_01", "hla_b_38_01", "hla_b_38_02", "hla_b_39_01", "hla_b_39_05", "hla_b_39_06", "hla_b_39_24", "hla_b_40_01", "hla_b_40_02", "hla_b_40_06", "hla_b_40_32", "hla_b_41_01", "hla_b_42_01", "hla_b_44_02", "hla_b_44_03", "hla_b_44_05", "hla_b_45_01", "hla_b_46_01", "hla_b_47_01", "hla_b_48_01", "hla_b_49_01", "hla_b_50_01", "hla_b_51_01", "hla_b_51_08", "hla_b_52_01", "hla_b_53_01", "hla_b_54_01", "hla_b_55_01", "hla_b_55_02", "hla_b_56_01", "hla_b_57_01", "hla_b_57_03", "hla_b_58_01", "hla_b_58_02", "hla_b_67_01", "hla_b_73_01", "hla_b_81_01", "hla_c_01_02", "hla_c_02_02", "hla_c_03_02", "hla_c_03_03", "hla_c_03_04", "hla_c_04_01", "hla_c_04_03", "hla_c_05_01", "hla_c_06_02", "hla_c_07_01", "hla_c_07_02", "hla_c_07_04", "hla_c_08_01", "hla_c_08_02", "hla_c_12_02", "hla_c_12_03", "hla_c_12_04", "hla_c_14_02", "hla_c_14_03", "hla_c_15_02", "hla_c_15_05", "hla_c_16_01", "hla_c_16_02", "hla_c_17_01", "hla_e_01_03", "hla_g_01_01", "hla_g_01_03", "hla_g_01_04"]
structure_alleles = ["hla_a_01_01", "hla_a_02_01", "hla_a_02_03", "hla_a_02_06", "hla_a_02_07", "hla_a_03_01", "hla_a_11_01", "hla_a_24_02", "hla_a_24_50", "hla_a_29_02", "hla_a_30_01", "hla_a_30_03", "hla_a_68_01", "hla_a_68_02", "hla_b_07_02", "hla_b_08_01", "hla_b_14_02", "hla_b_15_01", "hla_b_15_02", "hla_b_18_01", "hla_b_27_03", "hla_b_27_05", "hla_b_27_06", "hla_b_27_09", "hla_b_27_247", "hla_b_35_01", "hla_b_35_08", "hla_b_37_01", "hla_b_39_01", "hla_b_40_01", "hla_b_40_02", "hla_b_41_04", "hla_b_42_01", "hla_b_42_02", "hla_b_44_02", "hla_b_44_03", "hla_b_44_05", "hla_b_46_01", "hla_b_51_01", "hla_b_52_01", "hla_b_53_01", "hla_b_57_01", "hla_b_57_03", "hla_b_57_06", "hla_b_58_01", "hla_b_58_11", "hla_b_81_01", "hla_c_01_02", "hla_c_03_04", "hla_c_04_01", "hla_c_05_01", "hla_c_06_02", "hla_c_07_02", "hla_c_08_01", "hla_c_08_02", "hla_c_14_02", "hla_e_01_01", "hla_e_01_03", "hla_g_01_01"]


motif_set = set(motif_alleles)
structure_set = set(structure_alleles)
print (f"Number of motifs alleles: {len(motif_set)}")
print (f"Number of structure alleles: {len(structure_set)}")

print (f"Alleles for which we have structures and motifs: {len(motif_set.intersection(structure_set))}")
print (f"Alleles for which we have motifs, but no structures: {len(motif_set.difference(structure_set))}")
print (f"Alleles for which we have structures, but no motifs: {len(structure_set.difference(motif_set))}")

from matplotlib_venn import venn2, venn2_circles
import matplotlib.pyplot as plt

venn2(subsets = (len(motif_set.difference(structure_set)), len(structure_set.difference(motif_set)), len(motif_set.intersection(structure_set))), set_labels = ('Motifs', 'Structures'))
plt.savefig('motifs_vs_structures_venn.svg')
plt.savefig('motifs_vs_structures_venn.png')