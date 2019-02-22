import matplotlib.pyplot as plt
from pyteomics import mass
import numpy as np

from spectrum_utils import spectrum, plot

# for i, spec in enumerate(reader.read_mgf('../data/iPRG2012.mgf')):
#     if i == 3809:
#         frags = list(fragments('AAGLTAAYAR', types='aby', maxcharge=2))
#         annotations = []
#         for m in spec.mz:
#             found = False
#             for frag in frags:
#                 if abs(frag[0] - m) <= 0.02:
#                     print(m, frag)
#                     annotations.append(frag[1])
#                     frags.remove(frag)
#                     found = True
#                     break
#             if not found:
#                 annotations.append(None)
#         print(len(annotations), len(spec.mz))
#         spec.annotation = np.asarray(annotations)

#         spec.peptide = 'AAGLTAAYAR'

#         fig, ax = plt.subplots(figsize=(14, 6))

#         spec.process_peaks(remove_precursor=True, max_peaks_used=50, min_intensity=0.01, scaling='sqrt')
#         plot.plot_spectrum(spec, ax, 100, 1400, 1)

#         plt.tight_layout()

#         plt.savefig('3_annotations.png')
#         # plt.show()
#         plt.close()

#         break


mz = [101.07122, 109.68925, 115.86999, 120.0811, 129.06595,
      129.1023, 130.06534, 130.08633, 130.95578, 131.11859,
      132.07968, 136.07567, 141.31848, 142.21034, 147.11273,
      155.08185, 159.07605, 159.09169, 160.07516, 160.09506,
      170.06009, 171.06319, 175.11871, 186.12335, 187.08667,
      197.12814, 199.10913, 199.17963, 200.10286, 201.10281,
      201.12341, 204.11336, 204.13414, 211.14496, 213.15979,
      214.09769, 214.1535, 215.10004, 216.64253, 217.13332,
      224.13977, 225.12349, 226.08267, 228.13486, 229.1201,
      233.16495, 239.08228, 242.14993, 243.10875, 244.10834,
      245.35168, 253.09628, 256.10834, 261.1596, 273.13486,
      282.17978, 284.10266, 285.10632, 287.17343, 289.11084,
      301.129, 302.1327, 305.1811, 310.1727, 313.18762,
      318.15506, 321.1563, 326.17813, 328.2328, 334.15466,
      344.00247, 344.19284, 345.1581, 347.17093, 352.19876,
      353.1817, 356.1928, 361.2311, 362.20374, 363.207,
      367.14435, 370.21045, 384.16763, 385.17123, 395.1348,
      402.18094, 409.20935, 412.16074, 413.16446, 417.17773,
      424.21753, 425.2271, 429.1876, 430.19342, 441.24667,
      454.23114, 462.21368, 463.1971, 480.2241, 481.22632,
      484.24506, 488.21954, 497.25027, 498.25397, 506.25055,
      508.22073, 512.25824, 516.2635, 517.2272, 521.2413,
      523.74426, 525.2426, 526.2479, 530.26794, 538.2638,
      542.2722, 545.24817, 548.2805, 549.2873, 550.7827,
      555.2859, 559.2943, 559.7948, 571.2866, 587.2942,
      599.28723, 599.7854, 600.4091, 608.29016, 629.78845,
      636.28156, 639.34784, 646.328, 647.32837, 653.3023,
      670.3162, 677.3326, 685.323, 692.32416, 695.3513,
      696.3513, 702.3403, 703.35065, 707.32135, 709.8466,
      710.34625, 717.88837, 718.3871, 718.88586, 724.3413,
      729.36005, 730.32825, 748.37524, 757.3715, 763.14374,
      766.2276, 766.38745, 767.39105, 767.7366, 783.813,
      831.4068, 835.4187, 858.4285, 859.41895, 869.44836,
      871.4104, 876.4352, 877.4197, 878.41345, 886.42194,
      887.41785, 894.44464, 895.44763, 904.4519, 905.4429,
      906.4541, 972.48206, 989.50586, 990.4984, 1007.53156,
      1008.5328, 1017.5188, 1018.51953, 1037.5978, 1046.4642,
      1082.4929, 1100.5747, 1117.5706, 1118.574, 1119.5591,
      1135.5819, 1136.5898, 1145.5717, 1214.5857, 1215.5576,
      1227.598, 1232.5897, 1233.6329, 1249.6306, 1250.6298,
      1260.6073, 1261.614, 1272.6572]
intensity = [81011.57, 4123.349, 4006.9321, 66933.17, 11000.844,
             67220.16, 88976.56, 24920.941, 4863.066, 5661.9194,
             41732.137, 31081.303, 4782.2046, 5040.6426, 49796.156,
             6270.5312, 18098.807, 546127.56, 6226.79, 29613.12,
             178521.7, 11150.161, 15695.941, 31083.436, 62828.98,
             41016.887, 5672.7583, 12625.682, 21674.656, 22450.273,
             17107.639, 15654.782, 120404.91, 11095.885, 11517.213,
             124366.77, 45346.207, 11787.429, 9688.401, 6285.3345,
             10114.264, 16933.383, 15288.553, 10701.1045, 18149.885,
             41805.56, 20613.857, 61143.707, 71694.67, 14007.289,
             11060.377, 17671.592, 18029.537, 35232.992, 19302.346,
             10412.318, 350351.7, 42339.85, 14224.52, 9759.676,
             222278.44, 15488.672, 116594.016, 15155.255, 26307.275,
             33678.273, 15527.332, 22253.271, 9714.508, 10548.278,
             5726.4375, 55536.63, 19797.115, 35993.074, 14344.74,
             78794.98, 38641.426, 17308.31, 321167., 49607.965,
             30500.693, 27519.66, 23826.633, 26812.957, 106583.16,
             11173.92, 15909.75, 146074.52, 23989.854, 19026.498,
             45752.113, 11420.22, 166889.69, 22282.277, 16894.553,
             17635.023, 18778.855, 40425.176, 75016.37, 15470.299,
             12459.969, 15727.637, 50514.742, 36844.938, 11248.006,
             18829.473, 16063.4, 13619.328, 13366.1045, 10708.629,
             10835.664, 52481.062, 19239.938, 31166.791, 12338.772,
             42192.125, 15621.709, 346269.72, 80146.69, 17259.057,
             14637.047, 66929.305, 34878.094, 26068.5, 12426.318,
             19593.424, 19109.434, 11150.9, 26131.781, 13620.475,
             21111.895, 10257.758, 16610.361, 12370.426, 16981.797,
             11065.478, 20208.62, 17041.19, 11709.368, 243265.67,
             71733.34, 10427.157, 13114.231, 13370.711, 17090.348,
             35983.938, 23507.3, 179288.53, 18139.242, 24420.232,
             11087.336, 12962.59, 25980.557, 18395.598, 10993.147,
             17204.465, 319705.3, 107382.64, 10592.768, 10974.094,
             14374.425, 14261.798, 18317.068, 25027.375, 10635.955,
             15698.649, 74991.85, 170194.28, 62863.05, 12419.671,
             12949.016, 240654.78, 115458.13, 13719.998, 13071.061,
             10587.502, 16986.84, 16629.197, 13235.243, 242701.9,
             136257.86, 12841.443, 11758.576, 11649.871, 11818.519,
             12165.414, 12222.698, 14005.287, 84136.805, 41772.547,
             132280.81, 81302.51, 11834.507, 12680., 12015.682,
             12064.93, 26659.176, 16421.662, 75185.07, 61431.11,
             22042.248, 18096.48, 12666.438]

retention_time = 44.649
precursor_mz = 718.360
precursor_charge = 2
identifier = '41840'

spec = spectrum.MsmsSpectrum(identifier, precursor_mz, precursor_charge,
                             mz, intensity, retention_time=retention_time,
                             peptide='WNQLQAFWGTGK')

fragment_tol_mass = 10
fragment_tol_mode = 'ppm'

# spec = (spec.set_mz_range(min_mz=100, max_mz=1400)
#         .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
#         .filter_intensity(min_intensity=0.05, max_num_peaks=150)
#         .scale_intensity(scaling='root')
#         .annotate_peaks(fragment_tol_mass, fragment_tol_mode,
#                         ion_types='aby'))

# fig, axes = plt.subplots(2, 2, True, True, figsize=(18, 12))
#
# # Panel 1: Raw spectrum.
# spec = spec.set_mz_range(min_mz=100, max_mz=1400)
#
# plot.spectrum(spec, ax=axes[0, 0])
#
# axes[0, 0].set_title('1) raw spectrum')
# axes[0, 0].xaxis.label.set_visible(False)
#
# # Panel 2: Remove noise peaks.
# spec = (spec.remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
#         .filter_intensity(min_intensity=0.05, max_num_peaks=150))
#
# plot.spectrum(spec, ax=axes[0, 1])
#
# axes[0, 1].set_title('2) remove noise peaks')
#
# axes[0, 1].xaxis.label.set_visible(False)
# axes[0, 1].yaxis.label.set_visible(False)
#
# # Panel 3: Scale peak intensities.
# spec = spec.scale_intensity(scaling='root')
#
# plot.spectrum(spec, ax=axes[1, 0])
#
# axes[1, 0].set_title('3) scale peak intensities')
#
# # Panel 4: Annotation fragment peaks.
# spec = spec.annotate_peaks(fragment_tol_mass, fragment_tol_mode,
#                            ion_types='aby')
#
# plot.spectrum(spec, ax=axes[1, 1])
#
# axes[1, 1].set_title('4) annotate peaks')
#
# axes[1, 1].yaxis.label.set_visible(False)
#
# fig.suptitle(f'{spec.peptide}/{spec.precursor_charge}',
#              y=.98, size='x-large', weight='bold')
#
# plt.tight_layout(pad=4, w_pad=1, h_pad=2)
#
# plt.savefig('spectrum_utils.png', dpi=100)
# # plt.show()
# plt.close()


fig, ax = plt.subplots(figsize=(12, 6))

import copy
spec0 = copy.deepcopy(spec)
spec = (spec.set_mz_range(min_mz=100, max_mz=1400)
        .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
        .filter_intensity(min_intensity=0.05, max_num_peaks=150)
        .scale_intensity(scaling='root')
        .annotate_peaks(fragment_tol_mass, fragment_tol_mode,
                        ion_types='aby'))

plot.mirror(spec0, spec, ax)

plt.tight_layout()

plt.savefig('temp.png', dpi=100)
# plt.show()
plt.close()
