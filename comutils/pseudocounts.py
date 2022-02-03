import numpy as np
import ccmpred.pseudocounts

# P_BLOSUM[a,b] = P[b|a]
# computed from CCMPredPy
P_BLOSUM = np.array([[0.29014845, 0.03103914, 0.02564103, 0.02968961, 0.02159244,
        0.02564103, 0.04048583, 0.0782726 , 0.0148448 , 0.04318489,
        0.05937922, 0.04453441, 0.01754386, 0.02159244, 0.02968961,
        0.08502024, 0.04993252, 0.00539811, 0.01754386, 0.06882591],
       [0.04457364, 0.34496124, 0.03875969, 0.03100775, 0.00775194,
        0.04844961, 0.05232558, 0.03294574, 0.02325581, 0.02325581,
        0.04651163, 0.12015504, 0.01550388, 0.01744186, 0.01937984,
        0.04457364, 0.03488372, 0.00581395, 0.01744186, 0.03100775],
       [0.04269663, 0.04494382, 0.31685393, 0.08314607, 0.00898876,
        0.03370787, 0.0494382 , 0.06516854, 0.03146067, 0.02247191,
        0.03146067, 0.05393258, 0.01123596, 0.01797753, 0.02022472,
        0.06966292, 0.0494382 , 0.00449438, 0.01573034, 0.02696629],
       [0.04104478, 0.02985075, 0.06902985, 0.39738806, 0.00746269,
        0.02985075, 0.09141791, 0.04664179, 0.01865672, 0.02238806,
        0.02798507, 0.04477612, 0.00932836, 0.01492537, 0.02238806,
        0.05223881, 0.03544776, 0.00373134, 0.01119403, 0.02425373],
       [0.06504065, 0.01626016, 0.01626016, 0.01626016, 0.48373984,
        0.01219512, 0.01626016, 0.03252033, 0.00813008, 0.04471545,
        0.06504065, 0.0203252 , 0.01626016, 0.0203252 , 0.01626016,
        0.04065041, 0.03658537, 0.00406504, 0.01219512, 0.05691057],
       [0.05588235, 0.07352941, 0.04411765, 0.04705882, 0.00882353,
        0.21470588, 0.10294118, 0.04117647, 0.02941176, 0.02647059,
        0.04705882, 0.09117647, 0.02058824, 0.01470588, 0.02352941,
        0.05588235, 0.04117647, 0.00588235, 0.02058824, 0.03529412],
       [0.05524862, 0.04972376, 0.04051565, 0.09023941, 0.00736648,
        0.06445672, 0.29650092, 0.03499079, 0.02578269, 0.02209945,
        0.03683241, 0.07550645, 0.01289134, 0.01657459, 0.02578269,
        0.05524862, 0.03683241, 0.00552486, 0.01657459, 0.03130755],
       [0.0782726 , 0.02294197, 0.0391363 , 0.03373819, 0.01079622,
        0.01889339, 0.02564103, 0.51012146, 0.01349528, 0.01889339,
        0.02834008, 0.03373819, 0.00944669, 0.01619433, 0.01889339,
        0.05128205, 0.02968961, 0.00539811, 0.01079622, 0.0242915 ],
       [0.04198473, 0.04580153, 0.05343511, 0.03816794, 0.00763359,
        0.03816794, 0.05343511, 0.03816794, 0.35496183, 0.02290076,
        0.03816794, 0.04580153, 0.01526718, 0.03053435, 0.01908397,
        0.04198473, 0.02671756, 0.00763359, 0.05725191, 0.02290076],
       [0.04712813, 0.01767305, 0.01472754, 0.01767305, 0.01620029,
        0.01325479, 0.01767305, 0.02061856, 0.00883652, 0.27098675,
        0.16789396, 0.02356406, 0.03681885, 0.04418262, 0.01472754,
        0.02503682, 0.03976436, 0.00589102, 0.02061856, 0.17673049],
       [0.04453441, 0.0242915 , 0.01417004, 0.01518219, 0.01619433,
        0.01619433, 0.02024291, 0.02125506, 0.01012146, 0.11538462,
        0.37550607, 0.02530364, 0.04959514, 0.05465587, 0.01417004,
        0.0242915 , 0.03340081, 0.00708502, 0.02226721, 0.09615385],
       [0.05699482, 0.10708117, 0.04145078, 0.04145078, 0.00863558,
        0.05354059, 0.07081174, 0.04317789, 0.02072539, 0.02763385,
        0.04317789, 0.27806563, 0.01554404, 0.01554404, 0.02763385,
        0.05354059, 0.03972366, 0.00518135, 0.01727116, 0.0328152 ],
       [0.05220884, 0.03212851, 0.02008032, 0.02008032, 0.01606426,
        0.02811245, 0.02811245, 0.02811245, 0.01606426, 0.10040161,
        0.19678715, 0.03614458, 0.16064257, 0.04819277, 0.01606426,
        0.03614458, 0.04016064, 0.00803213, 0.02409639, 0.09236948],
       [0.03382664, 0.01902748, 0.01691332, 0.01691332, 0.01057082,
        0.01057082, 0.01902748, 0.02536998, 0.01691332, 0.06342495,
        0.1141649 , 0.01902748, 0.02536998, 0.38689218, 0.01057082,
        0.02536998, 0.02536998, 0.01691332, 0.08879493, 0.05496829],
       [0.05684755, 0.02583979, 0.02325581, 0.03100775, 0.01033592,
        0.02067183, 0.03617571, 0.03617571, 0.0129199 , 0.02583979,
        0.03617571, 0.04134367, 0.01033592, 0.0129199 , 0.49354005,
        0.04392765, 0.03617571, 0.00258398, 0.0129199 , 0.03100775],
       [0.10994764, 0.04013962, 0.05410122, 0.04886562, 0.01745201,
        0.03315881, 0.05235602, 0.06631763, 0.01919721, 0.02966841,
        0.04188482, 0.05410122, 0.01570681, 0.02094241, 0.02966841,
        0.21989529, 0.08202443, 0.0052356 , 0.01745201, 0.04188482],
       [0.0729783 , 0.03550296, 0.0433925 , 0.03747535, 0.01775148,
        0.02761341, 0.03944773, 0.0433925 , 0.01380671, 0.05325444,
        0.06508876, 0.04536489, 0.01972387, 0.02366864, 0.02761341,
        0.09270217, 0.24654832, 0.00591716, 0.01775148, 0.07100592],
       [0.03030303, 0.02272727, 0.01515152, 0.01515152, 0.00757576,
        0.01515152, 0.02272727, 0.03030303, 0.01515152, 0.03030303,
        0.0530303 , 0.02272727, 0.01515152, 0.06060606, 0.00757576,
        0.02272727, 0.02272727, 0.49242424, 0.06818182, 0.03030303],
       [0.04049844, 0.02803738, 0.02180685, 0.01869159, 0.00934579,
        0.02180685, 0.02803738, 0.02492212, 0.04672897, 0.04361371,
        0.06853583, 0.03115265, 0.01869159, 0.13084112, 0.01557632,
        0.03115265, 0.02803738, 0.02803738, 0.31775701, 0.04672897],
       [0.06995885, 0.02194787, 0.01646091, 0.01783265, 0.01920439,
        0.01646091, 0.02331962, 0.02469136, 0.00823045, 0.16460905,
        0.1303155 , 0.0260631 , 0.03155007, 0.03566529, 0.01646091,
        0.03292181, 0.04938272, 0.00548697, 0.02057613, 0.26886145]])



def get_cond_proba(a, knowing):
    return P_BLOSUM[knowing, a]
