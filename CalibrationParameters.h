// TIGRESS/GRIFFIN Ge gains

/*float non_lin[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
5.08228e-07, -4.82052e-08, -6.18976e-07, -2.04252e-09, 
-1.52255e-07, -1.43147e-07, 3.84638e-07, 1.20819e-07, 
3.69298e-07, -4.24956e-07, -2.6627e-07, -3.3745e-07, 
-3.26023e-08, -5.07603e-07, -9.8375e-08, -4.97298e-07, 
1.6349e-07, -6.46231e-07, -1.56677e-08, 2.87455e-07, 
6.43196e-07, 9.42003e-07, 3.41351e-07, 5.45562e-07, 
6.26241e-08, -2.00187e-07, -1.9018e-07, -4.67897e-07, 
4.64305e-07, 5.77788e-07, -1.08002e-07, 8.95037e-08, 
-6.49338e-08, -1.60127e-07, -1.52535e-07, -1.67949e-07, 
8.39588e-08, 4.16906e-07, -5.8567e-07, -3.37886e-08, 
-3.83614e-07, 1.70321e-07, 1.61107e-07, 1.69194e-07, 
-4.5703e-07, 3.30654e-08, -1.87088e-07, 7.88521e-09};
float gain[64] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1.30028, 1.28837, 1.70327, 1.28174, 
1.29706, 1.28156, 1.29066, 1.25439, 
1.27338, 1.31997, 1.33142, 1.28226, 
1.34632, 1.38626, 1.28489, 1.2551, 
1.3171, 1.321, 1.08626, 1.33325, 
1.35087, 1.32508, 1.3511, 1.35688, 
1.29207, 1.2862, 1.30334, 1.28084, 
1.35352, 1.36695, 1.40272, 1.3716, 
1.28648, 1.30738, 1.33873, 1.30532, 
1.34277, 1.35464, 1.31176, 0.949038, 
1.29298, 1.31377, 1.31107, 1.33904, 
1.31943, 1.34217, 1.37831, 1.36135};
float offset[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0.0129006, 0.187672, -1.86355, -1.17198, 
-0.337478, -0.305201, -0.392721, -0.17677, 
-0.139695, -0.560315, -0.807924, -0.366275, 
-0.156494, -0.822417, -0.26079, -0.787941, 
-0.675975, -1.29244, 0.25322, 0.0820229, 
-0.417743, -0.0229732, 0.55548, -0.345156, 
0.16319, -0.401742, -1.87201, -1.28556, 
0.160531, -0.0789081, -0.138985, 0.0575451, 
-0.223433, -0.158921, -0.0803983, 0.204968, 
-0.52123, 0.450782, 0.0110789, -0.413933, 
-0.883553, -0.41119, -0.0226218, 0.524527, 
-0.202664, -0.66735, -0.0894732, -0.805515};
float lin_gain[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.30087, 1.28841, 1.70288, 1.2817, 1.29692, 1.28145, 1.29113, 1.25456, 1.27397, 1.31934, 1.33102, 1.28174, 1.34627, 1.38552, 1.28474, 1.25436, 1.31735, 1.32005, 1.08623, 1.33365, 1.35183, 1.32634, 1.35159, 1.35765, 1.29228, 1.28601, 1.30318, 1.28033, 1.35403, 1.36762, 1.40267, 1.37169, 1.28647, 1.3069, 1.33859, 1.3052, 1.34289, 1.35525, 1.31087, 0.948974, 1.29238, 1.31402, 1.31131, 1.33928, 1.31897, 1.34221, 1.37812, 1.36136};
float lin_offset[64] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -0.0756639, 0.162384, -1.86127, -1.15846, -0.322693, -0.301253, -0.466999, -0.205858, -0.246014, -0.459248, -0.739909, -0.275278, -0.149026, -0.70629, -0.23651, -0.666433, -0.718591, -1.14038, 0.258446, 0.0196571, -0.574329, -0.207601, 0.478532, -0.462971, 0.10018, -0.382297, -1.86441, -1.21967, 0.0903841, -0.180591, -0.146822, 0.0457484, -0.234643, -0.0211939, -0.0671971, 0.208166, -0.538588, 0.355154, 0.167279, -0.400373, -0.782271, -0.45052, -0.0605921, 0.488574, -0.150351, -0.673065, -0.0651614, -0.80525};
*/

/*float non_lin[64] = {-4.33985e-13, -4.11572e-13, -1.91511e-13, -3.46902e-13, 
-1.03391e-13, -6.40734e-15, -4.21315e-12, -9.00197e-13, 
-4.69342e-14, 3.28456e-13, 2.40462e-14, -2.67286e-13, 
-2.27582e-14, -3.13272e-14, 1.47449e-13, -1.91356e-13, 
7.40235e-12, -2.60834e-14, -7.60202e-13, -1.21616e-13, 
0, 0, 0, 0, 
-1.23306e-12, -6.97267e-13, -1.09428e-12, -1.58905e-13, 
6.44738e-13, 2.37066e-14, -8.78668e-13, -3.94523e-13, 
4.34075e-12, 1.63151e-12, 3.85991e-12, 3.2549e-12, 
0, 0, 0, 0, 
-5.12215e-13, -7.75826e-13, -7.83725e-14, 4.95669e-13, 
-6.53137e-14, -3.87778e-13, -2.00704e-11, -6.50996e-13, 
2.65887e-13, -1.90465e-13, -4.67978e-13, -3.88648e-13, 
-4.27582e-13, -4.04759e-13, -6.74639e-14, -2.71155e-13, 
2.94753e-13, -5.59241e-13, 1.68098e-12, -1.65704e-13, 
-4.8818e-13, 3.14235e-14, -8.40148e-14, -9.90666e-13};

float gain[64] = {0.00124693, 0.00122323, 0.00128181, 0.00125681, 
0.00128276, 0.00130098, 0.00131408, 0.00132802, 
0.00129451, 0.00126057, 0.00132559, 0.00125164, 
0.00126012, 0.00127729, 0.00130037, 0.00123209, 
0.00130083, 0.00124948, 0.00138401, 0.00122456, 
0, 0, 0, 0, 
0.00131385, 0.00129657, 0.00132362, 0.00132808, 
0.00131462, 0.00134047, 0.00122668, 0.00123061, 
0.00122869, 0.00128896, 0.00132073, 0.00127069, 
0, 0, 0, 0, 
0.00127967, 0.00130118, 0.00128329, 0.00131038, 
0.00128321, 0.0012717, 0.00128311, 0.00129538, 
0.00128228, 0.00127936, 0.0010416, 0.00126186, 
0.0013154, 0.00132028, 0.00128045, 0.00125726, 
0.00133492, 0.00124204, 0.00134666, 0.0013307, 
0.0012475, 0.00130637, 0.00129872, 0.00129605};

float offset[64] = {0.12889, -0.0990848, 0.0312276, -0.156058, 
0.0644573, 0.0833964, -1.41268, -1.31676, 
-0.584787, 0.130665, -0.666296, -0.487789, 
0.172625, 0.131745, -0.172239, -0.246176, 
3.29367, -0.177151, -0.341942, 0.0223437, 
-1, -1, -1, -1, 
-0.429573, -0.235645, 0.295308, -0.142769, 
0.0234017, -0.375087, -0.42788, -0.475948, 
1.06516, -0.171322, 0.868352, 0.721381, 
-1, -1, -1, -1, 
-0.110568, -0.256732, -0.121521, 0.418946, 
-0.371162, -0.260258, 0.377676, -0.392372, 
0.426418, 0.0329988, -0.0502978, -0.213185, 
-0.121009, -0.0411156, -0.489989, -0.298866, 
-0.132055, -0.160205, 0.570791, 0.0801563, 
-0.0959226, -0.245665, -0.28175, -0.513453};*/

float non_lin[64];
float gain[64] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
float offset[64];
float lin_gain[64];
float lin_offset[64];

// GRIFFIN Only
float paces_gain[5];
float paces_offset[5];
float paces_nonlin[5];

float LaBr_gain[8];
float LaBr_offset[8];
float LaBr_nonlin[8];

// TIGRESS Only
float s3gains[56] = {3.63649, 3.91071, 3.9746, 4.00467, 4.02152, 3.92327, 3.84077, 3.70667, 3.61808, 3.86967, 3.84195, 3.73979, 3.70377, 3.65657, 3.66666, 3.99898, 3.76519, 3.98185, 4.01095, 4.05842, 4.07346, 3.68905, 3.9904, 4.1232, 3.72379, 3.75035, 3.94548, 4.16266, 4.10576, 3.87673, 4.31536, 3.92705, 3.90813, 3.9533, 4.00154, 3.89349, 4.15728, 3.85721, 3.89569, 3.73031, 3.78473, 3.9041, 3.87616, 3.86734, 3.88288, 3.8784, 3.74162, 3.80638, 3.9263, 4.18367, 3.90189, 3.84411, 3.79734, 3.79303, 3.80986, 3.86041};
float s3offsets[56] = {206.203, 139.417, 160.275, 154.121, 145.603, 152.77, 150.4, 133.286, 135.012, 138.753, 140.563, 120.378, 138.224, 169.963, 136.235, 129.698, 142.864, 132.031, 154.276, 157.382, 162.381, 175.551, 169.201, 159.442, 150.83, 187.947, 182.371, 172.099, 164.152, 126.482, 124.755, 175.988, 202.047, 160.313, 188.739, 214.331, 181.777, 160.879, 155.104, 177.557, 140.598, 163.294, 158.101, 146.121, 150.661, 126.482, 162.801, 139.258, 142.417, 130.851, 108.929, 131.826, 124.937, 141.03, 149.27, 122.84};

float sharcgains[960];
float sharcoffsets[960];

float csigains[128] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
float csioffsets[128] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

float trificgains[64];
float trificoffsets[64];

float seggains[512]; 
float segoffsets[512];
