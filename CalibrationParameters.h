// TIGRESS/GRIFFIN Ge gains

float non_lin[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2.48315e-07, -3.01499e-07, 2.03059e-06, -4.18081e-07, -1.38592e-07, -9.33839e-08, 8.72549e-08, -1.27197e-07, 3.87659e-08, -6.8497e-07, -2.68466e-07, -1.17992e-07, -1.42218e-07, 8.3554e-08, -1.05125e-06, 1.03612e-07, -1.81968e-07, -1.39806e-07, -1.76273e-07, -2.88925e-07, 1.48833e-07, -7.91312e-08, 3.899e-09, 3.83256e-07, 3.64447e-08, 3.66722e-07, -4.9484e-07, -1.78657e-07, -3.16504e-07, 1.03819e-07, -5.31851e-07, -2.95469e-07, -2.91015e-07, -3.83884e-07, -5.89389e-08, 1.85899e-07, 6.32101e-07, 1.47406e-07, -7.14032e-07, 1.92781e-08, -1.01788e-06, -1.87902e-07, 6.93737e-08, 1.22437e-07, -5.18553e-08, -2.09315e-07, -7.8508e-07, -4.57289e-07};
float gain[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.29824, 1.28741, 1.68592, 1.27918, 1.29036, 1.25617, 1.2825, 1.25592, 1.36305, 1.36864, 1.33405, 0.936814, 1.29927, 1.44866, 1.34019, 1.28923, 1.3323, 1.32125, 1.08959, 1.33709, 1.26544, 1.26551, 1.29949, 1.29991, 1.28803, 1.30102, 1.34206, 1.29085, 1.37187, 1.37281, 1.40174, 1.36961, 1.28566, 1.3141, 1.34688, 1.30675, 1.31599, 1.31037, 1.31262, 1.34386, 1.30952, 1.44078, 1.26994, 1.34796, 1.3621, 1.34535, 1.39637, 1.381};
float offset[64] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1.68195, -0.373984, -10.1209, -0.124188, -0.501908, -0.381316, 0.216181, -0.443385, -0.0245366, -0.360332, -0.830114, -0.527173, -0.924022, -0.158051, -2.16199, -0.36441, -0.339923, -1.10286, -0.175503, 0.0256084, 0.47484, -0.309727, -0.363141, -0.259711, -0.00393496, 0.22273, -0.643067, 0.230376, -0.252369, 0.292646, 0.139656, -0.331143, 0.0742457, -0.34732, -0.0443097, 0.445289, -0.272722, 0.0895158, -0.494853, -0.620237, -0.704072, -0.144047, 0.0661586, 0.175775, -0.742415, -0.820456, -0.509444, -0.608925};
float lin_gain[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.29778, 1.28688, 1.68873, 1.27839, 1.29011, 1.256, 1.28265, 1.25567, 1.36312, 1.36744, 1.33359, 0.936523, 1.299, 1.44879, 1.3384, 1.28942, 1.33197, 1.32101, 1.08926, 1.33658, 1.26572, 1.26537, 1.2995, 1.30059, 1.28809, 1.30165, 1.34122, 1.29052, 1.37134, 1.37298, 1.40087, 1.36912, 1.28516, 1.31346, 1.34678, 1.30705, 1.31709, 1.31062, 1.31139, 1.34389, 1.30776, 1.44047, 1.27006, 1.34816, 1.36201, 1.34499, 1.39511, 1.38025};
float lin_offset[64] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1.58504, -0.277754, -10.553, 0.0405208, -0.444167, -0.351358, 0.189869, -0.383314, -0.0401338, -0.129269, -0.747447, -0.443627, -0.870148, -0.177725, -1.8409, -0.404259, -0.273148, -1.05867, -0.109757, 0.130907, 0.414738, -0.284276, -0.364484, -0.390078, -0.016522, 0.110038, -0.499589, 0.296988, -0.158211, 0.262137, 0.284905, -0.241258, 0.161487, -0.239478, -0.0267516, 0.394407, -0.465941, 0.0469726, -0.283757, -0.625441, -0.397669, -0.0891372, 0.0448686, 0.141838, -0.72778, -0.756607, -0.303541, -0.48681};

/*float non_lin[64];
float gain[64] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
float offset[64];
float lin_gain[64];
float lin_offset[64];*/

// GRIFFIN Only
float paces_gain[5];
float paces_offset[5];
float paces_nonlin[5];

float LaBr_gain[8];
float LaBr_offset[8];
float LaBr_nonlin[8];

// TIGRESS Only
float s3gains[56] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
float s3offsets[56];

/*float s3gains[56] = {3.73688, 3.97644, 3.97209, 4.06567, 4.07417, 3.72729, 3.99868, 4.14024, 3.65777, 3.70051, 3.8939, 4.01821, 4.11233, 3.83754, 4.2669, 3.91142, 
3.74409, 3.90966, 3.95283, 3.9848, 3.97579, 3.91388, 3.83781, 3.68854, 3.57907, 3.79474, 3.80197, 3.73692, 3.6547, 3.6464, 3.6785, 3.96741, 
  3.91081, 3.91085, 3.94589, 3.93481, 4.18267, 3.84952, 3.87708, 3.80978, 3.78623, 3.93487, 3.92981, 3.84301, 3.73069, 3.83246, 3.81863, 3.82885, 3.81034, 4.1219, 3.84844, 3.83402, 3.83547, 3.81999, 3.87006, 3.93667};
float s3offsets[56] = {163.919, 164.166, 160.275, 135.236, 153.593, 144.202, 159.505, 156.324, 153.355, 158.61, 166.344, 150.418, 158.264, 151.029, 148.582, 144.393,
146.607, 162.96, 151.189, 160.89, 163.497, 170.852, 158.729, 155.921, 155.509, 187.151, 174.887, 157.223, 168.924, 179.212, 164.941, 179.936, 
 193.402, 179.428, 198.154, 202.369, 202.388, 176.442, 186.69, 121.445, 171.03, 208.899, 160.162, 166.975, 281.534, 161.815, 164.099, 87.0373, 143.557, 41.5523, 145.487, 110.628, 130.831, 129.202, 139.787, 119.528};*/

float sharcgains[960];
float sharcoffsets[960];

float csigains[128] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
float csioffsets[128] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

float trificgains[64];
float trificoffsets[64];

float seggains[512]; 
float segoffsets[512];
