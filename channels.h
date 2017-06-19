#ifndef CHANNELS_H
#define CHANNELS_H
#include <initializer_list>
#include <vector>
#include <string>
#include <cmath>

#ifndef UTILITIES_H
#include "utilities.h"
#endif

using namespace std;

// ----- typedefs

// symplectic_map is a function func(double *in, double *out, int nqubits)
typedef int (*symplectic_map)(double *, double *, int);

// ----- Adjoint representation of noisy T channels

// T + depolarizing noise
int noisyT_1Q(double p, double *y) {
	y[1] = 1 - (4.*p)/3.;
	y[2] = 0;
	y[3] = 0; 
	y[4] = 0; 
	y[5] = (1 - p)/sqrt(2) - p/(3.*sqrt(2));
	y[6] = -((1 - p)/sqrt(2)) + p/(3.*sqrt(2));
	y[7] = 0;
	y[8] = (1 - p)/sqrt(2) - p/(3.*sqrt(2));
	y[9] = (1 - p)/sqrt(2) - p/(3.*sqrt(2));
}

// T + random Pauli noise
int noisyT_1Q(double *p, double *y) {
	y[1] = p[1] - p[2] - p[3] + p[4];
	y[2] = y[3] = y[4] = 0;
	y[5] = p[1]/sqrt(2) + p[2]/sqrt(2) - p[3]/sqrt(2) - p[4]/sqrt(2);
	y[6] = -(p[1]/sqrt(2)) + p[2]/sqrt(2) - p[3]/sqrt(2) + p[4]/sqrt(2);
	y[7] = 0;
	y[8] = p[1]/sqrt(2) + p[2]/sqrt(2) - p[3]/sqrt(2) - p[4]/sqrt(2);
	y[9] = p[1]/sqrt(2) - p[2]/sqrt(2) + p[3]/sqrt(2) - p[4]/sqrt(2);
}

// T on first qubit + (depolarizing noise) otimes (depolarizing noise)
int noisyT_2Q(double p, double *y) {
	vector<double> ytemp = {pow(1 - p,2) - pow(p,2)/9.,0,0,0,0,(-2*(1 - p)*p)/(3.*sqrt(3)) + (2*pow(p,2))/(9.*sqrt(3)),
   0,0,0,0,(2*sqrt(0.6666666666666666)*(1 - p)*p)/3. - (2*sqrt(0.6666666666666666)*pow(p,2))/9.,0,0,
   0,0,0,pow(1 - p,2) - pow(p,2)/9.,0,0,0,0,0,0,0,0,0,(2*(1 - p)*p)/3. - (2*pow(p,2))/9.,0,0,0,0,
   0,pow(1 - p,2)/sqrt(2) - pow(p,2)/(9.*sqrt(2)),0,0,0,0,
   (sqrt(2)*(1 - p)*p)/3. - (sqrt(2)*pow(p,2))/9.,
   -(pow(1 - p,2)/sqrt(2)) + pow(p,2)/(9.*sqrt(2)),0,0,0,0,
   -(sqrt(2)*(1 - p)*p)/3. + (sqrt(2)*pow(p,2))/9.,0,0,0,0,
   pow(1 - p,2)/sqrt(2) - (sqrt(2)*(1 - p)*p)/3. + pow(p,2)/(9.*sqrt(2)),0,0,0,0,0,0,0,0,
   -(pow(1 - p,2)/sqrt(2)) + (sqrt(2)*(1 - p)*p)/3. - pow(p,2)/(9.*sqrt(2)),0,0,0,0,0,0,
   pow(1 - p,2) - pow(p,2)/9.,0,0,0,0,0,0,0,0,0,(2*(1 - p)*p)/3. - (2*pow(p,2))/9.,
   (-2*(1 - p)*p)/(3.*sqrt(3)) + (2*pow(p,2))/(9.*sqrt(3)),0,0,0,0,
   pow(1 - p,2) + (4*(1 - p)*p)/9. - (7*pow(p,2))/27.,0,0,0,0,
   (2*sqrt(2)*(1 - p)*p)/9. - (2*sqrt(2)*pow(p,2))/27.,0,0,0,0,0,0,0,0,0,0,
   pow(1 - p,2)/sqrt(2) - (sqrt(2)*(1 - p)*p)/3. + pow(p,2)/(9.*sqrt(2)),0,0,
   -(pow(1 - p,2)/sqrt(2)) + (sqrt(2)*(1 - p)*p)/3. - pow(p,2)/(9.*sqrt(2)),0,0,0,0,0,0,0,
   (sqrt(2)*(1 - p)*p)/3. - (sqrt(2)*pow(p,2))/9.,0,0,0,0,
   pow(1 - p,2)/sqrt(2) - pow(p,2)/(9.*sqrt(2)),-(sqrt(2)*(1 - p)*p)/3. + (sqrt(2)*pow(p,2))/9.,
   0,0,0,0,-(pow(1 - p,2)/sqrt(2)) + pow(p,2)/(9.*sqrt(2)),0,0,0,
   pow(1 - p,2)/sqrt(2) - pow(p,2)/(9.*sqrt(2)),0,0,0,0,
   (sqrt(2)*(1 - p)*p)/3. - (sqrt(2)*pow(p,2))/9.,pow(1 - p,2)/sqrt(2) - pow(p,2)/(9.*sqrt(2)),0,
   0,0,0,(sqrt(2)*(1 - p)*p)/3. - (sqrt(2)*pow(p,2))/9.,0,0,0,0,0,0,0,
   pow(1 - p,2)/sqrt(2) - (sqrt(2)*(1 - p)*p)/3. + pow(p,2)/(9.*sqrt(2)),0,0,
   pow(1 - p,2)/sqrt(2) - (sqrt(2)*(1 - p)*p)/3. + pow(p,2)/(9.*sqrt(2)),0,0,0,0,0,
   (2*sqrt(0.6666666666666666)*(1 - p)*p)/3. - (2*sqrt(0.6666666666666666)*pow(p,2))/9.,0,0,0,0,
   (2*sqrt(2)*(1 - p)*p)/9. - (2*sqrt(2)*pow(p,2))/27.,0,0,0,0,
   pow(1 - p,2) + (2*(1 - p)*p)/9. - (5*pow(p,2))/27.,0,0,0,0,0,
   (2*(1 - p)*p)/3. - (2*pow(p,2))/9.,0,0,0,0,0,0,0,0,0,pow(1 - p,2) - pow(p,2)/9.,0,0,0,0,0,0,
   pow(1 - p,2)/sqrt(2) - (sqrt(2)*(1 - p)*p)/3. + pow(p,2)/(9.*sqrt(2)),0,0,0,0,0,0,0,0,
   pow(1 - p,2)/sqrt(2) - (sqrt(2)*(1 - p)*p)/3. + pow(p,2)/(9.*sqrt(2)),0,0,0,0,
   (sqrt(2)*(1 - p)*p)/3. - (sqrt(2)*pow(p,2))/9.,0,0,0,0,
   pow(1 - p,2)/sqrt(2) - pow(p,2)/(9.*sqrt(2)),(sqrt(2)*(1 - p)*p)/3. - (sqrt(2)*pow(p,2))/9.,0,
   0,0,0,pow(1 - p,2)/sqrt(2) - pow(p,2)/(9.*sqrt(2)),0,0,0,0,0,
   (2*(1 - p)*p)/3. - (2*pow(p,2))/9.,0,0,0,0,0,0,0,0,0,pow(1 - p,2) - pow(p,2)/9.};

   copy(ytemp.begin(), ytemp.end(), y+1);

   return 0;
}

// T on first qubit + random Pauli noise on 2 qubits
int noisyT_2Q(double *p, double *y) {

  y[0] = 0;
  y[1] = p[0] - p[1] + p[12] - p[13] - p[14] + p[15] - p[2] + p[3];
  y[2] = 0;
  y[3] = 0;
  y[4] = 0;
  y[5] = 0;
  y[6] = p[10]/sqrt(3) - p[11]/sqrt(3) - p[4]/sqrt(3) + p[5]/sqrt(3) + p[6]/sqrt(3) - p[7]/sqrt(3) - p[8]/sqrt(3) + p[9]/sqrt(3);
  y[7] = 0;
  y[8] = 0;
  y[9] = 0;
  y[10] = 0;
  y[11] = -(sqrt(0.6666666666666666)*p[10]) + sqrt(0.6666666666666666)*p[11] + sqrt(0.6666666666666666)*p[4] - sqrt(0.6666666666666666)*p[5] - sqrt(0.6666666666666666)*p[6] + sqrt(0.6666666666666666)*p[7] + sqrt(0.6666666666666666)*p[8] - sqrt(0.6666666666666666)*p[9];
  y[12] = 0;
  y[13] = 0;
  y[14] = 0;
  y[15] = 0;
  y[16] = 0;
  y[17] = p[0] + p[1] + p[12] + p[13] - p[14] - p[15] - p[2] - p[3];
  y[18] = 0;
  y[19] = 0;
  y[20] = 0;
  y[21] = 0;
  y[22] = 0;
  y[23] = 0;
  y[24] = 0;
  y[25] = 0;
  y[26] = 0;
  y[27] = -p[10] - p[11] + p[4] + p[5] - p[6] - p[7] + p[8] + p[9];
  y[28] = 0;
  y[29] = 0;
  y[30] = 0;
  y[31] = 0;
  y[32] = 0;
  y[33] = p[0]/sqrt(2) - p[11]/sqrt(2) - p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) + p[4]/sqrt(2) + p[7]/sqrt(2) - p[8]/sqrt(2);
  y[34] = 0;
  y[35] = 0;
  y[36] = 0;
  y[37] = 0;
  y[38] = p[1]/sqrt(2) - p[10]/sqrt(2) - p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) + p[5]/sqrt(2) + p[6]/sqrt(2) - p[9]/sqrt(2);
  y[39] = -(p[0]/sqrt(2)) - p[11]/sqrt(2) + p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) + p[4]/sqrt(2) + p[7]/sqrt(2) - p[8]/sqrt(2);
  y[40] = 0;
  y[41] = 0;
  y[42] = 0;
  y[43] = 0;
  y[44] = -(p[1]/sqrt(2)) - p[10]/sqrt(2) + p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) + p[5]/sqrt(2) + p[6]/sqrt(2) - p[9]/sqrt(2);
  y[45] = 0;
  y[46] = 0;
  y[47] = 0;
  y[48] = 0;
  y[49] = p[0]/sqrt(2) + p[10]/sqrt(2) - p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) + p[5]/sqrt(2) - p[6]/sqrt(2) - p[9]/sqrt(2);
  y[50] = 0;
  y[51] = 0;
  y[52] = p[1]/sqrt(2) + p[11]/sqrt(2) - p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) + p[4]/sqrt(2) - p[7]/sqrt(2) - p[8]/sqrt(2);
  y[53] = 0;
  y[54] = 0;
  y[55] = -(p[1]/sqrt(2)) + p[11]/sqrt(2) + p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) + p[4]/sqrt(2) - p[7]/sqrt(2) - p[8]/sqrt(2);
  y[56] = 0;
  y[57] = 0;
  y[58] = -(p[0]/sqrt(2)) + p[10]/sqrt(2) + p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) + p[5]/sqrt(2) - p[6]/sqrt(2) - p[9]/sqrt(2);
  y[59] = 0;
  y[60] = 0;
  y[61] = 0;
  y[62] = 0;
  y[63] = 0;
  y[64] = 0;
  y[65] = p[0] - p[1] + p[12] - p[13] + p[14] - p[15] + p[2] - p[3];
  y[66] = 0;
  y[67] = 0;
  y[68] = 0;
  y[69] = 0;
  y[70] = 0;
  y[71] = 0;
  y[72] = 0;
  y[73] = 0;
  y[74] = 0;
  y[75] = p[10] - p[11] + p[4] - p[5] + p[6] - p[7] + p[8] - p[9];
  y[76] = p[10]/sqrt(3) - p[11]/sqrt(3) - p[4]/sqrt(3) + p[5]/sqrt(3) + p[6]/sqrt(3) - p[7]/sqrt(3) - p[8]/sqrt(3) + p[9]/sqrt(3);
  y[77] = 0;
  y[78] = 0;
  y[79] = 0;
  y[80] = 0;
  y[81] = p[0] + p[1]/3. - (2*p[10])/3. - (2*p[11])/3. + p[12] + p[13]/3. + p[14]/3. + p[15] + p[2]/3. + p[3] - (2*p[4])/3. - (2*p[5])/3. - (2*p[6])/3. - (2*p[7])/3. - (2*p[8])/3. - (2*p[9])/3.;
  y[82] = 0;
  y[83] = 0;
  y[84] = 0;
  y[85] = 0;
  y[86] = (2*sqrt(2)*p[1])/3. - (sqrt(2)*p[10])/3. - (sqrt(2)*p[11])/3. + (2*sqrt(2)*p[13])/3. + (2*sqrt(2)*p[14])/3. + (2*sqrt(2)*p[2])/3. - (sqrt(2)*p[4])/3. - (sqrt(2)*p[5])/3. - (sqrt(2)*p[6])/3. - (sqrt(2)*p[7])/3. - (sqrt(2)*p[8])/3. - (sqrt(2)*p[9])/3.;
  y[87] = 0;
  y[88] = 0;
  y[89] = 0;
  y[90] = 0;
  y[91] = 0;
  y[92] = 0;
  y[93] = 0;
  y[94] = p[1]/sqrt(2) + p[11]/sqrt(2) - p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) + p[4]/sqrt(2) - p[7]/sqrt(2) - p[8]/sqrt(2);
  y[95] = 0;
  y[96] = 0;
  y[97] = p[0]/sqrt(2) + p[10]/sqrt(2) - p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) + p[5]/sqrt(2) - p[6]/sqrt(2) - p[9]/sqrt(2);
  y[98] = 0;
  y[99] = 0;
  y[100] = -(p[0]/sqrt(2)) + p[10]/sqrt(2) + p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) + p[5]/sqrt(2) - p[6]/sqrt(2) - p[9]/sqrt(2);
  y[101] = 0;
  y[102] = 0;
  y[103] = -(p[1]/sqrt(2)) + p[11]/sqrt(2) + p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) + p[4]/sqrt(2) - p[7]/sqrt(2) - p[8]/sqrt(2);
  y[104] = 0;
  y[105] = 0;
  y[106] = 0;
  y[107] = 0;
  y[108] = p[1]/sqrt(2) - p[10]/sqrt(2) - p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) + p[5]/sqrt(2) + p[6]/sqrt(2) - p[9]/sqrt(2);
  y[109] = 0;
  y[110] = 0;
  y[111] = 0;
  y[112] = 0;
  y[113] = p[0]/sqrt(2) - p[11]/sqrt(2) - p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) + p[4]/sqrt(2) + p[7]/sqrt(2) - p[8]/sqrt(2);
  y[114] = -(p[1]/sqrt(2)) - p[10]/sqrt(2) + p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) + p[5]/sqrt(2) + p[6]/sqrt(2) - p[9]/sqrt(2);
  y[115] = 0;
  y[116] = 0;
  y[117] = 0;
  y[118] = 0;
  y[119] = -(p[0]/sqrt(2)) - p[11]/sqrt(2) + p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) + p[4]/sqrt(2) + p[7]/sqrt(2) - p[8]/sqrt(2);
  y[120] = 0;
  y[121] = 0;
  y[122] = 0;
  y[123] = p[0]/sqrt(2) - p[11]/sqrt(2) - p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) + p[4]/sqrt(2) + p[7]/sqrt(2) - p[8]/sqrt(2);
  y[124] = 0;
  y[125] = 0;
  y[126] = 0;
  y[127] = 0;
  y[128] = p[1]/sqrt(2) - p[10]/sqrt(2) - p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) + p[5]/sqrt(2) + p[6]/sqrt(2) - p[9]/sqrt(2);
  y[129] = p[0]/sqrt(2) + p[11]/sqrt(2) - p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) - p[4]/sqrt(2) - p[7]/sqrt(2) + p[8]/sqrt(2);
  y[130] = 0;
  y[131] = 0;
  y[132] = 0;
  y[133] = 0;
  y[134] = p[1]/sqrt(2) + p[10]/sqrt(2) - p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) - p[5]/sqrt(2) - p[6]/sqrt(2) + p[9]/sqrt(2);
  y[135] = 0;
  y[136] = 0;
  y[137] = 0;
  y[138] = 0;
  y[139] = p[1]/sqrt(2) + p[11]/sqrt(2) - p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) + p[4]/sqrt(2) - p[7]/sqrt(2) - p[8]/sqrt(2);
  y[140] = 0;
  y[141] = 0;
  y[142] = p[0]/sqrt(2) + p[10]/sqrt(2) - p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) + p[5]/sqrt(2) - p[6]/sqrt(2) - p[9]/sqrt(2);
  y[143] = 0;
  y[144] = 0;
  y[145] = p[0]/sqrt(2) - p[10]/sqrt(2) - p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) - p[5]/sqrt(2) + p[6]/sqrt(2) + p[9]/sqrt(2);
  y[146] = 0;
  y[147] = 0;
  y[148] = p[1]/sqrt(2) - p[11]/sqrt(2) - p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) - p[4]/sqrt(2) + p[7]/sqrt(2) + p[8]/sqrt(2);
  y[149] = 0;
  y[150] = 0;
  y[151] = -(sqrt(0.6666666666666666)*p[10]) + sqrt(0.6666666666666666)*p[11] + sqrt(0.6666666666666666)*p[4] - sqrt(0.6666666666666666)*p[5] - sqrt(0.6666666666666666)*p[6] + sqrt(0.6666666666666666)*p[7] + sqrt(0.6666666666666666)*p[8] - sqrt(0.6666666666666666)*p[9];
  y[152] = 0;
  y[153] = 0;
  y[154] = 0;
  y[155] = 0;
  y[156] = (2*sqrt(2)*p[1])/3. - (sqrt(2)*p[10])/3. - (sqrt(2)*p[11])/3. + (2*sqrt(2)*p[13])/3. + (2*sqrt(2)*p[14])/3. + (2*sqrt(2)*p[2])/3. - (sqrt(2)*p[4])/3. - (sqrt(2)*p[5])/3. - (sqrt(2)*p[6])/3. - (sqrt(2)*p[7])/3. - (sqrt(2)*p[8])/3. - (sqrt(2)*p[9])/3.;
  y[157] = 0;
  y[158] = 0;
  y[159] = 0;
  y[160] = 0;
  y[161] = p[0] - p[1]/3. - p[10]/3. - p[11]/3. + p[12] - p[13]/3. - p[14]/3. + p[15] - p[2]/3. + p[3] - p[4]/3. - p[5]/3. - p[6]/3. - p[7]/3. - p[8]/3. - p[9]/3.;
  y[162] = 0;
  y[163] = 0;
  y[164] = 0;
  y[165] = 0;
  y[166] = 0;
  y[167] = -p[10] - p[11] + p[4] + p[5] - p[6] - p[7] + p[8] + p[9];
  y[168] = 0;
  y[169] = 0;
  y[170] = 0;
  y[171] = 0;
  y[172] = 0;
  y[173] = 0;
  y[174] = 0;
  y[175] = 0;
  y[176] = 0;
  y[177] = p[0] + p[1] + p[12] + p[13] - p[14] - p[15] - p[2] - p[3];
  y[178] = 0;
  y[179] = 0;
  y[180] = 0;
  y[181] = 0;
  y[182] = 0;
  y[183] = 0;
  y[184] = p[0]/sqrt(2) + p[10]/sqrt(2) - p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) + p[5]/sqrt(2) - p[6]/sqrt(2) - p[9]/sqrt(2);
  y[185] = 0;
  y[186] = 0;
  y[187] = p[1]/sqrt(2) + p[11]/sqrt(2) - p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) + p[4]/sqrt(2) - p[7]/sqrt(2) - p[8]/sqrt(2);
  y[188] = 0;
  y[189] = 0;
  y[190] = p[1]/sqrt(2) - p[11]/sqrt(2) - p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) - p[4]/sqrt(2) + p[7]/sqrt(2) + p[8]/sqrt(2);
  y[191] = 0;
  y[192] = 0;
  y[193] = p[0]/sqrt(2) - p[10]/sqrt(2) - p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) - p[5]/sqrt(2) + p[6]/sqrt(2) + p[9]/sqrt(2);
  y[194] = 0;
  y[195] = 0;
  y[196] = 0;
  y[197] = 0;
  y[198] = p[1]/sqrt(2) - p[10]/sqrt(2) - p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) + p[5]/sqrt(2) + p[6]/sqrt(2) - p[9]/sqrt(2);
  y[199] = 0;
  y[200] = 0;
  y[201] = 0;
  y[202] = 0;
  y[203] = p[0]/sqrt(2) - p[11]/sqrt(2) - p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) + p[4]/sqrt(2) + p[7]/sqrt(2) - p[8]/sqrt(2);
  y[204] = p[1]/sqrt(2) + p[10]/sqrt(2) - p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) - p[5]/sqrt(2) - p[6]/sqrt(2) + p[9]/sqrt(2);
  y[205] = 0;
  y[206] = 0;
  y[207] = 0;
  y[208] = 0;
  y[209] = p[0]/sqrt(2) + p[11]/sqrt(2) - p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) - p[4]/sqrt(2) - p[7]/sqrt(2) + p[8]/sqrt(2);
  y[210] = 0;
  y[211] = 0;
  y[212] = 0;
  y[213] = 0;
  y[214] = 0;
  y[215] = p[10] - p[11] + p[4] - p[5] + p[6] - p[7] + p[8] - p[9];
  y[216] = 0;
  y[217] = 0;
  y[218] = 0;
  y[219] = 0;
  y[220] = 0;
  y[221] = 0;
  y[222] = 0;
  y[223] = 0;
  y[224] = 0;
  y[225] = p[0] - p[1] + p[12] - p[13] + p[14] - p[15] + p[2] - p[3];

  return 0;
}

// ----- Noise propagation

// probability distribution of depolarizing noise on 1 qubit
int dn_pdist_1Q(double p, double *pdist) {
  pdist[0] = (1.-p);
  pdist[1] = p/3.;
  pdist[2] = p/3.;
  pdist[3] = p/3.;
  return 0;
}

// probability distribution of depolarizing noise on n qubits
int dn_pdist_2Q(double p, double *pdist) {
  int i,j,n;

  // i == 0, j == 0
  pdist[0] = (1.-p)*(1.-p);

  // i == 0, j > 0
  for(j=1; j<4; j++) {
    n = get_linear_index(2, 4, {0,j});
    pdist[n] = (1.-p)*p/3.;
  }   
  for(i=1; i<4; i++) {
    // i > 0, j == 0
    n = get_linear_index(2, 4, {i,0});
    pdist[n] = (1.-p)*p/3.;

    // i,j > 0
    for(j=1; j<4; j++) {
      n = get_linear_index(2, 4, {i,j});
      pdist[n] = p*p/9.;
    }     
  }
}

// This is the convolution in (Z_2)^2n
int convolve_mod2(double *p1, double *p2, double *p3, int nqubits) {
  int len = pow(4,nqubits);
  int *is = new int[nqubits];
  int *js = new int[nqubits];
  int *sum = new int[nqubits];
  int i,j,n;

  for(i=0; i<len; i++) {
    get_multi_index(nqubits, 4, i, is);
    p3[i] = 0;

    for(j=0; j<len; j++) {
      // n
      get_multi_index(nqubits, 4, j, js);
      for(n=0; n<nqubits; n++) {
        sum[n] = is[n]^js[n];
      }
      n = get_linear_index(nqubits, 4, sum);
      p3[i] += p1[j] * p2[n];
    } 
  }
  return 0;
}

// transform probibility distribution
int apply_symplectic_map(double *pin, double *pout, int nqubits, symplectic_map S) {
  return 0;
}


// int Hadamard1(double *pin, double *pout, int nqubits);

// int Hadamard2(double *pin, double *pout, int nqubits);

// int S1(double *pin, double *pout, int nqubits);

// int S2(double *pin, double *pout, int nqubits);

// int CNOT12(double *pin, double *pout, int nqubits);

#endif