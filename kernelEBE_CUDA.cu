#include <cuda_runtime.h>
#include <stdio.h>

#define m00238 (0.002380952380952381)
#define mm0238 (-0.002380952380952381)
#define m00039 (0.0003968253968253968)
#define mm0158 (-0.0015873015873015873)
#define m01269 (0.012698412698412698)
#define m00634 (0.006349206349206349)
#define m00317 (0.0031746031746031746)
#define NIP 6

__global__
void clear_zero_cuda_kernel(
 int n,  double* rg)
{
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if( i >= n ) return;
  rg[i]=0.0;
}

__global__
void kernelEBE_CUDA_kernel(
 int ne,
 double* ug, double* rg,
 int* cny,
 double* d1, double* d2, double* d4, double* d5, double* d6,
 double* c1, double* c2, double* c3, double* c4, double* c5,
 double* c6, double* c7, double* c8, double* c9,
 double* tmp1, double* tmp2, double* coe)
{
  int ie = (blockDim.x * blockIdx.x + threadIdx.x)/NVEC;
  int ivec = (blockDim.x * blockIdx.x + threadIdx.x)%NVEC;
  
  if( ie >= ne ) return;
  if( ivec >= NVEC ) return;
  
  double u1=ug[ivec+(3*(cny[ie*10+0]-1)+0)*NVEC];
  double u2=ug[ivec+(3*(cny[ie*10+0]-1)+1)*NVEC];
  double u3=ug[ivec+(3*(cny[ie*10+0]-1)+2)*NVEC];
  double u4=ug[ivec+(3*(cny[ie*10+1]-1)+0)*NVEC];
  double u5=ug[ivec+(3*(cny[ie*10+1]-1)+1)*NVEC];
  double u6=ug[ivec+(3*(cny[ie*10+1]-1)+2)*NVEC];
  double u7=ug[ivec+(3*(cny[ie*10+2]-1)+0)*NVEC];
  double u8=ug[ivec+(3*(cny[ie*10+2]-1)+1)*NVEC];
  double u9=ug[ivec+(3*(cny[ie*10+2]-1)+2)*NVEC];
  double u10=ug[ivec+(3*(cny[ie*10+3]-1)+0)*NVEC];
  double u11=ug[ivec+(3*(cny[ie*10+3]-1)+1)*NVEC];
  double u12=ug[ivec+(3*(cny[ie*10+3]-1)+2)*NVEC];
  double u13=ug[ivec+(3*(cny[ie*10+4]-1)+0)*NVEC];
  double u14=ug[ivec+(3*(cny[ie*10+4]-1)+1)*NVEC];
  double u15=ug[ivec+(3*(cny[ie*10+4]-1)+2)*NVEC];
  double u16=ug[ivec+(3*(cny[ie*10+5]-1)+0)*NVEC];
  double u17=ug[ivec+(3*(cny[ie*10+5]-1)+1)*NVEC];
  double u18=ug[ivec+(3*(cny[ie*10+5]-1)+2)*NVEC];
  double u19=ug[ivec+(3*(cny[ie*10+6]-1)+0)*NVEC];
  double u20=ug[ivec+(3*(cny[ie*10+6]-1)+1)*NVEC];
  double u21=ug[ivec+(3*(cny[ie*10+6]-1)+2)*NVEC];
  double u22=ug[ivec+(3*(cny[ie*10+7]-1)+0)*NVEC];
  double u23=ug[ivec+(3*(cny[ie*10+7]-1)+1)*NVEC];
  double u24=ug[ivec+(3*(cny[ie*10+7]-1)+2)*NVEC];
  double u25=ug[ivec+(3*(cny[ie*10+8]-1)+0)*NVEC];
  double u26=ug[ivec+(3*(cny[ie*10+8]-1)+1)*NVEC];
  double u27=ug[ivec+(3*(cny[ie*10+8]-1)+2)*NVEC];
  double u28=ug[ivec+(3*(cny[ie*10+9]-1)+0)*NVEC];
  double u29=ug[ivec+(3*(cny[ie*10+9]-1)+1)*NVEC];
  double u30=ug[ivec+(3*(cny[ie*10+9]-1)+2)*NVEC];

  double detj = coe[ie];
  double dadx11 = c1[ie];
  double dadx21 = c2[ie];
  double dadx31 = c3[ie];
  double dadx12 = c4[ie];
  double dadx22 = c5[ie];
  double dadx32 = c6[ie];
  double dadx13 = c7[ie];
  double dadx23 = c8[ie];
  double dadx33 = c9[ie];

  double BDBu1  = 0.0;
  double BDBu2  = 0.0;
  double BDBu3  = 0.0;
  double BDBu4  = 0.0;
  double BDBu5  = 0.0;
  double BDBu6  = 0.0;
  double BDBu7  = 0.0;
  double BDBu8  = 0.0;
  double BDBu9  = 0.0;
  double BDBu10 = 0.0;
  double BDBu11 = 0.0;
  double BDBu12 = 0.0;
  double BDBu13 = 0.0;
  double BDBu14 = 0.0;
  double BDBu15 = 0.0;
  double BDBu16 = 0.0;
  double BDBu17 = 0.0;
  double BDBu18 = 0.0;
  double BDBu19 = 0.0;
  double BDBu20 = 0.0;
  double BDBu21 = 0.0;
  double BDBu22 = 0.0;
  double BDBu23 = 0.0;
  double BDBu24 = 0.0;
  double BDBu25 = 0.0;
  double BDBu26 = 0.0;
  double BDBu27 = 0.0;
  double BDBu28 = 0.0;
  double BDBu29 = 0.0;
  double BDBu30 = 0.0;
  
  double D_11 = 2.0*(d1[ie]-d2[ie])/3.0;
  double D_22 = D_11;
  double D_33 = D_11;
  double D_12 = (-d1[ie]+d2[ie])/3.0;
  double D_13 = D_12;
  double D_23 = D_12;
  double D_44 = d4[ie];
  double D_55 = d5[ie];
  double D_66 = d6[ie];
  double alpha_iip[3*NIP],weight[NIP];
  weight[0]=-4.0/30.0;
  weight[1]=3.0/40.0;
  weight[2]=3.0/40.0;
  weight[3]=3.0/40.0;
  weight[4]=3.0/40.0;
  weight[5]=1.0/6.0;
  alpha_iip[0]=0.25;
  alpha_iip[1]=0.25;
  alpha_iip[2]=0.25;
  alpha_iip[3]=1.0/6.0;
  alpha_iip[4]=1.0/6.0;
  alpha_iip[5]=1.0/6.0;
  alpha_iip[6]=0.5;
  alpha_iip[7]=1.0/6.0;
  alpha_iip[8]=1.0/6.0;
  alpha_iip[9]=1.0/6.0;
  alpha_iip[10]=0.5;
  alpha_iip[11]=1.0/6.0;
  alpha_iip[12]=1.0/6.0;
  alpha_iip[13]=1.0/6.0;
  alpha_iip[14]=0.5;
  alpha_iip[15]=0.25;
  alpha_iip[16]=0.25;
  alpha_iip[17]=0.25;

  for( int iip = 0; iip < NIP; ++iip )
  {
    if(iip==(NIP-1))
    {
      D_11 = (d1[ie]+2.0*d2[ie])/3.0;
      D_22 = D_11;
      D_33 = D_11;
      D_12 = D_11;
      D_13 = D_11;
      D_23 = D_11;
      D_44 = 0.0;
      D_55 = 0.0;
      D_66 = 0.0;
    }
    double dNda31 = 4.0*alpha_iip[3*iip+2] - 1.0;
    double dNda22 = 4.0*alpha_iip[3*iip+1] - 1.0;
    double dNda13 = 4.0*alpha_iip[3*iip+0] - 1.0;
    double dNda14 = 4.0*(alpha_iip[3*iip+0]+alpha_iip[3*iip+1]
                      +alpha_iip[3*iip+2]) - 3.0;
    double dNda24 = 4.0*(alpha_iip[3*iip+0]+alpha_iip[3*iip+1]
                      +alpha_iip[3*iip+2]) - 3.0;
    double dNda34 = 4.0*(alpha_iip[3*iip+0]+alpha_iip[3*iip+1]
                      +alpha_iip[3*iip+2]) - 3.0;
    double dNda25 = 4.0*alpha_iip[3*iip+2];
    double dNda35 = 4.0*alpha_iip[3*iip+1];
    double dNda16 = 4.0*alpha_iip[3*iip+1];
    double dNda26 = 4.0*alpha_iip[3*iip+0];
    double dNda17 = 4.0*alpha_iip[3*iip+2];
    double dNda37 = 4.0*alpha_iip[3*iip+0];
    double dNda18 = -4.0*alpha_iip[3*iip+2];
    double dNda28 = -4.0*alpha_iip[3*iip+2];
    double dNda38 = -4.0*alpha_iip[3*iip+0] -4.0*alpha_iip[3*iip+1]
                -8.0*alpha_iip[3*iip+2] + 4.0;
    double dNda19 = -4.0*alpha_iip[3*iip+1];
    double dNda29 = -4.0*alpha_iip[3*iip+0] -8.0*alpha_iip[3*iip+1]
                -4.0*alpha_iip[3*iip+2] + 4.0;
    double dNda39 = -4.0*alpha_iip[3*iip+1];
    double dNda1a = -8.0*alpha_iip[3*iip+0] -4.0*alpha_iip[3*iip+1]
                -4.0*alpha_iip[3*iip+2] + 4.0;
    double dNda2a = -4.0*alpha_iip[3*iip+0];
    double dNda3a = -4.0*alpha_iip[3*iip+0];

    double dNdx11 = dNda31*dadx31;
    double dNdx12 = dNda22*dadx21;
    double dNdx13 = dNda13*dadx11;
    double dNdx14 = dNda14*dadx11+dNda24*dadx21+dNda34*dadx31;
    double dNdx15 = dNda25*dadx21+dNda35*dadx31;
    double dNdx16 = dNda16*dadx11+dNda26*dadx21;
    double dNdx17 = dNda17*dadx11+dNda37*dadx31;
    double dNdx18 = dNda18*dadx11+dNda28*dadx21+dNda38*dadx31;
    double dNdx19 = dNda19*dadx11+dNda29*dadx21+dNda39*dadx31;
    double dNdx1a = dNda1a*dadx11+dNda2a*dadx21+dNda3a*dadx31;
    double dNdx21 = dNda31*dadx32;
    double dNdx22 = dNda22*dadx22;
    double dNdx23 = dNda13*dadx12;
    double dNdx24 = dNda14*dadx12+dNda24*dadx22+dNda34*dadx32;
    double dNdx25 = dNda25*dadx22+dNda35*dadx32;
    double dNdx26 = dNda16*dadx12+dNda26*dadx22;
    double dNdx27 = dNda17*dadx12+dNda37*dadx32;
    double dNdx28 = dNda18*dadx12+dNda28*dadx22+dNda38*dadx32;
    double dNdx29 = dNda19*dadx12+dNda29*dadx22+dNda39*dadx32;
    double dNdx2a = dNda1a*dadx12+dNda2a*dadx22+dNda3a*dadx32;
    double dNdx31 = dNda31*dadx33;
    double dNdx32 = dNda22*dadx23;
    double dNdx33 = dNda13*dadx13;
    double dNdx34 = dNda14*dadx13+dNda24*dadx23+dNda34*dadx33;
    double dNdx35 = dNda25*dadx23+dNda35*dadx33;
    double dNdx36 = dNda16*dadx13+dNda26*dadx23;
    double dNdx37 = dNda17*dadx13+dNda37*dadx33;
    double dNdx38 = dNda18*dadx13+dNda28*dadx23+dNda38*dadx33;
    double dNdx39 = dNda19*dadx13+dNda29*dadx23+dNda39*dadx33;
    double dNdx3a = dNda1a*dadx13+dNda2a*dadx23+dNda3a*dadx33;

    double Bu1 = u1*dNdx11;
    double Bu2 = u2*dNdx21;
    double Bu3 = u3*dNdx31;
    double Bu4 = u1*dNdx21 + u2*dNdx11;
    double Bu5 = u2*dNdx31 + u3*dNdx21;
    double Bu6 = u3*dNdx11 + u1*dNdx31;
    Bu1 = Bu1 + u4*dNdx12;
    Bu2 = Bu2 + u5*dNdx22;
    Bu3 = Bu3 + u6*dNdx32;
    Bu4 = Bu4 + u4*dNdx22 + u5*dNdx12;
    Bu5 = Bu5 + u5*dNdx32 + u6*dNdx22;
    Bu6 = Bu6 + u6*dNdx12 + u4*dNdx32;
    Bu1 = Bu1 + u7*dNdx13;
    Bu2 = Bu2 + u8*dNdx23;
    Bu3 = Bu3 + u9*dNdx33;
    Bu4 = Bu4 + u7*dNdx23 + u8*dNdx13;
    Bu5 = Bu5 + u8*dNdx33 + u9*dNdx23;
    Bu6 = Bu6 + u9*dNdx13 + u7*dNdx33;
    Bu1 = Bu1 + u10*dNdx14;
    Bu2 = Bu2 + u11*dNdx24;
    Bu3 = Bu3 + u12*dNdx34;
    Bu4 = Bu4 + u10*dNdx24 + u11*dNdx14;
    Bu5 = Bu5 + u11*dNdx34 + u12*dNdx24;
    Bu6 = Bu6 + u12*dNdx14 + u10*dNdx34;
    Bu1 = Bu1 + u13*dNdx15;
    Bu2 = Bu2 + u14*dNdx25;
    Bu3 = Bu3 + u15*dNdx35;
    Bu4 = Bu4 + u13*dNdx25 + u14*dNdx15;
    Bu5 = Bu5 + u14*dNdx35 + u15*dNdx25;
    Bu6 = Bu6 + u15*dNdx15 + u13*dNdx35;
    Bu1 = Bu1 + u16*dNdx16;
    Bu2 = Bu2 + u17*dNdx26;
    Bu3 = Bu3 + u18*dNdx36;
    Bu4 = Bu4 + u16*dNdx26 + u17*dNdx16;
    Bu5 = Bu5 + u17*dNdx36 + u18*dNdx26;
    Bu6 = Bu6 + u18*dNdx16 + u16*dNdx36;
    Bu1 = Bu1 + u19*dNdx17;
    Bu2 = Bu2 + u20*dNdx27;
    Bu3 = Bu3 + u21*dNdx37;
    Bu4 = Bu4 + u19*dNdx27 + u20*dNdx17;
    Bu5 = Bu5 + u20*dNdx37 + u21*dNdx27;
    Bu6 = Bu6 + u21*dNdx17 + u19*dNdx37;
    Bu1 = Bu1 + u22*dNdx18;
    Bu2 = Bu2 + u23*dNdx28;
    Bu3 = Bu3 + u24*dNdx38;
    Bu4 = Bu4 + u22*dNdx28 + u23*dNdx18;
    Bu5 = Bu5 + u23*dNdx38 + u24*dNdx28;
    Bu6 = Bu6 + u24*dNdx18 + u22*dNdx38;
    Bu1 = Bu1 + u25*dNdx19;
    Bu2 = Bu2 + u26*dNdx29;
    Bu3 = Bu3 + u27*dNdx39;
    Bu4 = Bu4 + u25*dNdx29 + u26*dNdx19;
    Bu5 = Bu5 + u26*dNdx39 + u27*dNdx29;
    Bu6 = Bu6 + u27*dNdx19 + u25*dNdx39;
    Bu1 = Bu1 + u28*dNdx1a;
    Bu2 = Bu2 + u29*dNdx2a;
    Bu3 = Bu3 + u30*dNdx3a;
    Bu4 = Bu4 + u28*dNdx2a + u29*dNdx1a;
    Bu5 = Bu5 + u29*dNdx3a + u30*dNdx2a;
    Bu6 = Bu6 + u30*dNdx1a + u28*dNdx3a;
    double DBu1 = D_11*Bu1+D_12*Bu2+D_13*Bu3;
    double DBu2 = D_12*Bu1+D_22*Bu2+D_23*Bu3;
    double DBu3 = D_13*Bu1+D_23*Bu2+D_33*Bu3;
    double DBu4 = D_44*Bu4;
    double DBu5 = D_55*Bu5;
    double DBu6 = D_66*Bu6;

    double BDBu_iip1 = dNdx11*DBu1 + dNdx21*DBu4 + dNdx31*DBu6;
    double BDBu_iip2 = dNdx11*DBu4 + dNdx21*DBu2 + dNdx31*DBu5;
    double BDBu_iip3 = dNdx11*DBu6 + dNdx21*DBu5 + dNdx31*DBu3;
    double BDBu_iip4 = dNdx12*DBu1 + dNdx22*DBu4 + dNdx32*DBu6;
    double BDBu_iip5 = dNdx12*DBu4 + dNdx22*DBu2 + dNdx32*DBu5;
    double BDBu_iip6 = dNdx12*DBu6 + dNdx22*DBu5 + dNdx32*DBu3;
    double BDBu_iip7 = dNdx13*DBu1 + dNdx23*DBu4 + dNdx33*DBu6;
    double BDBu_iip8 = dNdx13*DBu4 + dNdx23*DBu2 + dNdx33*DBu5;
    double BDBu_iip9 = dNdx13*DBu6 + dNdx23*DBu5 + dNdx33*DBu3;
    double BDBu_iip10 = dNdx14*DBu1 + dNdx24*DBu4 + dNdx34*DBu6;
    double BDBu_iip11 = dNdx14*DBu4 + dNdx24*DBu2 + dNdx34*DBu5;
    double BDBu_iip12 = dNdx14*DBu6 + dNdx24*DBu5 + dNdx34*DBu3;
    double BDBu_iip13 = dNdx15*DBu1 + dNdx25*DBu4 + dNdx35*DBu6;
    double BDBu_iip14 = dNdx15*DBu4 + dNdx25*DBu2 + dNdx35*DBu5;
    double BDBu_iip15 = dNdx15*DBu6 + dNdx25*DBu5 + dNdx35*DBu3;
    double BDBu_iip16 = dNdx16*DBu1 + dNdx26*DBu4 + dNdx36*DBu6;
    double BDBu_iip17 = dNdx16*DBu4 + dNdx26*DBu2 + dNdx36*DBu5;
    double BDBu_iip18 = dNdx16*DBu6 + dNdx26*DBu5 + dNdx36*DBu3;
    double BDBu_iip19 = dNdx17*DBu1 + dNdx27*DBu4 + dNdx37*DBu6;
    double BDBu_iip20 = dNdx17*DBu4 + dNdx27*DBu2 + dNdx37*DBu5;
    double BDBu_iip21 = dNdx17*DBu6 + dNdx27*DBu5 + dNdx37*DBu3;
    double BDBu_iip22 = dNdx18*DBu1 + dNdx28*DBu4 + dNdx38*DBu6;
    double BDBu_iip23 = dNdx18*DBu4 + dNdx28*DBu2 + dNdx38*DBu5;
    double BDBu_iip24 = dNdx18*DBu6 + dNdx28*DBu5 + dNdx38*DBu3;
    double BDBu_iip25 = dNdx19*DBu1 + dNdx29*DBu4 + dNdx39*DBu6;
    double BDBu_iip26 = dNdx19*DBu4 + dNdx29*DBu2 + dNdx39*DBu5;
    double BDBu_iip27 = dNdx19*DBu6 + dNdx29*DBu5 + dNdx39*DBu3;
    double BDBu_iip28 = dNdx1a*DBu1 + dNdx2a*DBu4 + dNdx3a*DBu6;
    double BDBu_iip29 = dNdx1a*DBu4 + dNdx2a*DBu2 + dNdx3a*DBu5;
    double BDBu_iip30 = dNdx1a*DBu6 + dNdx2a*DBu5 + dNdx3a*DBu3;
    BDBu1 = BDBu1 + weight[iip]*BDBu_iip1;
    BDBu2 = BDBu2 + weight[iip]*BDBu_iip2;
    BDBu3 = BDBu3 + weight[iip]*BDBu_iip3;
    BDBu4 = BDBu4 + weight[iip]*BDBu_iip4;
    BDBu5 = BDBu5 + weight[iip]*BDBu_iip5;
    BDBu6 = BDBu6 + weight[iip]*BDBu_iip6;
    BDBu7 = BDBu7 + weight[iip]*BDBu_iip7;
    BDBu8 = BDBu8 + weight[iip]*BDBu_iip8;
    BDBu9 = BDBu9 + weight[iip]*BDBu_iip9;
    BDBu10 = BDBu10 + weight[iip]*BDBu_iip10;
    BDBu11 = BDBu11 + weight[iip]*BDBu_iip11;
    BDBu12 = BDBu12 + weight[iip]*BDBu_iip12;
    BDBu13 = BDBu13 + weight[iip]*BDBu_iip13;
    BDBu14 = BDBu14 + weight[iip]*BDBu_iip14;
    BDBu15 = BDBu15 + weight[iip]*BDBu_iip15;
    BDBu16 = BDBu16 + weight[iip]*BDBu_iip16;
    BDBu17 = BDBu17 + weight[iip]*BDBu_iip17;
    BDBu18 = BDBu18 + weight[iip]*BDBu_iip18;
    BDBu19 = BDBu19 + weight[iip]*BDBu_iip19;
    BDBu20 = BDBu20 + weight[iip]*BDBu_iip20;
    BDBu21 = BDBu21 + weight[iip]*BDBu_iip21;
    BDBu22 = BDBu22 + weight[iip]*BDBu_iip22;
    BDBu23 = BDBu23 + weight[iip]*BDBu_iip23;
    BDBu24 = BDBu24 + weight[iip]*BDBu_iip24;
    BDBu25 = BDBu25 + weight[iip]*BDBu_iip25;
    BDBu26 = BDBu26 + weight[iip]*BDBu_iip26;
    BDBu27 = BDBu27 + weight[iip]*BDBu_iip27;
    BDBu28 = BDBu28 + weight[iip]*BDBu_iip28;
    BDBu29 = BDBu29 + weight[iip]*BDBu_iip29;
    BDBu30 = BDBu30 + weight[iip]*BDBu_iip30;
  }
  BDBu1 = BDBu1*detj;
  BDBu2 = BDBu2*detj;
  BDBu3 = BDBu3*detj;
  BDBu4 = BDBu4*detj;
  BDBu5 = BDBu5*detj;
  BDBu6 = BDBu6*detj;
  BDBu7 = BDBu7*detj;
  BDBu8 = BDBu8*detj;
  BDBu9 = BDBu9*detj;
  BDBu10 = BDBu10*detj;
  BDBu11 = BDBu11*detj;
  BDBu12 = BDBu12*detj;
  BDBu13 = BDBu13*detj;
  BDBu14 = BDBu14*detj;
  BDBu15 = BDBu15*detj;
  BDBu16 = BDBu16*detj;
  BDBu17 = BDBu17*detj;
  BDBu18 = BDBu18*detj;
  BDBu19 = BDBu19*detj;
  BDBu20 = BDBu20*detj;
  BDBu21 = BDBu21*detj;
  BDBu22 = BDBu22*detj;
  BDBu23 = BDBu23*detj;
  BDBu24 = BDBu24*detj;
  BDBu25 = BDBu25*detj;
  BDBu26 = BDBu26*detj;
  BDBu27 = BDBu27*detj;
  BDBu28 = BDBu28*detj;
  BDBu29 = BDBu29*detj;
  BDBu30 = BDBu30*detj;

  BDBu1=(
  m00238*u1+
  m00039*(u4+u7+u10)+
  mm0158*(u13+u19+u22)+
  mm0238*(u16+u25+u28)
  )*tmp1[ie]+BDBu1*tmp2[ie];
  BDBu4=(
  m00039*(u1+u7+u10)+
  m00238*u4+
  mm0158*(u13+u16+u25)+
  mm0238*(u19+u22+u28)
  )*tmp1[ie]+BDBu4*tmp2[ie];
  BDBu7=(
  m00039*(u1+u4+u10)+
  m00238*u7+
  mm0238*(u13+u22+u25)+
  mm0158*(u16+u19+u28)
  )*tmp1[ie]+BDBu7*tmp2[ie];
  BDBu10=(
  m00039*(u1+u4+u7)+
  m00238*u10+
  mm0238*(u13+u16+u19)+
  mm0158*(u22+u25+u28)
  )*tmp1[ie]+BDBu10*tmp2[ie];
  BDBu13=(
  mm0158*(u1+u4)+
  mm0238*(u7+u10)+
  m01269*u13+
  m00634*(u16+u19+u22+u25)+
  m00317*u28
  )*tmp1[ie]+BDBu13*tmp2[ie];
  BDBu16=(
  mm0238*(u1+u10)+
  mm0158*(u4+u7)+
  m00634*(u13+u19+u25+u28)+
  m01269*u16+
  m00317*u22
  )*tmp1[ie]+BDBu16*tmp2[ie];
  BDBu19=(
  mm0158*(u1+u7)+
  mm0238*(u4+u10)+
  m00634*(u13+u16+u22+u28)+
  m01269*u19+
  m00317*u25
  )*tmp1[ie]+BDBu19*tmp2[ie];
  BDBu22=(
  mm0158*(u1+u10)+
  mm0238*(u4+u7)+
  m00634*(u13+u19+u25+u28)+
  m00317*u16+
  m01269*u22
  )*tmp1[ie]+BDBu22*tmp2[ie];
  BDBu25=(
  mm0238*(u1+u7)+
  mm0158*(u4+u10)+
  m00634*(u13+u16+u22+u28)+
  m00317*u19+
  m01269*u25
  )*tmp1[ie]+BDBu25*tmp2[ie];
  BDBu28=(
  mm0238*(u1+u4)+
  mm0158*(u7+u10)+
  m00317*u13+
  m00634*(u16+u19+u22+u25)+
  m01269*u28
  )*tmp1[ie]+BDBu28*tmp2[ie];
  BDBu2=(
  m00238*u2+
  m00039*(u5+u8+u11)+
  mm0158*(u14+u20+u23)+
  mm0238*(u17+u26+u29)
  )*tmp1[ie]+BDBu2*tmp2[ie];
  BDBu5=(
  m00039*(u2+u8+u11)+
  m00238*u5+
  mm0158*(u14+u17+u26)+
  mm0238*(u20+u23+u29)
  )*tmp1[ie]+BDBu5*tmp2[ie];
  BDBu8=(
  m00039*(u2+u5+u11)+
  m00238*u8+
  mm0238*(u14+u23+u26)+
  mm0158*(u17+u20+u29)
  )*tmp1[ie]+BDBu8*tmp2[ie];
  BDBu11=(
  m00039*(u2+u5+u8)+
  m00238*u11+
  mm0238*(u14+u17+u20)+
  mm0158*(u23+u26+u29)
  )*tmp1[ie]+BDBu11*tmp2[ie];
  BDBu14=(
  mm0158*(u2+u5)+
  mm0238*(u8+u11)+
  m01269*u14+
  m00634*(u17+u20+u23+u26)+
  m00317*u29
  )*tmp1[ie]+BDBu14*tmp2[ie];
  BDBu17=(
  mm0238*(u2+u11)+
  mm0158*(u5+u8)+
  m00634*(u14+u20+u26+u29)+
  m01269*u17+
  m00317*u23
  )*tmp1[ie]+BDBu17*tmp2[ie];
  BDBu20=(
  mm0158*(u2+u8)+
  mm0238*(u5+u11)+
  m00634*(u14+u17+u23+u29)+
  m01269*u20+
  m00317*u26
  )*tmp1[ie]+BDBu20*tmp2[ie];
  BDBu23=(
  mm0158*(u2+u11)+
  mm0238*(u5+u8)+
  m00634*(u14+u20+u26+u29)+
  m00317*u17+
  m01269*u23
  )*tmp1[ie]+BDBu23*tmp2[ie];
  BDBu26=(
  mm0238*(u2+u8)+
  mm0158*(u5+u11)+
  m00634*(u14+u17+u23+u29)+
  m00317*u20+
  m01269*u26
  )*tmp1[ie]+BDBu26*tmp2[ie];
  BDBu29=(
  mm0238*(u2+u5)+
  mm0158*(u8+u11)+
  m00317*u14+
  m00634*(u17+u20+u23+u26)+
  m01269*u29
  )*tmp1[ie]+BDBu29*tmp2[ie];
  BDBu3=(
  m00238*u3+
  m00039*(u6+u9+u12)+
  mm0158*(u15+u21+u24)+
  mm0238*(u18+u27+u30)
  )*tmp1[ie]+BDBu3*tmp2[ie];
  BDBu6=(
  m00039*(u3+u9+u12)+
  m00238*u6+
  mm0158*(u15+u18+u27)+
  mm0238*(u21+u24+u30)
  )*tmp1[ie]+BDBu6*tmp2[ie];
  BDBu9=(
  m00039*(u3+u6+u12)+
  m00238*u9+
  mm0238*(u15+u24+u27)+
  mm0158*(u18+u21+u30)
  )*tmp1[ie]+BDBu9*tmp2[ie];
  BDBu12=(
  m00039*(u3+u6+u9)+
  m00238*u12+
  mm0238*(u15+u18+u21)+
  mm0158*(u24+u27+u30)
  )*tmp1[ie]+BDBu12*tmp2[ie];
  BDBu15=(
  mm0158*(u3+u6)+
  mm0238*(u9+u12)+
  m01269*u15+
  m00634*(u18+u21+u24+u27)+
  m00317*u30
  )*tmp1[ie]+BDBu15*tmp2[ie];
  BDBu18=(
  mm0238*(u3+u12)+
  mm0158*(u6+u9)+
  m00634*(u15+u21+u27+u30)+
  m01269*u18+
  m00317*u24
  )*tmp1[ie]+BDBu18*tmp2[ie];
  BDBu21=(
  mm0158*(u3+u9)+
  mm0238*(u6+u12)+
  m00634*(u15+u18+u24+u30)+
  m01269*u21+
  m00317*u27
  )*tmp1[ie]+BDBu21*tmp2[ie];
  BDBu24=(
  mm0158*(u3+u12)+
  mm0238*(u6+u9)+
  m00634*(u15+u21+u27+u30)+
  m00317*u18+
  m01269*u24
  )*tmp1[ie]+BDBu24*tmp2[ie];
  BDBu27=(
  mm0238*(u3+u9)+
  mm0158*(u6+u12)+
  m00634*(u15+u18+u24+u30)+
  m00317*u21+
  m01269*u27
  )*tmp1[ie]+BDBu27*tmp2[ie];
  BDBu30=(
  mm0238*(u3+u6)+
  mm0158*(u9+u12)+
  m00317*u15+
  m00634*(u18+u21+u24+u27)+
  m01269*u30
  )*tmp1[ie]+BDBu30*tmp2[ie];
  atomicAdd(rg+ivec+(3*(cny[10*ie+0]-1)+0)*NVEC, BDBu1);
  atomicAdd(rg+ivec+(3*(cny[10*ie+0]-1)+1)*NVEC, BDBu2);
  atomicAdd(rg+ivec+(3*(cny[10*ie+0]-1)+2)*NVEC, BDBu3);
  atomicAdd(rg+ivec+(3*(cny[10*ie+1]-1)+0)*NVEC, BDBu4);
  atomicAdd(rg+ivec+(3*(cny[10*ie+1]-1)+1)*NVEC, BDBu5);
  atomicAdd(rg+ivec+(3*(cny[10*ie+1]-1)+2)*NVEC, BDBu6);
  atomicAdd(rg+ivec+(3*(cny[10*ie+2]-1)+0)*NVEC, BDBu7);
  atomicAdd(rg+ivec+(3*(cny[10*ie+2]-1)+1)*NVEC, BDBu8);
  atomicAdd(rg+ivec+(3*(cny[10*ie+2]-1)+2)*NVEC, BDBu9);
  atomicAdd(rg+ivec+(3*(cny[10*ie+3]-1)+0)*NVEC, BDBu10);
  atomicAdd(rg+ivec+(3*(cny[10*ie+3]-1)+1)*NVEC, BDBu11);
  atomicAdd(rg+ivec+(3*(cny[10*ie+3]-1)+2)*NVEC, BDBu12);
  atomicAdd(rg+ivec+(3*(cny[10*ie+4]-1)+0)*NVEC, BDBu13);
  atomicAdd(rg+ivec+(3*(cny[10*ie+4]-1)+1)*NVEC, BDBu14);
  atomicAdd(rg+ivec+(3*(cny[10*ie+4]-1)+2)*NVEC, BDBu15);
  atomicAdd(rg+ivec+(3*(cny[10*ie+5]-1)+0)*NVEC, BDBu16);
  atomicAdd(rg+ivec+(3*(cny[10*ie+5]-1)+1)*NVEC, BDBu17);
  atomicAdd(rg+ivec+(3*(cny[10*ie+5]-1)+2)*NVEC, BDBu18);
  atomicAdd(rg+ivec+(3*(cny[10*ie+6]-1)+0)*NVEC, BDBu19);
  atomicAdd(rg+ivec+(3*(cny[10*ie+6]-1)+1)*NVEC, BDBu20);
  atomicAdd(rg+ivec+(3*(cny[10*ie+6]-1)+2)*NVEC, BDBu21);
  atomicAdd(rg+ivec+(3*(cny[10*ie+7]-1)+0)*NVEC, BDBu22);
  atomicAdd(rg+ivec+(3*(cny[10*ie+7]-1)+1)*NVEC, BDBu23);
  atomicAdd(rg+ivec+(3*(cny[10*ie+7]-1)+2)*NVEC, BDBu24);
  atomicAdd(rg+ivec+(3*(cny[10*ie+8]-1)+0)*NVEC, BDBu25);
  atomicAdd(rg+ivec+(3*(cny[10*ie+8]-1)+1)*NVEC, BDBu26);
  atomicAdd(rg+ivec+(3*(cny[10*ie+8]-1)+2)*NVEC, BDBu27);
  atomicAdd(rg+ivec+(3*(cny[10*ie+9]-1)+0)*NVEC, BDBu28);
  atomicAdd(rg+ivec+(3*(cny[10*ie+9]-1)+1)*NVEC, BDBu29);
  atomicAdd(rg+ivec+(3*(cny[10*ie+9]-1)+2)*NVEC, BDBu30);
}

extern "C" void kernelebe_cuda_(
 int* n_size, int* ne_size, int* n,int* ne,
 double* ug, double* rg,
 int* cny,
 double* d1, double* d2, double* d4, double* d5, double* d6,
 double* c1, double* c2, double* c3, double* c4, double* c5,
 double* c6, double* c7, double* c8, double* c9,
 double* tmp1, double* tmp2, double* coe)
{
  dim3 griddim, blockdim;
  int netmp(ne[0]);

  blockdim = dim3(128,1,1);
  griddim = dim3((3*n[0]*NVEC)/blockdim.x + 1,1,1);
  clear_zero_cuda_kernel<<<griddim,blockdim>>>(n[0]*3*NVEC,rg);
  cudaDeviceSynchronize();

  blockdim = dim3(32,1,1);
  griddim = dim3((netmp*NVEC)/blockdim.x + 1,1,1);
  kernelEBE_CUDA_kernel<<<griddim,blockdim>>>(netmp,ug,rg,cny,d1,d2,d4,d5,d6,c1,c2,c3,c4,c5,c6,c7,c8,c9,tmp1,tmp2,coe);
  cudaDeviceSynchronize();
}
