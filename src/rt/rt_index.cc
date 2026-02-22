#include <cmath>
#include <string.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <mpi.h>
#include "mem.h"
#include "grid.H"
#include "run.H"
#include "physics.H"
#include "rt.h"

using std::min;
using std::max;
using std::cout;
using std::endl;

RTS *rt_new(GridData &Grid,RunData &Run,PhysicsData &Physics)
{
  int rttype;
  if(Run.rank==0){ // check solver type only
    rttype = Physics.rt[i_rt_type];
  }

  MPI_Bcast(&rttype,1,MPI_INT,0,MPI_COMM_WORLD);

  switch(rttype){
      case(RT_DEFAULT): return new RTS(Grid,Run,Physics);
      default: return 0;
  }
}

RTS::RTS(GridData &Grid, RunData &Run, PhysicsData &Physics) {
  NDIM = Grid.NDIM;
  rttype = Physics.rt[i_rt_type];
  eps_const = Physics.rt[i_rt_epsilon];

  verbose=Run.verbose;
  call_count=0;

  xl=Grid.lbeg[1];
  xh=Grid.lend[1];
  yl=Grid.lbeg[2];
  yh=Grid.lend[2];
  zl=Grid.lbeg[0];
  zh=Grid.lend[0];

  // Ghost adjustment
  xl-=(xo=(Grid.lbeg[1]-Grid.lend[1])?1:0);
  yl-=(yo=(Grid.lbeg[2]-Grid.lend[2])?1:0);
  zl-=(zo=1);

  nx=xh-xl+1;
  ny=yh-yl+1;
  nz=zh-zl+1;
  myrank=Run.rank;

  for(int v=0;v<3;v++)
    next[v] = Grid.stride[v];

  int ndim=3,cart_periods[3];  

  MPI_Cart_get(cart_comm,ndim,cart_sizes,cart_periods,lrank);
  MPI_Cart_coords(cart_comm, Run.rank, ndim, lrank);
 
  for(int nleft,nright,nd=0;nd<ndim;nd++){
    MPI_Cart_shift(cart_comm,nd,1,&nleft,&nright);
    leftr[nd] = nleft;
    rightr[nd] = nright;
    
    isgbeg[nd]=Grid.is_gbeg[nd];
    isgend[nd]=Grid.is_gend[nd];
  }

  // colranks logic
  int ***colranks=i3dim(0,cart_sizes[2]-1,0,cart_sizes[1]-1,0,cart_sizes[0]-1);
  for(int j=0; j<cart_sizes[2];j++)
    for(int k=0; k<cart_sizes[1];k++)
      for(int i=0; i<cart_sizes[0];i++){
        int coords[3]={i,k,j};
        MPI_Cart_rank(cart_comm,coords,&(colranks[j][k][i]));
      }

  char carlson_fp[16];
  if (NMU == 3) strcpy(carlson_fp,"./carlson3.dat");
  else if (NMU == 6) strcpy(carlson_fp,"./carlson6.dat");
  else if (NMU == 10) strcpy(carlson_fp,"./carlson10.dat");
  else {
    fprintf(stdout," NMU %i not supported.\n ",NMU);
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  std::ifstream fptr(carlson_fp,std::ios::in);
  if(fptr){
    fptr.precision(16);
    for (int i=0; i < NMU; i++)
      fptr >> xmu[0][i] >> xmu[1][i] >> xmu[2][i] >> wmu[i];
    fptr.close();
  } else {
    fprintf(stdout,"carlson%i.dat not found. Aborting \n",NMU);
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  int* ranks=new int [cart_sizes[0]];
  MPI_Group MPI_GROUP_WORLD;

  MPI_Comm_group(MPI_COMM_WORLD,&MPI_GROUP_WORLD);
  MPI_Group ** grp_col = new MPI_Group * [cart_sizes[2]];
  comm_col = new MPI_Comm * [cart_sizes[2]];

  for(int j=0; j<cart_sizes[2];j++){
    grp_col[j]=new MPI_Group [cart_sizes[1]];
    comm_col[j] = new MPI_Comm [cart_sizes[1]];
  }

  for(int j=0; j< cart_sizes[2]; j++){
    for(int k=0; k< cart_sizes[1]; k++){
      for(int i=0; i< cart_sizes[0]; i++)
        ranks[i]=colranks[j][k][cart_sizes[0]-1-i];
      MPI_Group_incl(MPI_GROUP_WORLD,cart_sizes[0],ranks,&(grp_col[j][k]));
      MPI_Comm_create(MPI_COMM_WORLD,grp_col[j][k],&(comm_col[j][k]));
    }
  }

  load_bins(Run.kap_name);

  // Allocate MD arrays
  
  I_o = new double[nx*ny]; 
  I_band = d2dim(0, ny-1, 0, nx-1);
  
  Fr_mean = new double[Nbands];
  gFr_mean = new double[Nbands];
  memset(Fr_mean,0,Nbands*sizeof(double));
  memset(gFr_mean,0,Nbands*sizeof(double));

  lgTe = d3dim(0, ny-1, 0, nx-1, 0, nz-1);
  lgPe = d3dim(0, ny-1, 0, nx-1, 0, nz-1);
  T_ind = i3dim(0, ny-1, 0, nx-1, 0, nz-1);
  P_ind = i3dim(0, ny-1, 0, nx-1, 0, nz-1);
  rho = d3dim(0, ny-1, 0, nx-1, 0, nz-1);
  tr_switch = i3dim(0, ny-1, 0, nx-1, 0, nz-1);
  
  Tau = d3dim(0, ny-1, 0, nx-1, 0, nz-1);
  
  // Qt: (ny-yo)*(nx-xo)*(nz-zo)
  Qt = d3dim(0, ny-yo-1, 0, nx-xo-1, 0, nz-zo-1);

  Jt = d3dim(0, ny-1, 0, nx-1, 0, nz-1);
  St = d3dim(0, ny-1, 0, nx-1, 0, nz-1);

  sbuf = d2dim(0, ny-1, 0, nx-1); // [y][x]
  rbuf = d2dim(0, ny-1, 0, nx-1);
  
  // Qtemp
  Qtemp = d4dim(0, ny-1, 0, nx-1, 0, nz-1, 0, 1);

  B = d3dim(0, ny-1, 0, nx-1, 0, nz-1);
  kap = d3dim(0, ny-1, 0, nx-1, 0, nz-1);
  I_n = d3dim(0, ny-1, 0, nx-1, 0, nz-1);
  I_n1 = d3dim(0, ny-1, 0, nx-1, 0, nz-1);

  if (rttype==0){
      sig = kap;
      abn = kap;
  } else {
      sig = d3dim(0, ny-1, 0, nx-1, 0, nz-1);
      abn = d3dim(0, ny-1, 0, nx-1, 0, nz-1);
  }

  J_band = d3dim(0, ny-1, 0, nx-1, 0, nz-1);
  Fx = d3dim(0, ny-1, 0, nx-1, 0, nz-1);
  Fy = d3dim(0, ny-1, 0, nx-1, 0, nz-1);
  Fz = d3dim(0, ny-1, 0, nx-1, 0, nz-1);

  // coeff: [2][ny][nx][nz]
  coeff = d4dim(0, 1, 0, ny-1, 0, nx-1, 0, nz-1);
  coeff1 = d3dim(0, ny-1, 0, nx-1, 0, nz-1);
  coeff2 = d3dim(0, ny-1, 0, nx-1, 0, nz-1);

  // Col_out
  col_offz=Grid.beg[0]-Grid.gbeg[0];
  col_nz = Grid.lsize[0];
  col_nzt = Grid.gsize[0];
  col_nvar = 9;
  Col_out = d3dim(0, Nbands-1, 0, col_nz-1, 0, col_nvar-1);
  for(int b=0; b<Nbands; ++b)
      for(int k=0; k<col_nz; ++k)
          for(int v=0; v<col_nvar; ++v)
              Col_out[b][k][v] = 0.0;

  numits = i5dim(0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1);

  // x_sbuf etc.
  if (NDIM==3){
      y_sbuf = d7dim(0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, nx-1, 0, nz-1);
      y_rbuf = d7dim(0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, nx-1, 0, nz-1);
      y_oldbuf = d7dim(0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, nx-1, 0, nz-1);
  }

  if (NDIM>1){
      x_sbuf = d7dim(0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, ny-1, 0, nz-1);
      x_rbuf = d7dim(0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, ny-1, 0, nz-1);
      x_oldbuf = d7dim(0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, ny-1, 0, nz-1);
  }

  z_sbuf = d7dim(0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, ny-1, 0, nx-1);
  z_rbuf = d7dim(0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, ny-1, 0, nx-1);
  z_oldbuf = d7dim(0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, ny-1, 0, nx-1);

  double DX=Grid.dx[1];
  double DY=Grid.dx[2];
  double DZ=Grid.dx[0];
  
  // More explicit and intuitive index rule for coefficient calculation
  for(int l=0; l<NMU; ++l){
    // Calculate ratios of grid spacing to direction cosine
    double aM[] = {DX/xmu[0][l], DY/xmu[1][l], DZ/xmu[2][l]};
    
    // Determine the driving axis (base) which minimizes aM
    int base = 0; // Default to X axis (0)
    
    if (NDIM == 1) {
        base = 2; // In 1D, Z is the only axis
    } else if (NDIM == 2) {
        // In 2D (X, Z), compare Z vs X
        // Note: aM[2] <= aM[0] is checked.
        if (aM[2] <= aM[0]) base = 2;
        else base = 0;
    } else {
        // In 3D, compare all axes to find minimum
        if (aM[2] <= aM[0] && aM[2] <= aM[1]) base = 2;
        else if (aM[1] <= aM[0]) base = 1;
        else base = 0;
    }
    
    ibase[l] = base;
    ds_upw[l] = aM[base];
    
    // Axes perpendicular to base axis
    int u = (base + 1) % 3;
    int v = (base + 2) % 3;
    
    // Calculate interpolation weights based on ratios relative to base axis
    double r_u = aM[base] / aM[u];
    double r_v = aM[base] / aM[v];
    
    // Store coefficients for the selected base axis
    a_00[l] = (1.0 - r_u) * (1.0 - r_v);
    a_01[l] = r_v * (1.0 - r_u);
    a_10[l] = r_u * (1.0 - r_v);
    a_11[l] = r_u * r_v;
  }

  // Initialization of sbufs/oldbufs
  if (NDIM>1){
    for(int b=0; b<Nbands; ++b)
      for(int yy=0; yy<2; ++yy)
        for(int xx=0; xx<2; ++xx)
          for(int zz=0; zz<2; ++zz)
            for(int l=0; l<NMU; ++l)
              for(int i=0; i<ny; ++i) 
                for(int j=0; j<nz; ++j){
                    x_sbuf[b][yy][xx][zz][l][i][j] = 0.0;
                    x_rbuf[b][yy][xx][zz][l][i][j] = 0.0;
                    x_oldbuf[b][yy][xx][zz][l][i][j] = 0.0;
                }
  }

  if (NDIM==3){
    for(int b=0; b<Nbands; ++b)
      for(int yy=0; yy<2; ++yy)
        for(int xx=0; xx<2; ++xx)
          for(int zz=0; zz<2; ++zz)
            for(int l=0; l<NMU; ++l)
              for(int i=0; i<nx; ++i) 
                for(int j=0; j<nz; ++j){
                    y_sbuf[b][yy][xx][zz][l][i][j] = 0.0;
                    y_rbuf[b][yy][xx][zz][l][i][j] = 0.0;
                    y_oldbuf[b][yy][xx][zz][l][i][j] = 0.0;
                }
  }
  
  // z_sbuf init
  for(int b=0; b<Nbands; ++b)
    for(int yy=0; yy<2; ++yy)
      for(int xx=0; xx<2; ++xx)
        for(int zz=0; zz<2; ++zz)
          for(int l=0; l<NMU; ++l)
            for(int i=0; i<ny; ++i) 
              for(int j=0; j<nx; ++j){
                  z_sbuf[b][yy][xx][zz][l][i][j] = 0.0;
                  z_rbuf[b][yy][xx][zz][l][i][j] = 0.0;
                  z_oldbuf[b][yy][xx][zz][l][i][j] = 0.0;
              }
}

RTS::~RTS()
{
  for(int j = 0;j < cart_sizes[2]; j++) delete [] comm_col[j];
  delete [] comm_col;

  delete [] I_o;
  del_d2dim(I_band, 0, ny-1, 0, nx-1);
  
  del_i3dim(tr_switch, 0, ny-1, 0, nx-1, 0, nz-1);
  del_d3dim(lgTe, 0, ny-1, 0, nx-1, 0, nz-1);
  del_d3dim(lgPe, 0, ny-1, 0, nx-1, 0, nz-1);
  del_i3dim(T_ind, 0, ny-1, 0, nx-1, 0, nz-1);
  del_i3dim(P_ind, 0, ny-1, 0, nx-1, 0, nz-1);
  del_d3dim(rho, 0, ny-1, 0, nx-1, 0, nz-1);
  
  del_d3dim(Tau, 0, ny-1, 0, nx-1, 0, nz-1);
  del_d3dim(Qt, 0, ny-yo-1, 0, nx-xo-1, 0, nz-zo-1);
  del_d3dim(Jt, 0, ny-1, 0, nx-1, 0, nz-1);
  del_d3dim(St, 0, ny-1, 0, nx-1, 0, nz-1);

  del_d3dim(B, 0, ny-1, 0, nx-1, 0, nz-1);
  del_d3dim(kap, 0, ny-1, 0, nx-1, 0, nz-1);
  
  if (rttype==1){
      del_d3dim(sig, 0, ny-1, 0, nx-1, 0, nz-1);
      del_d3dim(abn, 0, ny-1, 0, nx-1, 0, nz-1);
  }
  
  del_d3dim(J_band, 0, ny-1, 0, nx-1, 0, nz-1);
  del_d3dim(I_n, 0, ny-1, 0, nx-1, 0, nz-1);
  del_d3dim(I_n1, 0, ny-1, 0, nx-1, 0, nz-1);

  del_d3dim(Fx, 0, ny-1, 0, nx-1, 0, nz-1);
  del_d3dim(Fy, 0, ny-1, 0, nx-1, 0, nz-1);
  del_d3dim(Fz, 0, ny-1, 0, nx-1, 0, nz-1);

  del_d4dim(coeff, 0, 1, 0, ny-1, 0, nx-1, 0, nz-1);
  del_d3dim(coeff1, 0, ny-1, 0, nx-1, 0, nz-1);
  del_d3dim(coeff2, 0, ny-1, 0, nx-1, 0, nz-1);

  del_d4dim(Qtemp, 0, ny-1, 0, nx-1, 0, nz-1, 0, 1);
  del_d2dim(sbuf, 0, ny-1, 0, nx-1);
  del_d2dim(rbuf, 0, ny-1, 0, nx-1);
  del_d3dim(Col_out, 0, Nbands-1, 0, col_nz-1, 0, col_nvar-1);
   
  del_i5dim(numits,0,Nbands-1,FWD,BWD,RIGHT,LEFT,UP,DOWN,0,NMU-1);
  
  if (NDIM>1){
    del_d7dim(x_sbuf, 0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, ny-1, 0, nz-1);
    del_d7dim(x_rbuf, 0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, ny-1, 0, nz-1);
    del_d7dim(x_oldbuf, 0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, ny-1, 0, nz-1);
  }
  if (NDIM==3){
    del_d7dim(y_sbuf, 0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, nx-1, 0, nz-1);
    del_d7dim(y_rbuf, 0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, nx-1, 0, nz-1);
    del_d7dim(y_oldbuf, 0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, nx-1, 0, nz-1);
  }
  
  del_d7dim(z_sbuf, 0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, ny-1, 0, nx-1);
  del_d7dim(z_rbuf, 0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, ny-1, 0, nx-1);
  del_d7dim(z_oldbuf, 0, Nbands-1, 0, 1, 0, 1, 0, 1, 0, NMU-1, 0, ny-1, 0, nx-1);
  
  ACCH::Delete(this, sizeof(RTS));
}

void RTS::load_bins(char* kap_name){
  std::ifstream fp_rt(kap_name, std::ios::in | std::ios::binary);

  int pt_rhot;
  int junk4;

  if (fp_rt.is_open()) {

    fp_rt.read((char*)&N5000, sizeof(int));
    fp_rt.read((char*)&NT, sizeof(int));
    fp_rt.read((char*)&Np, sizeof(int));
    fp_rt.read((char*)&Nbands, sizeof(int));
    fp_rt.read((char*)&pt_rhot, sizeof(int));
    fp_rt.read((char*)&fullodf, sizeof(int));
    fp_rt.read((char*)&scatter, sizeof(int));
    fp_rt.read((char*)&junk4, sizeof(int)); 
     
    if (pt_rhot!=0){
      fprintf(stdout,"pt_rhot = %i, but rho-T bins are not currently implemented, aborting",pt_rhot);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    if (myrank == 0) {
      fprintf(stdout,"Reading in RTS Header, current settings:\n");
      fprintf(stdout,"Full ODF is %d and scatter is %d \n",fullodf,scatter);
      fprintf(stdout,"RT bins are NT %d Np %d Nbands %d \n",NT,Np,Nbands);
      fprintf(stdout,"Reference bin is %d and coronal back heating bin %d \n",N5000,junk4);
    }
    
    if ((scatter==0)&&(eps_const==0)&&(rttype==1)){
      fprintf(stdout, "rt_type =%i, but scattering bins = %i and photon destruction probability is %e. \n",scatter, rttype,eps_const);
      fprintf(stdout, "For scattering either scattering bins, or constant photon destruction probability are required, aborting. \n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    tab_T = new double[NT];
    tab_p = new double[Np];

    invT_tab = new double[NT];
    invP_tab = new double[Np];
    
    fp_rt.read((char*)&tab_T[0],NT*sizeof(double));
    fp_rt.read((char*)&tab_p[0],Np*sizeof(double));
    
    for (int i=0;i<NT; i++)
      tab_T[i] = tab_T[i]*TENLOG;
     
    for (int j=0;j<Np; j++)
      tab_p[j] = tab_p[j]*TENLOG;

    for (int l=0; l<=NT-2; l++)
      invT_tab[l]= 1./(tab_T[l+1] - tab_T[l]);
    for (int m=0; m<=Np-2; m++)
      invP_tab[m] = 1./ (tab_p[m+1]- tab_p[m]);

    if (N5000){
      kap_5000_tab = new float[NT*Np];
      B_5000_tab = new float[NT];
       
      fp_rt.read((char*)&kap_5000_tab[0],NT*Np*sizeof(float));
      fp_rt.read((char*)&B_5000_tab[0],NT*sizeof(float));
    }
      
    kap_tab = new float[Nbands*NT*Np];
    B_tab = new float[Nbands*NT];
    fp_rt.read((char*)&kap_tab[0],Nbands*NT*Np*sizeof(float));
    fp_rt.read((char*)&B_tab[0],Nbands*NT*sizeof(float));

    if(fullodf){
      nu_tab = f1dim(0,Nlam-1);
      fp_rt.read((char*)&nu_tab[0],Nlam*sizeof(float));

      if (scatter>0){
        acont_pT = f3dim(0,Nlam-1,0,NT-1,0,Np-1);
        kcont_pT = f3dim(0,Nlam-1,0,NT-1,0,Np-1);
        fp_rt.read((char*)&acont_pT[0][0][0],Nlam*NT*Np*sizeof(float));
        fp_rt.read((char*)&kcont_pT[0][0][0],Nlam*NT*Np*sizeof(float));
        
        Npp=1;
        tau_pp_tab = d1dim(0,Npp-1);
        invtau_pp_tab = d1dim(0,Npp-1);

        tau_pp_tab[0] = -99.0;
        invtau_pp_tab[0] = 1.0;

        if (rttype==0)
          for (int lam = 0 ;lam < Nlam;lam++)
            for (int bin = 0;bin<Nbin;bin++)
              for (int tt = 0;tt<NT;tt++)
                for (int pp = 0;pp<Np;pp++)
                  kap_tab[((Nbin*lam+bin)*NT+tt)*Np+pp]
                        = log(exp(kap_tab[((Nbin*lam+bin)*NT+tt)*Np+pp])+acont_pT[lam][tt][pp]);
      }
    } else if (scatter>0){
        Npp = scatter;
        scatter=1;

        sig_tab = f3dim(0,Nbands-1,0,NT-1,0,Np-1);
        abn_tab = f3dim(0,Nbands-1,0,NT-1,0,Np-1);
        fp_rt.read((char*)&abn_tab[0][0][0],Nbands*NT*Np*sizeof(float));
        fp_rt.read((char*)&sig_tab[0][0][0],Nbands*NT*Np*sizeof(float));

        tau_pp_tab = d1dim(0,Npp-1);
        invtau_pp_tab = d1dim(0,Npp-1);

        fp_rt.read((char*)&tau_pp_tab[0],Npp*sizeof(double));

        for (int z=0;z<=Npp-2; z++)
          invtau_pp_tab[z] = 1./(tau_pp_tab[z+1]- tau_pp_tab[z]);
        
        kap_pp_tab = f2dim(0,Nbands-1,0,Npp-1);
        abn_pp_tab = f2dim(0,Nbands-1,0,Npp-1);
        sig_pp_tab = f2dim(0,Nbands-1,0,Npp-1);

        fp_rt.read((char*)&abn_pp_tab[0][0],Nbands*Npp*sizeof(float));
        fp_rt.read((char*)&kap_pp_tab[0][0],Nbands*Npp*sizeof(float));
        fp_rt.read((char*)&sig_pp_tab[0][0],Nbands*Npp*sizeof(float));
    }
        
    if ((scatter==0)&&(rttype==1)){
      Npp=1;
      tau_pp_tab = d1dim(0,Npp-1);
      invtau_pp_tab = d1dim(0,Npp-1);

      tau_pp_tab[0] = -99;
      invtau_pp_tab[0] = 1;
    }

    fp_rt.close();
     
  } else {
    cout << "rt_init: kappa file not found: " << kap_name << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

double RTS::tau(int z,int x,int y){
  double Tau_local = Tau[y][x][z] +
                     Tau[y-yo][x][z] +
                     Tau[y][x-xo][z] +
                     Tau[y][x][z-zo] +
                     Tau[y-yo][x-xo][z] +
                     Tau[y-yo][x][z-zo] +
                     Tau[y][x-xo][z-zo] +
                     Tau[y-yo][x-xo][z-zo];
  return Tau_local * 0.125;
}

double RTS::Qtot(int z,int x,int y)
{  
  return Qt[y-yl-yo][x-xl-xo][z-zl-zo];
}

double RTS::Jtot(int z,int x,int y){
  double J_local = Jt[y][x][z] +
                   Jt[y-yo][x][z] +
                   Jt[y][x-xo][z] +
                   Jt[y][x][z-zo] +
                   Jt[y-yo][x-xo][z] +
                   Jt[y-yo][x][z-zo] +
                   Jt[y][x-xo][z-zo] +
                   Jt[y-yo][x-xo][z-zo];
  return J_local * 0.125;
}

double RTS::Stot(int z,int x,int y){
  double S_local = St[y][x][z] +
                   St[y-yo][x][z] +
                   St[y][x-xo][z] +
                   St[y][x][z-zo] +
                   St[y-yo][x-xo][z] +
                   St[y-yo][x][z-zo] +
                   St[y][x-xo][z-zo] +
                   St[y-yo][x-xo][z-zo];
  return S_local * 0.125;
}

void RTS::UpdateIout()
{
}

double RTS::Iout(int x,int y)
{
  double io = I_o[(y-yl)*nx+(x-xl)] + I_o[(y-yo-yl)*nx+(x-xl)] + I_o[(y-yl)*nx+(x-xo-xl)] + I_o[(y-yo-yl)*nx+(x-xo-xl)];
  io *= 0.25;
  return io;
}

double RTS::wrapper(int rt_upd,GridData &Grid,RunData &Run,const PhysicsData &Physics){
  for(int y = 0; y < ny; y++)
    for(int x = 0; x < nx; x++) {
      I_band[y][x] = 0.0; 
      I_o[y*nx+x] = 0.0;
    }

  double DX=Grid.dx[1],DZ=Grid.dx[0],DY=Grid.dx[2];
  int cont_bin = Physics.rt[i_rt_iout];

  const double Temp_TR = Physics.rt[i_rt_tr_tem];
  const double Pres_TR = Physics.rt[i_rt_tr_pre];

  need_I = 0;

  if (Run.NeedsSlice())
    need_I = 1;
  
  if((Run.iteration%rt_upd)&&(dt_rad>0)){
    if (cont_bin == 2){
      for (int band=Nbands-1;band>=0;--band){
        get_Tau_and_Iout(Grid, Run, Physics,DZ,&B_tab[band*NT],&kap_tab[band*NT*Np],(double***)I_band,need_I);
        for (int y=0;y<ny;y++)
          for (int x=0;x<nx;x++)
            I_o[y*nx+x] +=I_band[y][x];
      }
    }

    if (cont_bin==0){
      get_Tau_and_Iout(Grid, Run, Physics,DZ,B_tab,kap_tab,(double***)I_band,need_I);
      for (int y=0;y<ny;y++)
        for (int x=0;x<nx;x++)
          I_o[y*nx+x] +=I_band[y][x];
    }

    if (N5000){
      int I5000_out = 0;
      if ((cont_bin==1)&&(need_I==1))
        I5000_out = 1;
      
      get_Tau_and_Iout(Grid, Run, Physics,DZ,B_5000_tab,kap_5000_tab,(double***)I_band,I5000_out);
     
      if ((cont_bin==1)&&(need_I==1)) {
        for (int y=0;y<ny;y++)
          for (int x=0;x<nx;x++)
            I_o[y*nx+x] += I_band[y][x];
      }
    }

    calc_Qtot_and_Tau(Grid, Run, Physics);
    return dt_rad;
  }

  cState *U=Grid.U;
  double N = pow(2,NDIM);

  for(int y=0;y<ny;++y){
    for(int x=0;x<nx;++x){
      int y0 = y+yl;
      int x0 = x+xl;
      int off0 = x0*next[1]+y0*next[2];
      int off1 = x0*next[1]+(y0+yo)*next[2];
      int off2 = (x0+xo)*next[1]+y0*next[2];
      int off3 = (x0+xo)*next[1]+(y0+yo)*next[2];
      
      for(int z=0;z<nz;++z){
        int z0 = z+zl;
        int inode[]={off0+z0,off0+z0+zo,off2+z0,off2+z0+zo,off1+z0,off1+z0+zo,off3+z0,off3+z0+zo};
        double Tm=0.0,pm=0.0,rm=0.0;
        for(int l=0;l<N;++l){
          Tm+=Grid.temp[inode[l]];
          pm+=Grid.pres[inode[l]];
          rm+=U[inode[l]].d;
        }

        Tm /= N;
        pm /= N;
        rm /= N;

        double pswitch = min(max(pm-Pres_TR,0.0),1.0);
        double tswitch = min(max(Temp_TR-Tm,0.0),1.0);
        tr_switch[y][x][z] = (int) max(pswitch,tswitch);
        
        lgTe[y][x][z]=log(Tm);
        lgPe[y][x][z]=log(pm);
        rho[y][x][z] =rm;
    
        lgTe[y][x][z] = min(max(lgTe[y][x][z],(double) tab_T[0]),(double) tab_T[NT-1]);
        lgPe[y][x][z] = min(max(lgPe[y][x][z],(double) tab_p[0]),(double) tab_p[Np-1]);
      }
      
      for(int z=0;z<nz;++z){
        int l=0;
        int m=0;
        if(lgTe[y][x][z]<tab_T[0])
          l=0;
        else if(lgTe[y][x][z]>tab_T[NT-1])
          l=NT-2;
        else {
          for (int li=0; li<=NT-2; li++){
            int lflag = 0;
            if ((lgTe[y][x][z] >= tab_T[li]) && (lgTe[y][x][z] <= tab_T[li+1]))
              lflag = lflag+ 1;
            if(lflag==1)
              l = li;
          }
        }

        if(lgPe[y][x][z]<tab_p[0])
          m=0;
        else if(lgPe[y][x][z]>tab_p[Np-1])
          m=Np-2;
        else {
          for (int mi=0; mi<=Np-2; mi++){
             int mflag = 0;
             if ((lgPe[y][x][z] >= tab_p[mi]) && (lgPe[y][x][z] <= tab_p[mi+1]))
               mflag = mflag+1;
             if(mflag==1)
               m=mi;
          }
        }

        T_ind[y][x][z] = l;
        P_ind[y][x][z] = m;

      }
    }
  }
  
  for(int y=0; y<ny; y++)
    for(int x=0; x<nx; x++)
      for(int z=0; z<nz; z++) {
        St[y][x][z] = 0.0;
        Jt[y][x][z] = 0.0;
      }

  for(int y = 0; y < ny-yo; y++)
    for(int x = 0; x < nx-xo; x++)
      for(int z = 0; z < nz-zo; z++) {
        Qt[y][x][z] = 0.0;
      }

  for(int band=Nbands-1;band>=0;--band){

    if(fullodf){
      nu_ind = band/Nbin;
      bin_ind = band%Nbin;
    }

    if((myrank==0) && (verbose>1)){
      fprintf(stdout,"rt running for band %i of %i \n",band+1,Nbands);
      if (fullodf)
        fprintf(stdout,"nu band = %i, %e and bin %i \n", nu_ind, nu_tab[nu_ind],bin_ind);
    }

    for(int y=0;y<ny;++y){
      for(int x=0;x<nx;++x){
        for(int z=0;z<nz;++z){
      
          int l = T_ind[y][x][z];
          int m = P_ind[y][x][z];
      
          double xt = (lgTe[y][x][z]-tab_T[l])*invT_tab[l];
          double xp = (lgPe[y][x][z]-tab_p[m])*invP_tab[m];
      
          B[y][x][z]=exp(xt*B_tab[band*NT+l+1]+(1.-xt)*B_tab[band*NT+l]);
      
          kap[y][x][z] = 
            exp(xt*(xp*kap_tab[(band*NT+l+1)*Np+m+1]+(1.-xp)*kap_tab[(band*NT+l+1)*Np+m])+
            (1.-xt)*(xp*kap_tab[(band*NT+l)*Np+m+1]+(1.-xp)*kap_tab[(band*NT+l)*Np+m]));
    
          kap[y][x][z] *= tr_switch[y][x][z];
          B[y][x][z]   *= tr_switch[y][x][z];
        }
      }
    }

  if(isgbeg[0]==1) {
    for(int YDIR=FWD;YDIR<=BWD;++YDIR) {
      for(int XDIR=RIGHT;XDIR<=LEFT;++XDIR) {
        for(int l=0;l<NMU;++l) {
          for(int y=0;y<ny;++y) {
            for(int x=0;x<nx;++x) {
              z_rbuf[band][YDIR][XDIR][UP][l][y][x]=B[y][x][0];
            }
	  }
	}
      }
    }
  }

  if(isgend[0]==1) {
    for(int YDIR=FWD;YDIR<=BWD;++YDIR) {
      for(int XDIR=RIGHT;XDIR<=LEFT;++XDIR) {
        for(int l=0; l<NMU; l++) {
            for(int y=0; y<ny; y++)
                for(int x=0; x<nx; x++)
                    z_rbuf[band][YDIR][XDIR][DOWN][l][y][x] = 0.0;
	    }
      }
    }
  }

  driver(DZ,DX,DY,band); 

  gFr_mean[band] = 0.0;

  if(isgend[0]==1){
    double gFr_mean_reduc = 0.0;
    for(int y=0; y<ny-yo; y++)
      for(int x=0; x<nx-xo; x++) {
        gFr_mean_reduc+=Fz[y][x][nz-1] +
                        Fz[y+yo][x][nz-1] +
                        Fz[y][x+xo][nz-1] +
                        Fz[y+yo][x+xo][nz-1];
      }
    gFr_mean[band]=gFr_mean_reduc*0.25;
  }
 
  if (((cont_bin==1)&&(band==0))||((cont_bin==0)&&(band==1))) {
    for(int y=0;y<ny;++y)
      for(int x=0;x<nx;++x)
        I_o[y*nx+x] = 0.0;
  }

  }

  MPI_Allreduce(&gFr_mean[0],&Fr_mean[0],Nbands,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  F_o = 0.0;
 

  for (int band=0;band<Nbands;++band){
    Fr_mean[band]/=(Grid.gsize[1]*Grid.gsize[2]);
    F_o+=Fr_mean[band];
  }
  
  if (N5000){
    int I5000_out = 0;
    if ((cont_bin==1)&&(need_I==1))
      I5000_out = 1;
    
    get_Tau_and_Iout(Grid, Run, Physics,DZ,B_5000_tab,kap_5000_tab,(double***)I_band,I5000_out);
   
    if ((cont_bin==1)&&(need_I==1)) {
      for (int y=0;y<ny;y++)
        for (int x=0;x<nx;x++)
          I_o[y*nx+x] += I_band[y][x];
    }
  }
  calc_Qtot_and_Tau(Grid, Run, Physics);
  
  if (Run.NeedsSlice() && Run.RT_HAVG)
    save_1D_avg(Run.path_2D,Run.globiter,Run.time); 

  return dt_rad;
}

void RTS::calc_Qtot_and_Tau(GridData &Grid, const RunData &Run, const PhysicsData &Physics){

 double tau_min = pow(Physics.rt[i_rt_tau_min],2);
 cState *U = Grid.U;

  dt_rad=0.0;
  double _dt_rad = 0.0;
  double qsum=0.0;
  double _qsum=0.0;

  for(int y = yo; y < ny; y++)
    for(int x = xo; x < nx; x++) {
      int off0 = (x+xl)*next[1]+(y+yl)*next[2];
      for(int z = zo; z < nz; z++) {
        int node = off0+z+zl;
	Grid.Tau[node] = tau(z, x, y);
	Grid.Jtot[node] = Jtot(z, x, y);
	Grid.Stot[node] = Stot(z, x, y);
	
	double scale = pow(Grid.Tau[node], 2);
	scale = scale/(scale + tau_min)*tr_switch[y][x][z];
        double Qt_step = Qt[y-yo][x-xo][z-zo]*scale;
	double inv_dt = fabs(Qt_step)/U[node].e;
	_dt_rad = max(_dt_rad, inv_dt);
	Grid.Qtot[node] = Qt_step;
	_qsum += Qt_step;
      }
    }
  dt_rad = _dt_rad;
  qsum = _qsum;

  exchange_single_acc(Grid,Grid.Tau);

  double Fqrad;

  MPI_Allreduce(&qsum,&Fqrad,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  Fqrad=-Fqrad*Grid.dx[0]/(Grid.gsize[1]*Grid.gsize[2]);

  if(myrank==0) fprintf(stdout,"RT energy flux: %21.15E %21.15E\n",F_o,Fqrad);
  
  F_o=Fqrad;
  dt_rad=Physics.rt[i_rt_cfl]/dt_rad;
  
}

void RTS::get_Tau_and_Iout(GridData &Grid, const RunData &Run, const PhysicsData &Physics, double DZ, float * B_Iout_tab, float * kap_Iout_tab, double *** I_band, int calc_int){
  
  double ttime=MPI_Wtime(),stime=0.0,atime=0.0;
  
  double N = pow(2,Grid.NDIM);

  const double Temp_TR = Physics.rt[i_rt_tr_tem];
  const double Pres_TR = Physics.rt[i_rt_tr_pre];

  for(int y = 0; y < ny; y++) {
    for(int x = 0; x < nx; x++) {
      int off0 = (x+xl)*next[1]+(y+yl)*next[2];
      int off1 = (x+xl)*next[1]+(y+yo+yl)*next[2];
      int off2 = (x+xo+xl)*next[1]+(y+yl)*next[2];
      int off3 = (x+xo+yl)*next[1]+(y+yo+yl)*next[2];
      for(int z = 0; z < nz; z++) {
        int inode[]={off0+z+zl,off0+z+zo+zl,off2+z+zl,off2+z+zo+zl,off1+z+zl,off1+z+zo+zl,off3+z+zl,off3+z+zo+zl};
        double lgP = 0, lgT = 0, rm = 0;
        for(int l=0;l<N;++l){
          lgP+=Grid.pres[inode[l]];
          lgT+=Grid.temp[inode[l]];
          rm +=Grid.U[inode[l]].d;
        }
        lgP /= N;
        lgT /= N;
        rm  /= N;

        double pswitch  = min(max(lgP-Pres_TR,0.0),1.0);
        double tswitch  = min(max(Temp_TR-lgT,0.0),1.0);
        double trswitch = max(pswitch,tswitch);

        lgP          = log(lgP);
        lgT          = log(lgT);
        rho[y][x][z] = rm;

        lgT = min(max(lgT, (double) tab_T[0]),(double) tab_T[NT-1]);
        lgP = min(max(lgP, (double) tab_p[0]),(double) tab_p[Np-1]);

        int l=0;
        int m=0;
        if(lgT<tab_T[0])
          l=0;
        else if(lgT>tab_T[NT-1])
          l=NT-2;
        else {
        for (int li=0; li<=NT-2; li++){
            int lflag = 0;
            if ((lgT >= tab_T[li]) && (lgT <= tab_T[li+1])){
              lflag = lflag+1;
            if(lflag==1)
              l = li;
            }
          }
        }

        if(lgP<tab_p[0])
          m=0;
        else if(lgP>tab_p[Np-1])
          m=Np-2;
        else {
        for (int mi=0; mi<=Np-2; mi++){
            int mflag = 0;
            if ((lgP >= tab_p[mi]) && (lgP <= tab_p[mi+1]))
              mflag = mflag+1;
            if(mflag==1)
              m = mi;
          }
        }

        double xt = (lgT-tab_T[l])*invT_tab[l];
        double xp = (lgP-tab_p[m])*invP_tab[m];

        B[y][x][z]=exp(xt*B_Iout_tab[l+1]+(1.-xt)*B_Iout_tab[l]);

        kap[y][x][z]=
                    exp(xt*(xp*kap_Iout_tab[Np*(l+1)+(m+1)]+
                           (1.-xp)*kap_Iout_tab[Np*(l+1)+m])+
                           (1.-xt)*(xp*kap_Iout_tab[Np*l+(m+1)]+
                           (1.-xp)*kap_Iout_tab[(Np*l+m)]));

        B[y][x][z]   *= trswitch;
        kap[y][x][z] *= trswitch;

      }
    }
  }

  for(int y = 0; y < ny; y++)
    for(int x = 0; x < nx; x++) {
      Tau[y][x][nz-1]=1.0e-12 * ((double) isgend[0]);
      for(int z = nz-2; z >= 0; z--) {
        double k0=kap[y][x][z],r0=rho[y][x][z],k_upw=kap[y][x][z+1],r_upw=rho[y][x][z+1];
        double delta_tau=DZ*((k0*r0+k_upw*r_upw)*inv3+(k0*r_upw+k_upw*r0)*inv6);
        Tau[y][x][z]=Tau[y][x][z+1]+delta_tau;
      }
    }

  for(int y=0;y<ny;++y){ // loop over RT grid
    for(int x=0;x<nx;++x){
      rbuf[y][x]=0.e0;
      sbuf[y][x]=Tau[y][x][0];
    }
  }

  double ctime=MPI_Wtime();
  MPI_Exscan(sbuf[0], rbuf[0], nx*ny, MPI_DOUBLE, MPI_SUM, comm_col[lrank[2]][lrank[1]]);
  stime+=MPI_Wtime()-ctime;
  
  if(isgend[0] == 0) {
    for(int y=0;y<ny;++y){ // loop over RT grid
      for(int x=0;x<nx;++x){
        double rbufyx = rbuf[y][x];
        for(int z=0;z<nz;++z){
          Tau[y][x][z]+=rbufyx;
        }
      }
    }
  }

  if (calc_int){
  //  Outgoing Intensity at top (Long Characteristics)
  for(int y=0;y<ny;++y){ // loop over RT grid
    for(int x=0;x<nx;++x){
      rbuf[y][x]=0.0;
      sbuf[y][x]=0.0;
      double tmp = 0.0;
      for(int z=1;z<nz;++z){
        double Ss1 = B[y][x][z];
        double Ss2 = B[y][x][z-1];
        double delta_tau=Tau[y][x][z-1]-Tau[y][x][z];
        if(delta_tau>dtau_min){
          double edt=exp(-delta_tau);
          double c1=(1.0-edt)/delta_tau;
          tmp+=(Ss1*(1.0-c1)+Ss2*(c1-edt))*exp(-Tau[y][x][z]);
        }else{
          tmp+=0.5*delta_tau*(Ss1+Ss2);
        }
      }
      sbuf[y][x] = tmp;
    }
  }
  ctime=MPI_Wtime();
  MPI_Allreduce(sbuf[0], rbuf[0], nx*ny, MPI_DOUBLE, MPI_SUM, comm_col[lrank[2]][lrank[1]]);
  atime+=MPI_Wtime()-ctime;
  for(int y=0;y<ny;++y){
    for(int x=0;x<nx;++x){
      I_band[y][x][0]+=rbuf[y][x];
    }
  }
  }

  ttime=MPI_Wtime()-ttime;
  if((myrank==0)&&(verbose>2)) printf("tau5000 time: %f %f %f %f \n",ttime,stime,atime,(stime+atime)/ttime);
  
}    

void RTS::tauscale_qrad(int band, double DX,double DY,double DZ, double *** Ss){

  double ttime=MPI_Wtime(),stime=0.0,atime=0.0;
  double idx=1.0/DX,idy=1.0/DY,idz=1.0/DZ;

  if(NDIM==1) {
    idx=0.;
    idy=0.;
  }
  if(NDIM==2) {
    idy=0.0;
  }

  for(int y=0;y<ny;y++){ 
    for(int x=0;x<nx;x++){
      Tau[y][x][nz-1]=1.0e-12  * ((double) isgend[0]);
      for(int z=nz-2;z>=0;--z){
	double k0=kap[y][x][z],r0=rho[y][x][z],k_upw=kap[y][x][z+1],r_upw=rho[y][x][z+1];
	Tau[y][x][z]=Tau[y][x][z+1] + (DZ*((k0*r0+k_upw*r_upw)*inv3+(k0*r_upw+k_upw*r0)*inv6));
      }
      
      rbuf[y][x]=0.e0;
      sbuf[y][x]=Tau[y][x][0];
    }
  }

  double ctime=MPI_Wtime();
  MPI_Exscan(sbuf[0], rbuf[0], nx*ny, MPI_DOUBLE, MPI_SUM, comm_col[lrank[2]][lrank[1]]);
  stime+=MPI_Wtime()-ctime;
  if (isgend[0] == 0){
    for(int y=0;y<ny;++y){ 
      for(int x=0;x<nx;++x){
        double rbufyx = rbuf[y][x];
        for(int z=0;z<nz;++z){
          Tau[y][x][z]+=rbufyx;
        }
      }
    }
  }

  if (need_I){
    for(int y=0;y<ny;++y){ 
      for(int x=0;x<nx;++x){
        rbuf[y][x]=0.0;
        sbuf[y][x]=0.0;
        double tmp = 0.0;
        for(int z=1;z<nz;++z){
          double Ss1 = Ss[y][x][z];
          double Ss2 = Ss[y][x][z-1];
          double delta_tau=Tau[y][x][z-1]-Tau[y][x][z];
          if(delta_tau>dtau_min){
            double edt=exp(-delta_tau);
            double c1=(1.0-edt)/delta_tau;
            tmp+=(Ss1*(1.0-c1)+Ss2*(c1-edt))*exp(-Tau[y][x][z]);
          }else{
            tmp+=0.5*delta_tau*(Ss1+Ss2);
          }  
        }
        sbuf[y][x] = tmp;
      }
    }
    ctime=MPI_Wtime(); 
    MPI_Allreduce(sbuf[0], rbuf[0], nx*ny, MPI_DOUBLE, MPI_SUM, comm_col[lrank[2]][lrank[1]]);
    ctime+=MPI_Wtime()-ctime;
    for(int y=0;y<ny;++y){
      for(int x=0;x<nx;++x){
        I_o[y*nx+x]+=rbuf[y][x];
      }
    }

  }
  for(int y=0; y<ny; y++)
    for(int x=0; x<nx; x++)
      for(int z=0; z<nz; z++) {
         I_n[y][x][z] = kap[y][x][z]*rho[y][x][z]*(J_band[y][x][z]-Ss[y][x][z]);
         St[y][x][z] += Ss[y][x][z];
         Jt[y][x][z] += J_band[y][x][z];
      }

  double inv_tau_0=1.0e1;
  for(int y=0;y<ny-yo;y++){
    for(int x=0;x<nx-xo;x++){
      for(int z=0;z<nz-zo;z++){
        double qf1=
          (
            (
              Fz[y][x][z+zo]+
              Fz[y][x+xo][z+zo]+
              Fz[y+yo][x][z+zo]+
              Fz[y+yo][x+xo][z+zo]
            )-(
              Fz[y][x][z]+
              Fz[y][x+xo][z]+
              Fz[y+yo][x][z]+
              Fz[y+yo][x+xo][z]
            )
          )*idz

          +

          (
            (
              Fx[y][x+xo][z]+
              Fx[y][x+xo][z+zo]+
              Fx[y+yo][x+xo][z]+
              Fx[y+yo][x+xo][z+zo]
            )-(
              Fx[y][x][z]+
              Fx[y][x][z+zo]+
              Fx[y+yo][x][z]+
              Fx[y+yo][x][z+zo]
            )
          )*idx

          +

          (
            (
              Fy[y+yo][x][z]+
              Fy[y+yo][x][z+zo]+
              Fy[y+yo][x+xo][z]+
              Fy[y+yo][x+xo][z+zo]
            )-(
              Fy[y][x][z]+
              Fy[y][x][z+zo]+
              Fy[y][x+xo][z]+
              Fy[y][x+xo][z+zo]
            )
          )*idy;
        qf1*=-0.25e0; 
        double qj1 = I_n[y][x][z] +
                     I_n[y][x][z+zo] +
                     I_n[y][x+xo][z] +
                     I_n[y][x+xo][z+zo] +
                     I_n[y+yo][x][z] +
                     I_n[y+yo][x][z+zo] +
                     I_n[y+yo][x+xo][z] +
                     I_n[y+yo][x+xo][z+zo];
        qj1*=0.5*PI; 
        double tau_local=tau(z+zo,x+xo,y+yo);
        double weight=exp(-tau_local*inv_tau_0);
        
        Qt[y][x][z]+=weight*qj1+(1.0-weight)*qf1;
      }
    }
  }

   if (save_col){
     int col_bnd2 = col_bnd[2];
     int col_bnd1 = col_bnd[1];
     
     if ( (col_bnd2>=yl) && (col_bnd2<=yh) && (col_bnd1>=xl) && (col_bnd1<=xh)){
       for (int k=0;k<col_nz;k++){
         Col_out[band][k][0] = J_band[col_bnd2-yl][col_bnd1-xl][k+col_offz];
         Col_out[band][k][1] = Ss[col_bnd2-yl][col_bnd1-xl][k+col_offz];
         Col_out[band][k][2] = kap[col_bnd2-yl][col_bnd1-xl][k+col_offz]; 
         Col_out[band][k][5] = B[col_bnd2-yl][col_bnd1-xl][k+col_offz];
         Col_out[band][k][6] = Tau[col_bnd2-yl][col_bnd1-xl][k+col_offz];
         Col_out[band][k][7] = I_n[col_bnd2-yl][col_bnd1-xl][k+col_offz];
         
       }
     }
   }

}

void RTS::driver(double DZ, double DX, double DY, int band){
  double etime=0.0,atime=0.0,cmp_time1=0.0,cmp_time2=0.0,buf_time=0.0,err_time=0.0,flx_time=0.0,tau_time=0.0; 
  double ttime=MPI_Wtime();
  
  int stepvec[3][4][3] = { {{1,0,0},{1,0,1},{1,1,0},{1,1,1}},
               {{0,1,0},{1,1,0},{0,1,1},{1,1,1}},
               {{0,0,1},{0,1,1},{1,0,1},{1,1,1}} };

  for(int y = 0; y < ny; y++)
    for(int x = 0; x < nx; x++)
      for(int z = 0; z < nz; z++) {
        I_n[y][x][z] = 0.0;
        Fz[y][x][z] = 0.0;
        Fx[y][x][z] = 0.0;
        Fy[y][x][z] = 0.0;
        J_band[y][x][z] = 0.0;
      }

  double maxerr_up=0.0,maxerr_down=0.0;
  double itavg=0.0;
  double aravg=0.0;

  for(int ZDIR=UP;ZDIR<=DOWN;++ZDIR){
    int zi_i=(ZDIR==UP)?zl+1:zh-1,zi_f=(ZDIR==UP)?zh:zl,zstep=(ZDIR==UP)?1:-1;
    for(int XDIR=RIGHT;XDIR<=LEFT;++XDIR){
      int xi_i=(XDIR==RIGHT)?xl+1:xh-1,xi_f=(XDIR==RIGHT)?xh:xl,xstep=(XDIR==RIGHT)?1:-1;
      for(int YDIR=FWD;YDIR<=BWD;++YDIR){
        int yi_i=(YDIR==FWD)?yl+1:yh-1,yi_f=(YDIR==FWD)?yh:yl,ystep=(YDIR==FWD)?1:-1;
          for(int l=0;l<NMU;++l){
            double I_min=max(1.0,threshold*Fr_mean[band]/(NMU*pow(2,NDIM)));
            double c[]={a_00[l],a_01[l],a_10[l],a_11[l]};
        for(int m=0;m<=3;++m){
          ixstep[m]=stepvec[ibase[l]][m][0]*xstep;
          iystep[m]=stepvec[ibase[l]][m][1]*ystep;
          izstep[m]=stepvec[ibase[l]][m][2]*zstep;
        }

        double stime=MPI_Wtime();
        interpol(zi_i,zi_f,zstep,xi_i,xi_f,xstep,yi_i,yi_f,ystep,l,coeff,B);
        cmp_time1+=MPI_Wtime()-stime;
        stime=MPI_Wtime(); 
            int rt_iter=0;
            int rt_min_iter=(call_count<2)?0:numits[band][YDIR][XDIR][ZDIR][l]-(!(call_count%3));            
        double gmaxerr=1.0E10;
        while(gmaxerr>=threshold){
          rt_iter+=1;
          itavg+=1.0;
          stime = MPI_Wtime();
          readbuf(band,l,ZDIR,XDIR,YDIR);
          buf_time += MPI_Wtime()-stime;
          stime=MPI_Wtime();
          integrate(c);
	  cmp_time2 += MPI_Wtime()-stime;
          stime = MPI_Wtime();
          writebuf(band,l,ZDIR,XDIR,YDIR); 
          buf_time += MPI_Wtime()-stime;
          stime = MPI_Wtime();
          exchange(band, l, ZDIR, XDIR,YDIR);
          etime += MPI_Wtime()-stime;
          if(rt_iter>=rt_min_iter){
            stime = MPI_Wtime();
            double err=error(band,l,ZDIR,XDIR,YDIR,I_min);
            err_time += MPI_Wtime()-stime;
            stime = MPI_Wtime();
            MPI_Allreduce(&err,&gmaxerr,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
            aravg += 1.0;
            atime += MPI_Wtime()-stime;
          }
        }
        if(ZDIR==UP){
          maxerr_up=max(maxerr_up,gmaxerr);
        }else{
          maxerr_down=max(maxerr_down,gmaxerr);
        }
        numits[band][YDIR][XDIR][ZDIR][l]=rt_iter;
        stime=MPI_Wtime();
        readbuf(band,l,ZDIR,XDIR,YDIR);
        buf_time += MPI_Wtime()-stime;
        stime=MPI_Wtime();
        flux(l,ZDIR,XDIR,YDIR);
        flx_time += MPI_Wtime()-stime;

      }
      }
    }
  }
  double stime=MPI_Wtime();
  tauscale_qrad(band,DX,DY,DZ,B);
  tau_time+=MPI_Wtime()-stime;
  if((myrank==0) && (verbose>1)){
    fprintf(stdout,"rt_driver iter : %f %f \n",aravg/(8.0*NMU),itavg/(8.0*NMU));
    fprintf(stdout,"rt_driver error: %e %e \n",maxerr_up,maxerr_down);
  }

  ttime=MPI_Wtime()-ttime;  
  if((myrank==0) && (verbose>2))
    fprintf(stdout,"rt_driver time : %f %f %f %f %f %f %f %f %f %f \n",ttime,cmp_time1,
       cmp_time2,buf_time,err_time,flx_time,tau_time,etime,atime,(etime+atime)/ttime);
   

  call_count+=1;
}

void RTS::interpol(int zi_i,int zi_f,int zstep,int xi_i,int xi_f,int xstep,
           int yi_i,int yi_f,int ystep,int l,double**** coeff, double *** B)
{
  double ds3=ds_upw[l]*inv3,ds6=ds_upw[l]*inv6;
  double c[]={a_00[l],a_01[l],a_10[l],a_11[l]};
  double c0 = c[0], c1 = c[1], c2 = c[2], c3 = c[3];

  if(NDIM==3){
    for(int y = 0; y < ny-1; y++) {
      for(int x = 0; x < nx-1; x++) {
        for(int z = 0; z < nz-1; z++) {
          int yi = (yi_i-yl) + y*ystep;
          int xi = (xi_i-xl) + x*xstep;
          int zi = (zi_i-zl) + z*zstep;
          
          double _r_upw=
            c0*rho[yi-iystep[0]][xi-ixstep[0]][zi-izstep[0]]+
            c1*rho[yi-iystep[1]][xi-ixstep[1]][zi-izstep[1]]+
            c2*rho[yi-iystep[2]][xi-ixstep[2]][zi-izstep[2]]+
            c3*rho[yi-iystep[3]][xi-ixstep[3]][zi-izstep[3]];
          
          double _k_upw=
            c0*kap[yi-iystep[0]][xi-ixstep[0]][zi-izstep[0]]+
            c1*kap[yi-iystep[1]][xi-ixstep[1]][zi-izstep[1]]+
            c2*kap[yi-iystep[2]][xi-ixstep[2]][zi-izstep[2]]+
            c3*kap[yi-iystep[3]][xi-ixstep[3]][zi-izstep[3]];

          double _S_upw=
            c0*B[yi-iystep[0]][xi-ixstep[0]][zi-izstep[0]]+
            c1*B[yi-iystep[1]][xi-ixstep[1]][zi-izstep[1]]+
            c2*B[yi-iystep[2]][xi-ixstep[2]][zi-izstep[2]]+
            c3*B[yi-iystep[3]][xi-ixstep[3]][zi-izstep[3]];

          double _r0=rho[yi][xi][zi];
          double _k0=kap[yi][xi][zi];
          double _S0=B[yi][xi][zi];

          double dt=ds3*(_k_upw*_r_upw+_k0*_r0)+ds6*(_k0*_r_upw+_k_upw*_r0);
          double expo=exp(-dt);
          double w0,w1;
          if (dt > dtau_min){
            w0=1.0-expo;
            w1=w0-dt*expo;
          }else{
            w0=dt-dt*dt/2.0+dt*dt*dt/6.0;
            w1=dt*dt/2.0-dt*dt*dt/3.0;
          }
          double source=_S0*(w0-w1/dt)+_S_upw*(w1/dt);

          if (dt > dtau_min2){
            coeff[0][y][x][z] = expo;
            coeff[1][y][x][z] = source;
          }else{
            coeff[0][y][x][z] = 1.0; 
            coeff[1][y][x][z] = 0.0;
          }
        }
      }
    }
  }
}

void RTS::integrate(const double c[4])
{
  int xstart, xend, xinc;
  int ystart, yend, yinc;
  int zstart, zend, zinc;
  
  xinc = 0;
  for(int m=0; m<4; m++) if(ixstep[m] != 0) xinc = ixstep[m];
  if(xinc == 0) xinc = 1; 
  else xinc = (xinc > 0) ? 1 : -1;
  
  yinc = 0;
  for(int m=0; m<4; m++) if(iystep[m] != 0) yinc = iystep[m];
  if(yinc == 0) yinc = 1;
  else yinc = (yinc > 0) ? 1 : -1;
  
  zinc = 0;
  for(int m=0; m<4; m++) if(izstep[m] != 0) zinc = izstep[m];
  if(zinc == 0) zinc = 1;
  else zinc = (zinc > 0) ? 1 : -1;
  
  if (xinc == 1) { xstart=0; xend=nx; } else { xstart=nx-1; xend=-1; }
  if (yinc == 1) { ystart=0; yend=ny; } else { ystart=ny-1; yend=-1; }
  if (zinc == 1) { zstart=0; zend=nz; } else { zstart=nz-1; zend=-1; }
  
  double c0=c[0], c1=c[1], c2=c[2], c3=c[3];
  
  int ys = ystart, ye = yend;
  int xs = xstart, xe = xend;
  int zs = zstart, ze = zend;
  
  if (yinc == 1) { ys=1; ye=ny; } else { ys=ny-2; ye=-1; }
  if (xinc == 1) { xs=1; xe=nx; } else { xs=nx-2; xe=-1; }
  if (zinc == 1) { zs=1; ze=nz; } else { zs=nz-2; ze=-1; }
  
  for(int cy=0; cy<ny-1; ++cy) {
      int yy = ys + cy*yinc;
      for(int cx=0; cx<nx-1; ++cx) {
          int xx = xs + cx*xinc;
          for(int cz=0; cz<nz-1; ++cz) {
              int zz = zs + cz*zinc;
              
              double I_upw = 
                c0 * I_n[yy-iystep[0]][xx-ixstep[0]][zz-izstep[0]] +
                c1 * I_n[yy-iystep[1]][xx-ixstep[1]][zz-izstep[1]] +
                c2 * I_n[yy-iystep[2]][xx-ixstep[2]][zz-izstep[2]] +
                c3 * I_n[yy-iystep[3]][xx-ixstep[3]][zz-izstep[3]];
              
              I_n[yy][xx][zz] = I_upw * coeff[0][cy][cx][cz] + coeff[1][cy][cx][cz];
          }
      }
  }
}

void RTS::readbuf(int band,int l,int  DIR,int XDIR,int YDIR){
  if (NDIM==3){
      if (YDIR==FWD){ 
          for(int x=0; x<nx; ++x)
            for(int z=0; z<nz; ++z)
               I_n[0][x][z] = y_rbuf[band][YDIR][XDIR][DIR][l][x][z];
      } else { 
          for(int x=0; x<nx; ++x)
            for(int z=0; z<nz; ++z)
               I_n[ny-1][x][z] = y_rbuf[band][YDIR][XDIR][DIR][l][x][z];
      }
      
      if (XDIR==RIGHT){ 
          for(int y=0; y<ny; ++y)
            for(int z=0; z<nz; ++z)
               I_n[y][0][z] = x_rbuf[band][YDIR][XDIR][DIR][l][y][z];
      } else {
          for(int y=0; y<ny; ++y)
            for(int z=0; z<nz; ++z)
               I_n[y][nx-1][z] = x_rbuf[band][YDIR][XDIR][DIR][l][y][z];
      }
      
      if (DIR==UP){ 
          for(int y=0; y<ny; ++y)
            for(int x=0; x<nx; ++x)
               I_n[y][x][0] = z_rbuf[band][YDIR][XDIR][DIR][l][y][x];
      } else {
          for(int y=0; y<ny; ++y)
            for(int x=0; x<nx; ++x)
               I_n[y][x][nz-1] = z_rbuf[band][YDIR][XDIR][DIR][l][y][x];
      }
  }
}

void RTS::writebuf(int band, int l,int DIR,int XDIR,int YDIR){
  if (NDIM==3){
      if (YDIR==FWD){
          for(int x=0; x<nx; ++x)
            for(int z=0; z<nz; ++z)
               y_sbuf[band][YDIR][XDIR][DIR][l][x][z] = I_n[ny-1][x][z];
      } else {
          for(int x=0; x<nx; ++x)
            for(int z=0; z<nz; ++z)
               y_sbuf[band][YDIR][XDIR][DIR][l][x][z] = I_n[0][x][z];
      }
      
      if (XDIR==RIGHT){
          for(int y=0; y<ny; ++y)
            for(int z=0; z<nz; ++z)
               x_sbuf[band][YDIR][XDIR][DIR][l][y][z] = I_n[y][nx-1][z];
      } else {
          for(int y=0; y<ny; ++y)
            for(int z=0; z<nz; ++z)
               x_sbuf[band][YDIR][XDIR][DIR][l][y][z] = I_n[y][0][z];
      }
      
      if (DIR==UP){
          for(int y=0; y<ny; ++y)
            for(int x=0; x<nx; ++x)
               z_sbuf[band][YDIR][XDIR][DIR][l][y][x] = I_n[y][x][nz-1];
      } else {
          for(int y=0; y<ny; ++y)
            for(int x=0; x<nx; ++x)
               z_sbuf[band][YDIR][XDIR][DIR][l][y][x] = I_n[y][x][0];
      }
  }
}

void RTS::flux(int l,int DIR,int XDIR,int YDIR){
  double w=wmu[l];
  double mux=xmu[0][l];
  double muy=xmu[1][l];
  double muz=xmu[2][l];
  
  for(int y=0; y<ny; y++)
    for(int x=0; x<nx; x++)
      for(int z=0; z<nz; z++) {
          double In = I_n[y][x][z];
          Fx[y][x][z] += w*mux*In;
          Fy[y][x][z] += w*muy*In;
          Fz[y][x][z] += w*muz*In;
          J_band[y][x][z] += w*In;
      }
}

double RTS::error(int band,int l,int ZDIR,int XDIR,int YDIR,double I_min){
  double max_err = 0.0;
  
  if(NDIM==3){
      if(YDIR==FWD){
          for(int x=0; x<nx; ++x) for(int z=0; z<nz; ++z){
             double val = y_sbuf[band][YDIR][XDIR][ZDIR][l][x][z];
             double old = y_oldbuf[band][YDIR][XDIR][ZDIR][l][x][z];
             max_err = max(max_err, fabs(val-old)/I_min);
             y_oldbuf[band][YDIR][XDIR][ZDIR][l][x][z] = val;
          }
      } else {
          for(int x=0; x<nx; ++x) for(int z=0; z<nz; ++z){
             double val = y_sbuf[band][YDIR][XDIR][ZDIR][l][x][z];
             double old = y_oldbuf[band][YDIR][XDIR][ZDIR][l][x][z];
             max_err = max(max_err, fabs(val-old)/I_min);
             y_oldbuf[band][YDIR][XDIR][ZDIR][l][x][z] = val;
          }
      }
      
      if(XDIR==RIGHT){
          for(int y=0; y<ny; ++y) for(int z=0; z<nz; ++z){
             double val = x_sbuf[band][YDIR][XDIR][ZDIR][l][y][z];
             double old = x_oldbuf[band][YDIR][XDIR][ZDIR][l][y][z];
             max_err = max(max_err, fabs(val-old)/I_min);
             x_oldbuf[band][YDIR][XDIR][ZDIR][l][y][z] = val;
          }
      } else {
          for(int y=0; y<ny; ++y) for(int z=0; z<nz; ++z){
             double val = x_sbuf[band][YDIR][XDIR][ZDIR][l][y][z];
             double old = x_oldbuf[band][YDIR][XDIR][ZDIR][l][y][z];
             max_err = max(max_err, fabs(val-old)/I_min);
             x_oldbuf[band][YDIR][XDIR][ZDIR][l][y][z] = val;
          }
      }
      
      if(ZDIR==UP){
          for(int y=0; y<ny; ++y) for(int x=0; x<nx; ++x){
             double val = z_sbuf[band][YDIR][XDIR][ZDIR][l][y][x];
             double old = z_oldbuf[band][YDIR][XDIR][ZDIR][l][y][x];
             max_err = max(max_err, fabs(val-old)/I_min);
             z_oldbuf[band][YDIR][XDIR][ZDIR][l][y][x] = val;
          }
      } else {
          for(int y=0; y<ny; ++y) for(int x=0; x<nx; ++x){
             double val = z_sbuf[band][YDIR][XDIR][ZDIR][l][y][x];
             double old = z_oldbuf[band][YDIR][XDIR][ZDIR][l][y][x];
             max_err = max(max_err, fabs(val-old)/I_min);
             z_oldbuf[band][YDIR][XDIR][ZDIR][l][y][x] = val;
          }
      }
  }
  return max_err;
}

void RTS::exchange(int band,int l,int DIR,int XDIR,int YDIR){
  MPI_Request req[6];
  MPI_Status stat[6];
  int nreq=0;
  
  int left = leftr[2];
  int right = rightr[2];
  int count = nx*nz;
  
  if (YDIR==FWD){
      MPI_Isend(y_sbuf[band][YDIR][XDIR][DIR][l][0], count, MPI_DOUBLE, right, 0, comm_col[lrank[2]][lrank[1]], &req[nreq++]);
      MPI_Irecv(y_rbuf[band][YDIR][XDIR][DIR][l][0], count, MPI_DOUBLE, left, 0, comm_col[lrank[2]][lrank[1]], &req[nreq++]);
  } else {
      MPI_Isend(y_sbuf[band][YDIR][XDIR][DIR][l][0], count, MPI_DOUBLE, left, 1, comm_col[lrank[2]][lrank[1]], &req[nreq++]);
      MPI_Irecv(y_rbuf[band][YDIR][XDIR][DIR][l][0], count, MPI_DOUBLE, right, 1, comm_col[lrank[2]][lrank[1]], &req[nreq++]);
  }
  
  left = leftr[1];
  right = rightr[1];
  count = ny*nz;
  if (XDIR==RIGHT){
      MPI_Isend(x_sbuf[band][YDIR][XDIR][DIR][l][0], count, MPI_DOUBLE, right, 2, comm_col[lrank[2]][lrank[1]], &req[nreq++]);
      MPI_Irecv(x_rbuf[band][YDIR][XDIR][DIR][l][0], count, MPI_DOUBLE, left, 2, comm_col[lrank[2]][lrank[1]], &req[nreq++]);
  } else {
      MPI_Isend(x_sbuf[band][YDIR][XDIR][DIR][l][0], count, MPI_DOUBLE, left, 3, comm_col[lrank[2]][lrank[1]], &req[nreq++]);
      MPI_Irecv(x_rbuf[band][YDIR][XDIR][DIR][l][0], count, MPI_DOUBLE, right, 3, comm_col[lrank[2]][lrank[1]], &req[nreq++]);
  }
  
  left = leftr[0];
  right = rightr[0];
  count = ny*nx;
  if (DIR==UP){
      MPI_Isend(z_sbuf[band][YDIR][XDIR][DIR][l][0], count, MPI_DOUBLE, right, 4, comm_col[lrank[2]][lrank[1]], &req[nreq++]);
      MPI_Irecv(z_rbuf[band][YDIR][XDIR][DIR][l][0], count, MPI_DOUBLE, left, 4, comm_col[lrank[2]][lrank[1]], &req[nreq++]);
  } else {
      MPI_Isend(z_sbuf[band][YDIR][XDIR][DIR][l][0], count, MPI_DOUBLE, left, 5, comm_col[lrank[2]][lrank[1]], &req[nreq++]);
      MPI_Irecv(z_rbuf[band][YDIR][XDIR][DIR][l][0], count, MPI_DOUBLE, right, 5, comm_col[lrank[2]][lrank[1]], &req[nreq++]);
  }
  
  MPI_Waitall(nreq, req, stat);
}

void RTS::IntegrateSetup(int yi_i, int xi_i, int zi_i, int ystep, int xstep, int zstep) {}
void RTS::interpol_and_integrate(const double c[4], double *** Ss, int l) {}
