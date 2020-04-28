#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define pi 3.141592


double dis(double*,double*,int,int);
double ang(double*,double*,double*,int,int,int);

const double avogadro=6.0221367e+23;
const double kb=1.380648e-23; 

enum {
  FCC_PARTICLE_NUM = 8,
  FCC_LATTICE_NUM = 3,
  PARTICLE_NUM= FCC_PARTICLE_NUM*FCC_LATTICE_NUM*FCC_LATTICE_NUM*FCC_LATTICE_NUM
};

int main()
{
  int i,j,k,l,m;
  int ix,iy,iz;
  int repeat_x=2;//蜴溷ｭ仙・譛滄・鄂ｮ郢ｰ繧願ｿ斐＠謨ｰx
  int repeat_y=2;//蜴溷ｭ仙・譛滄・鄂ｮ郢ｰ繧願ｿ斐＠謨ｰy
  int repeat_z=2;//蜴溷ｭ仙・譛滄・鄂ｮ郢ｰ繧願ｿ斐＠謨ｰz
  int loop=2000;//繝｡繧､繝ｳ繝ｫ繝ｼ繝礼ｹｰ繧願ｿ斐＠謨ｰ
  int count=0;//繝｡繧､繝ｳ繝ｫ繝ｼ繝励き繧ｦ繝ｳ繝・
  int par;//邊貞ｭ先焚
  int num_j[PARTICLE_NUM*18];
  int num_k[PARTICLE_NUM*18];

  double si[PARTICLE_NUM*3*18],cl[3];//螟我ｽ搾ｼｻﾃ・ｼｽ
  double v_si[PARTICLE_NUM*3*18],v_cl[3];//騾溷ｺｦ・ｻm/s・ｽ
  double v;
  double f_si[PARTICLE_NUM*3*18],f_cl[3];//蜉幢ｼｻN・ｽ
  double t,dt=0.5e-15;//譎る俣
  double angs=1e-10;//ﾃ・｣懈ｭ｣
  double ev=1.60218e-19;//N陬懈ｭ｣
  double T0=293.15,T,x;
  double co;
  double ri,rj,rk;
  double fi,fj,fk;
  //tersoff繝代Λ繝｡繝ｼ繧ｿ
  double lam1=2.4799;
  double lam2=1.7322;
  double beta=1.0999e-6;
  double n=7.8737e-1;
  double c=1.0039e+5;
  double d=1.6218e+1;
  double h=-5.9826e-1;
  double R=2.85;
  double D=0.15;
  double A=1.8308e+3;
  double B=4.7118e+2;
  double dz[3],dz_i[3],dz_j[3],dz_k[3];
  double fr,fa,fc_ij,fc_ik,dfc_ij,dfc_ik,b,db,g,dg,z;
  double r_ij,r_ik,r_jk;
  //蜴溷ｭ舌ヱ繝ｩ繝｡繝ｼ繧ｿ
  double m_si=28.0855/avogadro*1e-3;//[kg]
  double m_cl=35.453/avogadro*1e-3;//[kg]


  //繝ｻ繝ｻ繝ｻ繝ｻ繝ｻ蛻晄悄譚｡莉ｶ繝ｻ繝ｻ繝ｻ繝ｻ繝ｻ
  //繝ｻ繝ｻ繝ｻ蛻晄悄菴咲ｽｮ繝ｻ繝ｻ繝ｻ
  //Si
  i=0;
  for(iz=0;iz<repeat_z;iz++){
    for(iy=0;iy<repeat_y;iy++){
      for(ix=0;ix<repeat_x;ix++){
	si[i]  =0+5.43*ix;
	si[i+1]=0+5.43*iy;
	si[i+2]=0+5.43*iz;
	i+=3;
	
	si[i]  =0+5.43*ix;
	si[i+1]=2.715+5.43*iy;
	si[i+2]=2.715+5.43*iz;
	i+=3;

	si[i]  =2.715+5.43*ix;
	si[i+1]=2.715+5.43*iy;
	si[i+2]=0+5.43*iz;
	i+=3;

	si[i]  =2.715+5.43*ix;
	si[i+1]=0+5.43*iy;
	si[i+2]=2.715+5.43*iz;
	i+=3;

	si[i]  =4.0725+5.43*ix;
	si[i+1]=1.3575+5.43*iy;
	si[i+2]=4.0725+5.43*iz;
	i+=3;

	si[i]  =1.3575+5.43*ix;
	si[i+1]=1.3575+5.43*iy;
	si[i+2]=1.3575+5.43*iz;
	i+=3;

	si[i]  =1.3575+5.43*ix;
	si[i+1]=4.0725+5.43*iy;
	si[i+2]=4.0725+5.43*iz;
	i+=3;

	si[i]  =4.0725+5.43*ix;
	si[i+1]=4.0725+5.43*iy;
	si[i+2]=1.3525+5.43*iz;
	i+=3;
      }
    }
  }

  par=i/3;//邊貞ｭ先焚

  //繝溘Λ繝ｼ繧ｻ繝ｫ縺ｫ繧ｳ繝斐・
  for(i=0;i<par*3;i+=3){
    si[i+(par*3)]  =si[i]+5.43*repeat_x;
    si[i+(par*3)+1]=si[i+1];
    si[i+(par*3)+2]=si[i+2];	

    si[i+(par*3*2)]  =si[i];
    si[i+(par*3*2)+1]=si[i+1]+5.43*repeat_y;
    si[i+(par*3*2)+2]=si[i+2];

    si[i+(par*3*3)]  =si[i]-5.43*repeat_x;
    si[i+(par*3*3)+1]=si[i+1];
    si[i+(par*3*3)+2]=si[i+2];

    si[i+(par*3*4)]  =si[i];
    si[i+(par*3*4)+1]=si[i+1]-5.43*repeat_y;
    si[i+(par*3*4)+2]=si[i+2];

    si[i+(par*3*5)]  =si[i]+5.43*repeat_x;
    si[i+(par*3*5)+1]=si[i+1]+5.43*repeat_y;
    si[i+(par*3*5)+2]=si[i+2];

    si[i+(par*3*6)]  =si[i]+5.43*repeat_x;
    si[i+(par*3*6)+1]=si[i+1]-5.43*repeat_y;
    si[i+(par*3*6)+2]=si[i+2];

    si[i+(par*3*7)]  =si[i]-5.43*repeat_x;
    si[i+(par*3*7)+1]=si[i+1]+5.43*repeat_y;
    si[i+(par*3*7)+2]=si[i+2];

    si[i+(par*3*8)]  =si[i]-5.43*repeat_x;
    si[i+(par*3*8)+1]=si[i+1]-5.43*repeat_y;
    si[i+(par*3*8)+2]=si[i+2];

    si[i+(par*3*9)]  =si[i]+5.43*repeat_x;
    si[i+(par*3*9)+1]=si[i+1];
    si[i+(par*3*9)+2]=si[i+2]-5.43*repeat_z;	

    si[i+(par*3*10)]  =si[i];
    si[i+(par*3*10)+1]=si[i+1]+5.43*repeat_y;
    si[i+(par*3*10)+2]=si[i+2]-5.43*repeat_z;

    si[i+(par*3*11)]  =si[i]-5.43*repeat_x;
    si[i+(par*3*11)+1]=si[i+1];
    si[i+(par*3*11)+2]=si[i+2]-5.43*repeat_z;

    si[i+(par*3*12)]  =si[i];
    si[i+(par*3*12)+1]=si[i+1]-5.43*repeat_y;
    si[i+(par*3*12)+2]=si[i+2]-5.43*repeat_z;

    si[i+(par*3*13)]  =si[i]+5.43*repeat_x;
    si[i+(par*3*13)+1]=si[i+1]+5.43*repeat_y;
    si[i+(par*3*13)+2]=si[i+2]-5.43*repeat_z;

    si[i+(par*3*14)]  =si[i]+5.43*repeat_x;
    si[i+(par*3*14)+1]=si[i+1]-5.43*repeat_y;
    si[i+(par*3*14)+2]=si[i+2]-5.43*repeat_z;

    si[i+(par*3*15)]  =si[i]-5.43*repeat_x;
    si[i+(par*3*15)+1]=si[i+1]+5.43*repeat_y;
    si[i+(par*3*15)+2]=si[i+2]-5.43*repeat_z;

    si[i+(par*3*16)]  =si[i]-5.43*repeat_x;
    si[i+(par*3*16)+1]=si[i+1]-5.43*repeat_y;
    si[i+(par*3*16)+2]=si[i+2]-5.43*repeat_z;

    si[i+(par*3*17)]  =si[i];
    si[i+(par*3*17)+1]=si[i+1];
    si[i+(par*3*17)+2]=si[i+2]-5.43*repeat_z;


  }

  //Cl
  cl[0]=5.43*repeat_x/2;
  cl[1]=5.43*repeat_y/2;
  cl[2]=13;

  //繝ｻ繝ｻ繝ｻ蛻晄悄騾溷ｺｦ繝ｻ繝ｻ繝ｻ
  //Si
  for(i=0;i<par*3;i++){
    v_si[i]=0;
  }

  //Cl
  v_cl[0]=0;
  v_cl[1]=0;
  v_cl[2]=-sqrt(3*kb*T0/m_cl);

  //繝ｻ繝ｻ繝ｻ蜉帙・蛻晄悄蛹悶・繝ｻ繝ｻ
  //Si
  for(i=0;i<par*3;i++){
    f_si[i]=0;
  }

  //Cl
  f_cl[0]=0;
  f_cl[1]=0;
  f_cl[2]=0;


  //繝ｻ繝ｻ繝ｻ繝ｻ繝ｻ繝｡繧､繝ｳ繝ｫ繝ｼ繝励・繝ｻ繝ｻ繝ｻ繝ｻ
  for(t=0;t<dt*loop;t+=dt,count++){

    //繝ｻ繝ｻ繝ｻ騾溷ｺｦ縺ｮ譖ｴ譁ｰ・茨ｼ托ｼ峨・繝ｻ繝ｻ
    //Si
    for(i=0;i<par*3;i+=3){
      v_si[i]=v_si[i]+dt/2*f_si[i]/m_si*ev/angs;
      v_si[i+1]=v_si[i+1]+dt/2*f_si[i+1]/m_si*ev/angs;
      v_si[i+2]=v_si[i+2]+dt/2*f_si[i+2]/m_si*ev/angs;
    }
    
    //Cl
    v_cl[0]=v_cl[0]+dt/2*f_cl[0]/m_cl*ev/angs;
    v_cl[1]=v_cl[1]+dt/2*f_cl[1]/m_cl*ev/angs;
    v_cl[2]=v_cl[2]+dt/2*f_cl[2]/m_cl*ev/angs;
    
    //繝ｻ繝ｻ繝ｻ菴咲ｽｮ縺ｮ譖ｴ譁ｰ繝ｻ繝ｻ繝ｻ
    //Si
    for(i=0;i<par*3;i+=3){
      si[i]=si[i]+dt*v_si[i]/angs;
      si[i+1]=si[i+1]+dt*v_si[i+1]/angs;
      si[i+2]=si[i+2]+dt*v_si[i+2]/angs;
    }
    //Cl
    cl[0]=cl[0]+dt*v_cl[0]/angs;
    cl[1]=cl[1]+dt*v_cl[1]/angs;
    cl[2]=cl[2]+dt*v_cl[2]/angs;
      
      
    //繝ｻ繝ｻ繝ｻ蜻ｨ譛滓擅莉ｶ繝ｻ繝ｻ繝ｻ
    //蝓ｺ譛ｬ繧ｻ繝ｫ縺ｮ螟悶↓蜃ｺ縺溽ｲ貞ｭ舌・菫ｮ豁｣
    for(i=0;i<par*3;i+=3){
      if(si[i]<0){
	while(si[i]<0){
	  si[i]+=5.43*repeat_x;
	}
      }
      if(si[i+1]<0){
	while(si[i+1]<0){
	  si[i+1]+=5.43*repeat_y;
	}
      }
      if(si[i]>5.43*repeat_x){
	while(si[i]>5.43*repeat_x){
	  si[i]-=5.43*repeat_x;
	}
      }
      if(si[i+1]>5.43*repeat_y){
	while(si[i+1]>5.43*repeat_y){
	  si[i+1]-=5.43*repeat_y;
	}
      }
    }

    //繝溘Λ繝ｼ繧ｻ繝ｫ縺ｫ繧ｳ繝斐・
    for(i=0;i<par*3;i+=3){
      si[i+(par*3)]  =si[i]+5.43*repeat_x;
      si[i+(par*3)+1]=si[i+1];
      si[i+(par*3)+2]=si[i+2];	

      si[i+(par*3*2)]  =si[i];
      si[i+(par*3*2)+1]=si[i+1]+5.43*repeat_y;
      si[i+(par*3*2)+2]=si[i+2];

      si[i+(par*3*3)]  =si[i]-5.43*repeat_x;
      si[i+(par*3*3)+1]=si[i+1];
      si[i+(par*3*3)+2]=si[i+2];

      si[i+(par*3*4)]  =si[i];
      si[i+(par*3*4)+1]=si[i+1]-5.43*repeat_y;
      si[i+(par*3*4)+2]=si[i+2];

      si[i+(par*3*5)]  =si[i]+5.43*repeat_x;
      si[i+(par*3*5)+1]=si[i+1]+5.43*repeat_y;
      si[i+(par*3*5)+2]=si[i+2];

      si[i+(par*3*6)]  =si[i]+5.43*repeat_x;
      si[i+(par*3*6)+1]=si[i+1]-5.43*repeat_y;
      si[i+(par*3*6)+2]=si[i+2];

      si[i+(par*3*7)]  =si[i]-5.43*repeat_x;
      si[i+(par*3*7)+1]=si[i+1]+5.43*repeat_y;
      si[i+(par*3*7)+2]=si[i+2];

      si[i+(par*3*8)]  =si[i]-5.43*repeat_x;
      si[i+(par*3*8)+1]=si[i+1]-5.43*repeat_y;
      si[i+(par*3*8)+2]=si[i+2];

      //si[i+(par*3*9)]  =si[i]+5.43*repeat_x;
      //si[i+(par*3*9)+1]=si[i+1];
      //si[i+(par*3*9)+2]=si[i+2]-5.43*repeat_z;	

      //si[i+(par*3*10)]  =si[i];
      //si[i+(par*3*10)+1]=si[i+1]+5.43*repeat_y;
      //si[i+(par*3*10)+2]=si[i+2]-5.43*repeat_z;

      //si[i+(par*3*11)]  =si[i]-5.43*repeat_x;
      //si[i+(par*3*11)+1]=si[i+1];
      //si[i+(par*3*11)+2]=si[i+2]-5.43*repeat_z;

      //si[i+(par*3*12)]  =si[i];
      //si[i+(par*3*12)+1]=si[i+1]-5.43*repeat_y;
      //si[i+(par*3*12)+2]=si[i+2]-5.43*repeat_z;

      //si[i+(par*3*13)]  =si[i]+5.43*repeat_x;
      //si[i+(par*3*13)+1]=si[i+1]+5.43*repeat_y;
      //si[i+(par*3*13)+2]=si[i+2]-5.43*repeat_z;

      //si[i+(par*3*14)]  =si[i]+5.43*repeat_x;
      //si[i+(par*3*14)+1]=si[i+1]-5.43*repeat_y;
      //si[i+(par*3*14)+2]=si[i+2]-5.43*repeat_z;

      //si[i+(par*3*15)]  =si[i]-5.43*repeat_x;
      //si[i+(par*3*15)+1]=si[i+1]+5.43*repeat_y;
      //si[i+(par*3*15)+2]=si[i+2]-5.43*repeat_z;

      //si[i+(par*3*16)]  =si[i]-5.43*repeat_x;
      //si[i+(par*3*16)+1]=si[i+1]-5.43*repeat_y;
      //si[i+(par*3*16)+2]=si[i+2]-5.43*repeat_z;

      //si[i+(par*3*17)]  =si[i];
      //si[i+(par*3*17)+1]=si[i+1];
      //si[i+(par*3*17)+2]=si[i+2]-5.43*repeat_z;


    }

    //繝ｻ繝ｻ繝ｻ蜉帙・蛻晄悄蛹悶・繝ｻ繝ｻ
    //Si
    for(i=0;i<par*3;i++){
      f_si[i]=0;
    }

    //Cl
    f_cl[0]=0;
    f_cl[1]=0;
    f_cl[2]=0;

    //繝ｻ繝ｻ繝ｻ蜉帙・險育ｮ励・繝ｻ繝ｻ
    //Si
    for(i=0;i<par*3*18;i+=3){

      for(j=0,l=0;j<par*3*18;j+=3){ //繧ｫ繝・ヨ繧ｪ繝戊ｷ晞屬蜀・・蜴溷ｭ舌・菫晉ｮ｡(j)
	if(j!=i){
	  if(dis(si,si,i,j)<R+D){
	    num_j[l]=j;
	    l++;
	  }
	}
      }
      
      for(k=0,m=0;k<par*3*18;k+=3){ //繧ｫ繝・ヨ繧ｪ繝戊ｷ晞屬蜀・・蜴溷ｭ舌・菫晉ｮ｡(k)
	if(k!=i && k!=j){
	  if(dis(si,si,i,k)<R+D){
	    num_k[m]=k;
	    m++;
	  }
	}
      }


      for(j=0;j<l;j++){
	r_ij=dis(si,si,i,num_j[j]);
	fr=A*exp(-lam1*r_ij);
	fa=-B*exp(-lam2*r_ij);

	if(R-D<r_ij){
	  fc_ij=1/2*(1-sin(pi/2*(r_ij-R)/D));
	  dfc_ij=-pi/(4*D)*cos(pi/2*(r_ij-R)/D);
	}
	
	else{
	  fc_ij=1;
	  dfc_ij=0;
	}
	z=0;
	for(k=0;k<m;k++){//ﾎｶ縺ｮ險育ｮ・
	  if(num_j[j]!=num_k[k]){
	    r_ik=dis(si,si,i,num_k[k]);
	    if(R-D<r_ik){
	      fc_ik=1/2*(1-sin(pi/2*(r_ik-R)/D));
	      dfc_ik=-pi/(4*D)*cos(pi/2*(r_ik-R)/D);
	    }
	    else{
	      fc_ik=1;
	      dfc_ik=0;
	    }
	    g=1+c*c/d/d-c*c/(d*d+pow(h-ang(si,si,si,i,num_j[j],num_k[k]),2));
	    z+=fc_ik*g;
	  }
	}

	b=1+pow(beta,n)*pow(z,n);

	if(z!=0)db=-1/2*pow(beta,n)*pow(z,n-1)*pow(b,-1/(2*n)-1);
	else db=0;

	b=pow(b,-1/(2*n));


	fi=-(dfc_ij*fr-lam1*fc_ij*fr)*(si[num_j[j]]-si[i])/r_ij;
	fi+=-b*(dfc_ij*fa-lam2*fc_ij*fa)*(si[num_j[j]]-si[i])/r_ij;
	fi*=0.5;

	fj=-(dfc_ij*fr-lam1*fc_ij*fr)*(si[num_j[j]+1]-si[i+1])/r_ij;
	fj+=-b*(dfc_ij*fa-lam2*fc_ij*fa)*(si[num_j[j]+1]-si[i+1])/r_ij;
	fj*=0.5;

	fk=-(dfc_ij*fr-lam1*fc_ij*fr)*(si[num_j[j]+2]-si[i+2])/r_ij;
	fk+=-b*(dfc_ij*fa-lam2*fc_ij*fa)*(si[num_j[j]+2]-si[i+2])/r_ij;
	fk*=0.5;
	  
	f_si[i]-=fi;
	f_si[i+1]-=fj;
	f_si[i+2]-=fk;
	  
	f_si[num_j[j]]+=fi;
	f_si[num_j[j]+1]+=fj;
	f_si[num_j[j]+2]+=fk;
	  
	for(k=0;k<m;k++){
	  if(num_j[j]!=num_k[k]){

	    r_ik=dis(si,si,i,num_k[k]);

	    if(R-D<r_ik){
	      fc_ik=1/2*(1-sin(pi/2*(r_ik-R)/D));
	      dfc_ik=-pi/(4*D)*cos(pi/2*(r_ik-R)/D);
	    }
	    else{
	      fc_ik=1;
	      dfc_ik=0;
	    }

	    co=ang(si,si,si,i,num_j[j],num_k[k]);
	    g=1+c*c/d/d-c*c/(d*d+pow(h-co,2));
	    dg=-2*c*c*(h-co)/pow(d*d+pow(h-co,2),2);

	    dz[0]=g*dfc_ik*(si[num_k[k]]-si[i])/r_ik;
	    dz[1]=g*dfc_ik*(si[num_k[k]+1]-si[i+1])/r_ik;
	    dz[2]=g*dfc_ik*(si[num_k[k]+2]-si[i+2])/r_ik;

	    //
	    dz_i[0]=dz[0]+fc_ik*dg*((co/r_ij-1/r_ik)*(si[i]-si[num_j[j]])/r_ij+(co/r_ik-1/r_ij)*(si[i]-si[num_k[k]])/r_ik);
	    dz_i[1]=dz[1]+fc_ik*dg*((co/r_ij-1/r_ik)*(si[i+1]-si[num_j[j]+1])/r_ij+(co/r_ik-1/r_ij)*(si[i+1]-si[num_k[k]+1])/r_ik);
	    dz_i[2]=dz[2]+fc_ik*dg*((co/r_ij-1/r_ik)*(si[i+2]-si[num_j[j]+2])/r_ij+(co/r_ik-1/r_ij)*(si[i+2]-si[num_k[k]+2])/r_ik);

	    f_si[i]+=fc_ij*fa*db*dz_i[0];
	    f_si[i+1]+=fc_ij*fa*db*dz_i[1];
	    f_si[i+2]+=fc_ij*fa*db*dz_i[2];

	    //
	    dz_j[0]=fc_ik*dg/r_ij*(-co*(si[i]-si[num_j[j]])/r_ij+(si[i]-si[num_k[k]])/r_ik);
	    dz_j[1]=fc_ik*dg/r_ij*(-co*(si[i+1]-si[num_j[j]+1])/r_ij+(si[i+1]-si[num_k[k]+1])/r_ik);
	    dz_j[2]=fc_ik*dg/r_ij*(-co*(si[i+2]-si[num_j[j]+2])/r_ij+(si[i+2]-si[num_k[k]+2])/r_ik);

	    f_si[num_j[j]]+=fc_ij*fa*db*dz_j[0];
	    f_si[num_j[j]+1]+=fc_ij*fa*db*dz_j[1];
	    f_si[num_j[j]+2]+=fc_ij*fa*db*dz_j[2];
      
	    //
	    dz_k[0]=-dz[0]+fc_ik*dg/r_ik*(-co*(si[i]-si[num_k[k]])/r_ik+(si[i]-si[num_j[j]])/r_ij);
	    dz_k[0]=-dz[0]+fc_ik*dg/r_ik*(-co*(si[i]-si[num_k[k]])/r_ik+(si[i]-si[num_j[j]])/r_ij);
	    dz_k[0]=-dz[0]+fc_ik*dg/r_ik*(-co*(si[i]-si[num_k[k]])/r_ik+(si[i]-si[num_j[j]])/r_ij);

	    f_si[num_k[k]]+=fc_ij*fa*db*dz_k[0];
	    f_si[num_k[k]+1]+=fc_ij*fa*db*dz_k[1];
	    f_si[num_k[k]+2]+=fc_ij*fa*db*dz_k[2];    

	  }
	}
      }
    }
  
    //Cl
    for(i=0;i<par*3*18;i+=3){
      r_ij=dis(cl,si,0,i);
      if(r_ij<3.5){

	fi=-(-3631.65+8776.68*r_ij-8700*pow(r_ij,2)+4527.67*pow(r_ij,3)-1303.58*pow(r_ij,4)+196.754*pow(r_ij,5)-12.1613*pow(r_ij,6)+0.4972)*(cl[0]-si[i])/r_ij;
	fj=-(-3631.65+8776.68*r_ij-8700*pow(r_ij,2)+4527.67*pow(r_ij,3)-1303.58*pow(r_ij,4)+196.754*pow(r_ij,5)-12.1613*pow(r_ij,6)+0.4972)*(cl[1]-si[i+1])/r_ij;
	fk=-(-3631.65+8776.68*r_ij-8700*pow(r_ij,2)+4527.67*pow(r_ij,3)-1303.58*pow(r_ij,4)+196.754*pow(r_ij,5)-12.1613*pow(r_ij,6)+0.4972)*(cl[2]-si[i+2])/r_ij;

	f_cl[0]+=fi;
	f_cl[1]+=fj;
	f_cl[2]+=fk;

	f_si[i]-=fi;
	f_si[i+1]-=fj;
	f_si[i+2]-=fk;

      }
    }

	//繝ｻ繝ｻ繝ｻ貂ｩ蠎ｦ縺ｮ險育ｮ励・繝ｻ繝ｻ
	for(i=0,T=0;i<par*3;i+=3){
		T+=m_si*(pow(v_si[i],2)+pow(v_si[i+1],2)+pow(v_si[i+2],2));
	}
	
	T+=m_cl*(pow(v_cl[0],2)+pow(v_cl[1],2)+pow(v_cl[2],2));

	T/=3*(par+1)*kb;

	//繝ｻ繝ｻ繝ｻ貂ｩ蠎ｦ縺ｮ陬懈ｭ｣繝ｻ繝ｻ繝ｻ
	x=sqrt(T0/T);


    //繝ｻ繝ｻ繝ｻ騾溷ｺｦ縺ｮ譖ｴ譁ｰ・茨ｼ抵ｼ峨・繝ｻ繝ｻ
    //Si
    for(i=0;i<par*3;i+=3){
      v_si[i]=x*(v_si[i]+dt/2*f_si[i]/m_si*ev/angs);
      v_si[i+1]=x*(v_si[i+1]+dt/2*f_si[i+1]/m_si*ev/angs);
      v_si[i+2]=x*(v_si[i+2]+dt/2*f_si[i+2]/m_si*ev/angs);
    }
    
    //Cl
    v_cl[0]=x*(v_cl[0]+dt/2*f_cl[0]/m_cl*ev/angs);
    v_cl[1]=x*(v_cl[1]+dt/2*f_cl[1]/m_cl*ev/angs);
    v_cl[2]=x*(v_cl[2]+dt/2*f_cl[2]/m_cl*ev/angs);

    //繝ｻ繝ｻ繝ｻ繝輔ぃ繧､繝ｫ蜃ｺ蜉帙・繝ｻ繝ｻ
    if(count%50==0){
      FILE *outputfile;         // 蜃ｺ蜉帙せ繝医Μ繝ｼ繝
      char filename[100];
      sprintf(filename, "MD%06d.txt",count/50);  
      outputfile = fopen(filename, "w");  // 繝輔ぃ繧､繝ｫ繧呈嶌縺崎ｾｼ縺ｿ逕ｨ縺ｫ繧ｪ繝ｼ繝励Φ(髢九￥)
      if (outputfile == NULL) {          // 繧ｪ繝ｼ繝励Φ縺ｫ螟ｱ謨励＠縺溷ｴ蜷・
	printf("cannot open\n");// 繧ｨ繝ｩ繝ｼ繝｡繝・そ繝ｼ繧ｸ繧貞・縺励※
	exit(1);                         // 逡ｰ蟶ｸ邨ゆｺ・
      }
      
      for(i=0,j=0;i<par;i++,j+=3){
	fprintf(outputfile,"%d 2 %lf %lf %lf\n",j/3,si[j],si[j+1],si[j+2]);  // 繝輔ぃ繧､繝ｫ縺ｫ譖ｸ縺・
      }                    
      
      // for(i=0,j=par*3;i<par*17;i++,j+=3){
      //	fprintf(outputfile,"%d 3 %lf %lf %lf\n",j/3,si[j],si[j+1],si[j+2]);// 繝輔ぃ繧､繝ｫ縺ｫ譖ｸ縺・蜻ｨ譛溷｢・阜譚｡莉ｶ)
      //  }    

      fprintf(outputfile,"%d 0 %lf %lf %lf\n",par*9,cl[0],cl[1],cl[2]);// 繝輔ぃ繧､繝ｫ縺ｫ譖ｸ縺・cl)                        
      
      
      fclose(outputfile);          // 繝輔ぃ繧､繝ｫ繧偵け繝ｭ繝ｼ繧ｺ(髢峨§繧・
    }
}
    
  //繝ｻ繝ｻ繝ｻ邨ゆｺ・ｾ後・繝ｻ繝ｻ
  printf("fin\n");
  
    }




//繝ｻ繝ｻ繝ｻ繝ｻ繝ｻ繧ｵ繝悶・繝ｭ繧ｰ繝ｩ繝繝ｻ繝ｻ繝ｻ繝ｻ繝ｻ
//繝ｻ繝ｻ繝ｻ霍晞屬縺ｮ險育ｮ励・繝ｭ繧ｰ繝ｩ繝繝ｻ繝ｻ繝ｻ
double dis(double* atm_i,double* atm_j,int i,int j)
{
  double r;

  r=pow(atm_i[i]-atm_j[j],2)+pow(atm_i[i+1]-atm_j[j+1],2)+pow(atm_i[i+2]-atm_j[j+2],2);
  r=sqrt(r);
  return r;
}

//繝ｻ繝ｻ繝ｻ隗貞ｺｦ(cos)縺ｮ繝励Ο繧ｰ繝ｩ繝繝ｻ繝ｻ繝ｻ
double ang(double* atm_i,double* atm_j,double* atm_k,int i,int j,int k)
{
  double r_ij,r_ik,y,cos;
  double ij[3],ik[3];
  int x;

  for(x=0;x<3;x++){
    ij[x]=atm_i[i+x]-atm_j[j+x];
    ik[x]=atm_i[i+x]-atm_k[k+x];
  }

  y=ij[0]*ik[0]+ij[1]*ik[1]+ij[2]*ik[2];

  r_ij=dis(atm_i,atm_j,i,j);
  r_ik=dis(atm_i,atm_k,i,k);

  cos=y/(r_ij*r_ik);

  return cos;

}
