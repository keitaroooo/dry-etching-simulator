 #include <stdio.h>

 #include <stdlib.h>

 #include <string.h>

 #include <math.h>

 #include <time.h>

 #include "MT.h" //乱数を生成するヘッダーファイル

 #define pi 3.141592



 double dis(double[][3],double[][3],int,int);

 double ang(double[][3],double[][3],double[][3],int,int,int);

 const double avogadro=6.0221367e+23;

 const double kb=1.380648e-23; 






 
int main()

 {
  
   int i,j,k,l,m;
  
   int ix,iy,iz;
  
   int loop=2001;//メインループ繰り返し数
  
   int count=0;//メインループカウント

   double t,dt=1.0e-15;//時間
  
   double angs=1e-10;//1[Å]=1e-10[m]
  
   double ev=1.60218e-19;//1[eV]=1.60218e-19[J]

   double ev_angs=ev/angs;//速度の更新の際に単位変換として利用
  
   double T0=293.15,T,x;//温度[K](T0:設定温度)
  
   double co;//角度(cosθ)


   //原子パラメータ
  
   int par1_in_unit=8;//単位格子内の原子数
   
   int unit_num_x=3,unit_num_y=3,unit_num_z=2;//x,y,z方向の単位格子の数
   
   int par1=par1_in_unit*unit_num_x*unit_num_y*unit_num_z;//基本セル中の原子(Si)の数
   int par1_all=par1*18;//ミラーセルも含めた計算に用いる全原子の数(Si)
   
   int par2=5;//基本セル１つ当たりに近づける原子(Cl)の数
   int par2_all=par2*9;//ミラーセルも含めた計算に用いる全原子の数(Cl)
   
   int num_j[par1_all];//カットオフ距離内の原子の保管用(j)
   
   int num_k[par1_all];//カットオフ距離内の原子の保管用(k)
     
   double c1[par1_all][3],c2[par2_all][3];//座標［Å］
  
   double v1[par1_all][3],v2[par2_all][3],v;//速度［m/s］
  
   double f1[par1_all][3],f2[par2_all][3];//力［N］

   double m1=28.0855/avogadro*1e-3;//Siの原子質量[kg]
  
   double m2=35.453/avogadro*1e-3;//Clの原子質量[kg]

   double a=5.4301;//Siの格子定数[Å]

   double a_x=a*unit_num_x;//基本セルのx方向の1辺の長さ
   
   double a_y=a*unit_num_y;//基本セルのy方向の1辺の長さ
   
   double a_z=a*unit_num_z;//基本セルのz方向の1辺の長さ
 

   //tersoffパラメータ(Si)
   
   double lam1=2.4799;//[Åe-1]
  
   double lam2=1.7322;//[Åe-1]
  
   double beta=1.0999e-6;
  
   double n=7.8737e-1;
  
   double c=1.0039e+5;
  
   double d=1.6218e+1;
  
   double h=-5.9826e-1;

   double R=2.85;//Å
 
   double D=0.15;//Å
 
   double A=1.8308e+3;//[eV]

   double B=4.7118e+2;//[eV]

   double f[3];
   
   double dz[3],dz_i[3],dz_j[3],dz_k[3];
 
   double fr,fa,fc_ij,fc_ik,dfc_ij,dfc_ik,b,db,g,dg,z;
  
   double r_ij,r_ik,r_jk;
   

   //morseパラメータ
   
   //Cl_Si
   double E1=4.35057;
   
   double a1=1.43163;
   
   double ra1=2.1;

   double cut1=5.25;//SiとClのカットオフ距離[Å]

   
   //lenard jonesパラメータ
   
   //Cl_Cl
   double e1=316.0*kb;//ε/kb=316.0
   
   double o1=4.217;//σ(Å)

   double cut2=2.5*o1;//Cl同士のカットオフ距離[Å]
   

   //・・・・・初期条件・・・・・
  
   //・・・初期位置・・・
  

   //Si
  
   i=0;
  
   for(iz=0;iz<unit_num_z;iz++){
    
     for(iy=0;iy<unit_num_y;iy++){

       for(ix=0;ix<unit_num_x;ix++){
	
         c1[i][0]=0.+a*ix;
       
         c1[i][1]=0.+a*iy;

	 c1[i][2]=0.+a*iz;
  
	 i++;

	
         c1[i][0]=0.+a*ix;

	 c1[i][1]=1./2.*a+a*iy;

	 c1[i][2]=1./2.*a+a*iz;
	
         i++;

	
         c1[i][0]=1./2.*a+a*ix;
	
         c1[i][1]=0.+a*iy;
	
         c1[i][2]=1./2.*a+a*iz;
	
         i++;

	
         c1[i][0]=1./2.*a+a*ix;
	
         c1[i][1]=1./2.*a+a*iy;
	
         c1[i][2]=0.+a*iz;
         
         i++;

	
         c1[i][0]=1./4.*a+a*ix;
	
         c1[i][1]=1./4.*a+a*iy;
	
         c1[i][2]=1./4.*a+a*iz;
	
         i++;

	
         c1[i][0]=1./4.*a+a*ix;
	
         c1[i][1]=3./4.*a+a*iy;
	 
         c1[i][2]=3./4.*a+a*iz;
	
         i++;


         c1[i][0]=3./4.*a+a*ix;
	
         c1[i][1]=1./4.*a+a*iy;
	
         c1[i][2]=3./4.*a+a*iz;
	
         i++;

		
         c1[i][0]=3./4.*a+a*ix;
	
         c1[i][1]=3./4.*a+a*iy;
	
         c1[i][2]=1./4.*a+a*iz;
	
         i++;
      
       }  
     } 
   }

   //Cl

   init_genrand((unsigned)time(NULL));//現在時刻で乱数の初期化
   for(i=0;i<par2;i++){L4:

     c2[i][0]=unit_num_x*a*genrand_res53();//genrand_res53():0~1の乱数

     c2[i][1]=unit_num_y*a*genrand_res53();//0~2aの範囲の値を与える
  
     c2[i][2]=(unit_num_z-1)*a+3./4.*a+cut1;//基本格子(上面のSi原子)からカットオフ距離だけ離す
     
     if (c2[i][0]==c2[i-1][0]&&c2[i][1]==c2[i-1][1])goto L4;
  
   }

   
   //ミラーセルにコピー
  
   j=1;
   
     for(iy=-1;iy<=1;iy++){

       for(ix=-1;ix<=1;ix++){if(ix==0&&iy==0)ix=1;

	 //Siの上段(原子番号=par1+1,2・・・)
	 k=par1*j;
         for(i=0;i<par1;i++){
           c1[k+i][0]=c1[i][0]+a_x*ix;
           c1[k+i][1]=c1[i][1]+a_y*iy;
           c1[k+i][2]=c1[i][2];
         }

	 //Cl(原子番号=par2+1,2・・・)
	 k=par2*j;
	 for(i=0;i<par2;i++){
           c2[k+i][0]=c2[i][0]+a_x*ix;
           c2[k+i][1]=c2[i][1]+a_y*iy;
           c2[k+i][2]=c2[i][2];
	 }
	 j++;
       }
     }

    
     for(iy=-1;iy<=1;iy++){

       for(ix=-1;ix<=1;ix++){

	 k=par1*j;

         for(i=0;i<par1;i++){
           c1[k+i][0]=c1[i][0]+a_x*ix;
           c1[k+i][1]=c1[i][1]+a_y*iy;
           c1[k+i][2]=c1[i][2]-a_z;
	 }
	 j++;
       }
     }


   //・・・初期速度・・・
 
   //Si

   for(i=0;i<par1;i++){
     for(j=0;j<3;j++){
  
     v1[i][j]=0;

     }
   }

 
   //Cl

   double r1,r2;
   r1=sqrt(2*kb*T0/m2);//熱確速度
   r2=2*pi;//2π

   for(i=0;i<par2;i++){
     v2[i][0]=r1*sqrt(-log(genrand_res53()))*cos(r2*genrand_res53());
     v2[i][1]=r1*sqrt(-log(genrand_res53()))*cos(r2*genrand_res53());
     v2[i][2]=-r1*sqrt(-log(genrand_res53()))*fabs(cos(r2*genrand_res53()));//基盤原子の方向となるように負にする
   }


   //・・・力の初期条件・・・

   for(j=0;j<3;j++){

     //Si
     for(i=0;i<par1;i++){
       f1[i][j]=0;
     }

     //Cl
     for(i=0;i<par2;i++){
     f2[i][j]=0;

     }
   }
   
   
   //・・・・・メインループ・・・・・

   for(t=0;t<dt*loop;t+=dt,count++){

     //・・・位置の更新・・・
     for(j=0;j<3;j++){

     //Si
       for(i=0;i<par1;i++){
	 v1[i][j]=v1[i][j]+dt/2*f1[i][j]/m1*ev_angs;//速度の更新でも使用するので計算回数を減らすために保管しておく(ev_angsでm/sに変換しているf:ev/Å→J/m)
         c1[i][j]=c1[i][j]+dt*v1[i][j]/angs;//angsでmをÅに変換している
       }

     //Cl
       for(i=0;i<par2;i++){
	 v2[i][j]=v2[i][j]+dt/2*f2[i][j]/m2*ev_angs;
         c2[i][j]=c2[i][j]+dt*v2[i][j]/angs;
       }
     }


     //・・・周期条件・・・

     //基本セルの外に出た粒子の修正

     //Si
     for(i=0;i<par1;i++){
       if(c1[i][0]<0){
	
       while(c1[i][0]<0){
           c1[i][0]+=a_x;
         }
       }
      
       else if(c1[i][0]>a_x){

         while(c1[i][0]>a_x){

           c1[i][0]-=a_x;
	 }
       }
      
        if(c1[i][1]<0){
	
          while(c1[i][1]<0){

           c1[i][1]+=a_y;
	  }
	}
     
       else if(c1[i][1]>a_y){

         while(c1[i][1]>a_y){

           c1[i][1]-=a_y;
	 }
       }
     }

     //Cl
     for(i=0;i<par2;i++){
       if(c2[i][0]<0){
	
       while(c2[i][0]<0){
           c2[i][0]+=a_x;
         }
       }
      
       else if(c2[i][0]>a_x){

         while(c2[i][0]>a_x){

           c2[i][0]-=a_x;
	 }
       }
       
        if(c2[i][1]<0){
	
          while(c2[i][1]<0){

           c2[i][1]+=a_y;
	  }
	}
	
       else if(c2[i][1]>a_y){

         while(c2[i][1]>a_y){

           c2[i][1]-=a_y;
	 }
       }
     }

     //ミラーセルにコピー

     j=1;
    
       for(iy=-1;iy<=1;iy++){

         for(ix=-1;ix<=1;ix++){if(ix==0&&iy==0)ix=1;

	   //Siの上段(原子番号=par1+1,2・・・)
	   k=par1*j;
           for(i=0;i<par1;i++){
             c1[k+i][0]=c1[i][0]+a_x*ix;
             c1[k+i][1]=c1[i][1]+a_y*iy;
             c1[k+i][2]=c1[i][2];
	   }

	   //Cl(原子番号=par2+1,2・・・)
	   k=par2*j;
	   for(i=0;i<par2;i++){
             c2[k+i][0]=c2[i][0]+a_x*ix;
             c2[k+i][1]=c2[i][1]+a_y*iy;
             c2[k+i][2]=c2[i][2];
	   }
	   j++;
	 }
       }


     //・・・力の初期化・・・

     for(j=0;j<3;j++){

       //Si
       for(i=0;i<par1;i++){
	 f1[i][j]=0;
       }

       //Cl
       for(i=0;i<par2;i++){
         f2[i][j]=0;

       }
     }

     
     //・・・力の計算・・・


     //Si

     for(i=0;i<par1_all;i++){

     //カットオフ距離内の原子の保管(j)

   
     for(j=0,l=0;j<par1_all;j++){       
       if(j!=i){

	 if(dis(c1,c1,i,j)<R+D){

	     num_j[l]=j;

	     l++;
	 
       } 
} 
}

       //カットオフ距離内の原子の保管(k)
	
       for(j=0;j<l;j++){ 
         num_k[j]=num_j[j];  
       }



       for(j=0;j<l;j++){
         //i原子とj原子間におけるfrとfaの計算
         r_ij=dis(c1,c1,i,num_j[j]);

         fr=A*exp(-lam1*r_ij);

         fa=-B*exp(-lam2*r_ij);

           //i原子とj原子間におけるfcとd(fc)/drの計算
           if(R-D<r_ij){

             fc_ij=1./2.*(1-sin(pi/2*(r_ij-R)/D));

             dfc_ij=-pi/(4*D)*cos(pi/2*(r_ij-R)/D);

           }

           else{
	
             fc_ij=1;
	 
             dfc_ij=0;
	
           }


         //ζ(=z)の計算(i原子とk原子間におけるfcの計算)
         z=0;
         for(k=0;k<l;k++){	 
           if(num_k[k]!=num_j[j]){

             r_ik=dis(c1,c1,i,num_k[k]);
	 
             if(R-D<r_ik){
	    
               fc_ik=1./2.*(1-sin(pi/2*(r_ik-R)/D));

       
	     }

	     else{

	       fc_ik=1;
	     }

	   g=1+c*c/d/d-c*c/(d*d+pow(h-ang(c1,c1,c1,i,num_j[j],num_k[k]),2));

	   z+=fc_ik*g;

	 } 
}


	 
         //bの計算
         b=1+pow(beta,n)*pow(z,n);


	 if(z!=0)db=-1./2.*pow(beta,n)*pow(z,n-1)*pow(b,-1/(2*n)-1);

	 else db=0;


	 b=pow(b,-1/(2*n));


         //Fi(力)の計算(第1項)
         for(k=0;k<3;k++){
	   f[k]=(-dfc_ij*fr-fc_ij*(-lam1*fr))*(c1[i][k]-c1[num_j[j]][k])/r_ij;//-lam1*fr=dfr

	   f[k]+=-b*(dfc_ij*fa+fc_ij*(-lam2*fa))*(c1[i][k]-c1[num_j[j]][k])/r_ij;//-lam2*fa=dfa
	   f[k]*=1./2.;//第1項に×1/2していることになっている

           f1[i][k]+=f[k];
	   f1[num_j[j]][k]-=f[k];


	   }
 
         for(k=0;k<l;k++){
           if(num_k[k]!=num_j[j]){

             r_ik=dis(c1,c1,i,num_k[k]);
	 
             if(R-D<r_ik){
	    
               fc_ik=1./2.*(1-sin(pi/2*(r_ik-R)/D));

               dfc_ik=-pi/(4*D)*cos(pi/2*(r_ik-R)/D);
       
	     }

	     else{

	       fc_ik=1;

	       dfc_ik=0;

	     }

             //gとdgの計算
             co=ang(c1,c1,c1,i,num_j[j],num_k[k]);

	     g=1+c*c/d/d-c*c/(d*d+pow(h-co,2));

	     dg=-2*c*c*(h-co)/pow(d*d+pow(h-co,2),2);


	     
             //∂ζ/∂riの計算(第1項)
             for(m=0;m<3;m++){
               dz[m]=g*dfc_ik*(c1[i][m]-c1[num_k[k]][m])/r_ik;

               
 	       //∂ζ/∂riの計算
	       dz_i[m]=dz[m]+fc_ik*dg*((co/r_ij-1/r_ik)*(c1[num_j[j]][m]-c1[i][m])/r_ij+(co/r_ik-1/r_ij)*(c1[num_k[k]][m]-c1[i][m])/r_ik);

	       
               //Fi(力)の計算           
               f1[i][m]-=1./2.*(fc_ij*fa*db*dz_i[m]);//第2項に×1/2していることになっている


	       //∂ζ/∂rjの計算
 
	       dz_j[m]=fc_ik*dg/r_ij*(-co*(c1[num_j[j]][m]-c1[i][m])/r_ij+(c1[num_k[k]][m]-c1[i][m])/r_ik);

	    	    
               //Fj(力)の計算
               f1[num_j[j]][m]-=1./2.*(fc_ij*fa*db*dz_j[m]);//第2項に×1/2していることになっている

	    
      
	   
               //∂ζ/∂rkの計算(-∂ζ/∂rk)
               dz_k[m]=-dz[m]+fc_ik*dg/r_ik*(-co*(c1[num_k[k]][m]-c1[i][m])/r_ik+(c1[num_j[j]][m]-c1[i][m])/r_ij);

	       	    
               //Fk(力)の計算
               f1[num_k[k]][m]-=1./2.*(fc_ij*fa*db*dz_k[m]);//第2項に×1/2していることになっている


	     }//for(m=0:m<3;m++){
           }//if(num_k[k]!=num_j[j]){

         }//for(k=0;k<l;k++){
       }//for(j=0;j<l;l++){
     }//for(i=0;i<par1_all;i++){


     //Cl

     for(i=0;i<par2_all;i++){

       //近づける原子と基盤原子間に働く力の計算
       for(j=0;j<par1_all;j++){

         r_ij=dis(c2,c1,i,j);

         if(r_ij<cut1){
           for(k=0;k<3;k++){


             f[k]=-(2*a1*E1*(exp(-a1*(r_ij-ra1))-exp(-2*a1*(r_ij-ra1))))*(c2[i][k]-c1[j][k])/r_ij;//-∂U/∂r*(c2-c1)/rでc1→c2方向が正

             f2[i][k]+=f[k];
             f1[j][k]-=f[k];
           }
         
}


       }

       //近づける原子同士の間に働く力の計算
       for(j=i+1;j<par2_all;j++){

         r_ij=dis(c2,c2,i,j);

         if(r_ij<cut2){
           for(k=0;k<3;k++){

	     f[k]=-4*e1*(-12*pow(o1,12)/pow(r_ij,13)+6*pow(o1,6)/pow(r_ij,7))*(c2[i][k]-c2[j][k])/r_ij/ev;//-∂U/∂r*(c2-c1)/rでj→i方向が正,/ev:J→ev

	     f2[i][k]+=f[k];//(斥力が働いていてi原子のほうが正方向にあるときf3[i][k]>0)
	     f2[j][k]-=f[k];
	   }
	 }
       }
     }


     //・・・速度の更新・・・
    
     for(i=0;i<3;i++){
       //Si
    
       for(j=0;j<par1;j++){
      
         v1[j][i]=v1[j][i]+dt/2*f1[j][i]/m1*ev_angs;
       }

       //Cl
       for(j=0;j<par2;j++){
      
         v2[j][i]=v2[j][i]+dt/2*f2[j][i]/m2*ev_angs;
       }
     }
     
     //・・・温度の計算・・・

     for(i=0,T=0;i<par1;i++){

       T+=m1*(pow(v1[i][0],2)+pow(v1[i][1],2)+pow(v1[i][2],2));
	
     }
     for(i=0;i<par2;i++){

       T+=m2*(pow(v2[i][0],2)+pow(v2[i][1],2)+pow(v2[i][2],2));

     }
     
     T/=3*(par1+par2)*kb;

     //・・・温度の補正・・・
	
     x=sqrt(T0/T);

     //・・・速度の更新(温度による補正)・・・
    
     for(i=0;i<3;i++){
       //Si
    
       for(j=0;j<par1;j++){
      
         v1[j][i]*=x;
       }

       //Cl
       for(j=0;j<par2;j++){
      
         v2[j][i]*=x;
       }
     }
     
       
     //・・・ファイル出力・・・
    
     if(count%50==0){
      
       FILE *outputfile;//出力ストリーム(ファイルを定義)
      
       char filename[100];
      
       sprintf(filename, "MD%06d.txt",count/50);//MD%06d.txtという文字列を配列filenameに挿入 
      
       outputfile = fopen(filename, "w");//ファイルを書き込み用にオープン(MD%06d.txtというファイルを開く)
      
       if (outputfile == NULL) {//オープンに失敗した場合
	
         printf("cannot open\n");// エラーメッセージを出して

       
	 exit(1);// 異常終了
      
       }

      
      
       for(i=0;i<par1;i++){
	 fprintf(outputfile,"%d 2 %lf %lf %lf\n",i,c1[i][0],c1[i][1],c1[i][2]);//MD%06d.txtというファイルに書く
       }                          
      
     
       // for(i=par1;i<par1*18;i++){
      
       //fprintf(outputfile,"%d 3 %lf %lf %lf\n",i,c1[i][0],c1[i][1],c1[i][2]);//ファイルに書く(周期境界条件)今回は必要ない
      
       //}
       
       for(i=0;i<par2;i++){
	 fprintf(outputfile,"%d 0 %lf %lf %lf\n",i,c2[i][0],c2[i][1],c2[i][2]);//MD%06d.txtというファイルに書く(cl),[通し番号][粒子種][x座標][y座標][z座標]
       }

       
       //   for(i=par2;i<par2*9;i++){
       //fprintf(outputfile,"%d 0 %lf %lf %lf\n",i,c2[i][0],c2[i][1],c2[i][2]);//ファイルに書く(周期境界条件)今回は必要ない 
	 //}                       
      
      
     
       fclose(outputfile);// ファイルをクローズ(閉じる)
     
   } } 
  
   //・・・終了後・・・
  
   printf("fin\n");
  
 }






 //・・・・・サブプログラム・・・・・


 //・・・距離の計算プログラム・・・

 double dis(double atm_i[][3],double atm_j[][3],int i,int j)

 {
   int k;
   double r=0;


   for(k=0;k<3;k++){
     r+=pow(atm_i[i][k]-atm_j[j][k],2);

   }
   r=sqrt(r);
  
   return r;

 }



 //・・・角度(cos)のプログラム・・・

 double ang(double atm_i[][3],double atm_j[][3],double atm_k[][3],int i,int j,int k)

 {
   double r_ij,r_ik,y,cos;

   double ij[3],ik[3];

   int l;


   for(l=0;l<3;l++){

     ij[l]=atm_i[i][l]-atm_j[j][l];

     ik[l]=atm_i[i][l]-atm_k[k][l];

     y+=ij[l]*ik[l];
     }




   r_ij=dis(atm_i,atm_j,i,j);

   r_ik=dis(atm_i,atm_k,i,k);


   cos=y/(r_ij*r_ik);


   return cos;


 }


