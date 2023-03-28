#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <string.h>

#define MIN(a,b) (((a)<(b)) ? (a) : (b) ) // the usual "min" function
#define MAX(a,b) (((a)<(b)) ? (b) : (a) ) // the usual "max" function



    double determine_w1_doublearray(double **Theory,double kinput,int N1,char *spacing)
   {
       double w1;
if(strcmp(spacing,"linear") == 0 ){w1=(Theory[N1+1][0]-kinput)/(Theory[N1+1][0]-Theory[N1][0]);}
if(strcmp(spacing,"log") == 0 ){w1=(log(Theory[N1+1][0])-log(kinput))/(log(Theory[N1+1][0])-log(Theory[N1][0]));}
if(strcmp(spacing,"log10") == 0 ){w1=(log10(Theory[N1+1][0])-log10(kinput))/(log10(Theory[N1+1][0])-log10(Theory[N1][0]));}

       return w1;
   }


double P_interpol_w1_doublearray(double **Theory,int index, int N1,double w1,char *spacing)
{
double a;
if(strcmp(spacing,"linear") == 0 ){a=w1*Theory[N1][index]+(1.-w1)*Theory[N1+1][index];}

if(strcmp(spacing,"log") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]>0  ){a=w1*log(Theory[N1][index])+(1.-w1)*log(Theory[N1+1][index]);a=exp(a);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]<0  ){a=w1*log(-Theory[N1][index])+(1.-w1)*log(-Theory[N1+1][index]);a=-exp(a);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]>0  ){a=w1*(Theory[N1][index])+(1.-w1)*(Theory[N1+1][index]);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]<0  ){a=w1*(Theory[N1][index])+(1.-w1)*(Theory[N1+1][index]);}

if(strcmp(spacing,"log10") == 0 ){a=w1*log10(Theory[N1][index])+(1.-w1)*log10(Theory[N1+1][index]);a=pow(10,a);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]<0  ){a=w1*log10(-Theory[N1][index])+(1.-w1)*log10(-Theory[N1+1][index]);a=-pow(10,a);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]>0  ){a=w1*(Theory[N1][index])+(1.-w1)*(Theory[N1+1][index]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]<0  ){a=w1*(Theory[N1][index])+(1.-w1)*(Theory[N1+1][index]);}

return a;
}

double P_interpol_w1_singlearray(double *Theory, int N1,double w1, char *spacing)
{
double a;
if(strcmp(spacing,"linear") == 0 ){a=w1*Theory[N1]+(1.-w1)*Theory[N1+1];}

if(strcmp(spacing,"log") == 0 && Theory[N1]>0 && Theory[N1+1]>0){a=w1*log(Theory[N1])+(1.-w1)*log(Theory[N1+1]);a=exp(a);}
if(strcmp(spacing,"log") == 0 && Theory[N1]<0 && Theory[N1+1]<0 ){a=w1*log(-Theory[N1])+(1.-w1)*log(-Theory[N1+1]);a=-exp(a);}
if(strcmp(spacing,"log") == 0 && Theory[N1]>0 && Theory[N1+1]<0 ){a=w1*(Theory[N1])+(1.-w1)*(Theory[N1+1]);}
if(strcmp(spacing,"log") == 0 && Theory[N1]<0 && Theory[N1+1]>0 ){a=w1*(Theory[N1])+(1.-w1)*(Theory[N1+1]);}

if(strcmp(spacing,"log10") == 0 && Theory[N1]>0 && Theory[N1+1]>0){a=w1*log10(Theory[N1])+(1.-w1)*log10(Theory[N1+1]);a=pow(10,a);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]<0 && Theory[N1+1]<0){a=w1*log10(-Theory[N1])+(1.-w1)*log10(-Theory[N1+1]);a=-pow(10,a);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]>0 && Theory[N1+1]<0){a=w1*(Theory[N1])+(1.-w1)*(Theory[N1+1]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]<0 && Theory[N1+1]>0){a=w1*(Theory[N1])+(1.-w1)*(Theory[N1+1]);}

return a;
}


double P_interpol_w012_singlearray(double *Theory, int N1,double w0, double w1, double w2, char *spacing)
{

double a;
if(strcmp(spacing,"linear") == 0){a=w0*Theory[N1]+w1*Theory[N1+1]+w2*Theory[N1+2];}

if(strcmp(spacing,"log") == 0 && Theory[N1]>0 && Theory[N1+1]>0 && Theory[N1+2]>0){a=w0*log(Theory[N1])+w1*log(Theory[N1+1])+w2*log(Theory[N1+2]);a=exp(a);}
if(strcmp(spacing,"log") == 0 && Theory[N1]<0 && Theory[N1+1]<0 && Theory[N1+2]<0){a=w0*log(-Theory[N1])+w1*log(-Theory[N1+1])+w2*log(-Theory[N1+2]);a=-exp(a);}

if(strcmp(spacing,"log") == 0 && Theory[N1]<0 && Theory[N1+1]<0 && Theory[N1+2]>0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log") == 0 && Theory[N1]<0 && Theory[N1+1]>0 && Theory[N1+2]<0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log") == 0 && Theory[N1]>0 && Theory[N1+1]<0 && Theory[N1+2]<0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}

if(strcmp(spacing,"log") == 0 && Theory[N1]>0 && Theory[N1+1]>0 && Theory[N1+2]<0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log") == 0 && Theory[N1]>0 && Theory[N1+1]<0 && Theory[N1+2]>0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log") == 0 && Theory[N1]<0 && Theory[N1+1]>0 && Theory[N1+2]>0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}


if(strcmp(spacing,"log10") == 0 && Theory[N1]>0 && Theory[N1+1]>0 && Theory[N1+2]>0){a=w0*log10(Theory[N1])+w1*log10(Theory[N1+1])+w2*log10(Theory[N1+2]);a=pow(10,a);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]<0 && Theory[N1+1]<0 && Theory[N1+2]<0){a=w0*log10(-Theory[N1])+w1*log10(-Theory[N1+1])+w2*log10(-Theory[N1+2]);a=-pow(10,a);}

if(strcmp(spacing,"log10") == 0 && Theory[N1]<0 && Theory[N1+1]<0 && Theory[N1+2]>0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]<0 && Theory[N1+1]>0 && Theory[N1+2]<0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]>0 && Theory[N1+1]<0 && Theory[N1+2]<0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}

if(strcmp(spacing,"log10") == 0 && Theory[N1]>0 && Theory[N1+1]>0 && Theory[N1+2]<0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]>0 && Theory[N1+1]<0 && Theory[N1+2]>0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1]<0 && Theory[N1+1]>0 && Theory[N1+2]>0){a=w0*(Theory[N1])+w1*(Theory[N1+1])+w2*(Theory[N1+2]);}


return a;
}

double P_interpol_w012_doublearray(double **Theory,int index, int N1,double w0, double w1, double w2,char *spacing)
{
double a;
if(strcmp(spacing,"linear") == 0){a=w0*Theory[N1][index]+w1*Theory[N1+1][index]+w2*Theory[N1+2][index];}

//if(strcmp(spacing,"log") == 0){a=w0*log(Theory[N1][index])+w1*log(Theory[N1+1][index])+w2*log(Theory[N1+2][index]);a=exp(a);}
//if(strcmp(spacing,"log10") == 0){a=w0*log10(Theory[N1][index])+w1*log10(Theory[N1+1][index])+w2*log10(Theory[N1+2][index]);a=pow(10,a);}

if(strcmp(spacing,"log") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]>0 && Theory[N1+2][index]>0){a=w0*log(Theory[N1][index])+w1*log(Theory[N1+1][index])+w2*log(Theory[N1+2][index]);a=exp(a);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]<0 && Theory[N1+2][index]<0){a=w0*log(-Theory[N1][index])+w1*log(-Theory[N1+1][index])+w2*log(-Theory[N1+2][index]);a=-exp(a);}

if(strcmp(spacing,"log") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]<0 && Theory[N1+2][index]>0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]>0 && Theory[N1+2][index]<0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]<0 && Theory[N1+2][index]<0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}

if(strcmp(spacing,"log") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]>0 && Theory[N1+2][index]<0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]<0 && Theory[N1+2][index]>0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]>0 && Theory[N1+2][index]>0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}


if(strcmp(spacing,"log10") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]>0 && Theory[N1+2][index]>0){a=w0*log10(Theory[N1][index])+w1*log10(Theory[N1+1][index])+w2*log10(Theory[N1+2][index]);a=pow(10,a);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]<0 && Theory[N1+2][index]<0){a=w0*log10(-Theory[N1][index])+w1*log10(-Theory[N1+1][index])+w2*log10(-Theory[N1+2][index]);a=-pow(10,a);}

if(strcmp(spacing,"log10") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]<0 && Theory[N1+2][index]>0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]>0 && Theory[N1+2][index]<0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]<0 && Theory[N1+2][index]<0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}

if(strcmp(spacing,"log10") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]>0 && Theory[N1+2][index]<0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]>0 && Theory[N1+1][index]<0 && Theory[N1+2][index]>0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}
if(strcmp(spacing,"log10") == 0 && Theory[N1][index]<0 && Theory[N1+1][index]>0 && Theory[N1+2][index]>0){a=w0*(Theory[N1][index])+w1*(Theory[N1+1][index])+w2*(Theory[N1+2][index]);}


return a;
}
         

  double determine_w0_2ndorder_doublearray(double **k,double kinput,int N1,char *spacing)
   {
       double w1;
       if(strcmp(spacing,"linear") == 0){w1=(kinput-k[N1+1][0])*(kinput-k[N1+2][0])/((k[N1][0]-k[N1+1][0] )*(k[N1][0]-k[N1+2][0]));}
       if(strcmp(spacing,"log") == 0){w1=(log(kinput)-log(k[N1+1][0]))*(log(kinput)-log(k[N1+2][0]))/((log(k[N1][0])-log(k[N1+1][0]))*(log(k[N1][0])-log(k[N1+2][0])));}
       if(strcmp(spacing,"log10") == 0){w1=(log10(kinput)-log10(k[N1+1][0]))*(log10(kinput)-log10(k[N1+2][0]))/((log10(k[N1][0])-log10(k[N1+1][0]))*(log10(k[N1][0])-log10(k[N1+2][0])));}


       return w1;
   }

    double determine_w1_2ndorder_doublearray(double **k,double kinput,int N1, char *spacing)
   {
       double w1;
       if(strcmp(spacing,"linear") == 0){w1=(kinput-k[N1][0])*(kinput-k[N1+2][0])/((k[N1+1][0]-k[N1][0] )*(k[N1+1][0]-k[N1+2][0]));}
       if(strcmp(spacing,"log") == 0){w1=(log(kinput)-log(k[N1][0]))*(log(kinput)-log(k[N1+2][0]))/((log(k[N1+1][0])-log(k[N1][0]))*(log(k[N1+1][0])-log(k[N1+2][0])));}
       if(strcmp(spacing,"log10") == 0){w1=(log10(kinput)-log10(k[N1][0]))*(log10(kinput)-log10(k[N1+2][0]))/((log10(k[N1+1][0])-log10(k[N1][0]))*(log10(k[N1+1][0])-log10(k[N1+2][0])));}


       return w1;
   }

  double determine_w2_2ndorder_doublearray(double **k,double kinput,int N1, char *spacing)
   {
       double w1;
       if(strcmp(spacing,"linear") == 0){w1=(kinput-k[N1][0])*(kinput-k[N1+1][0])/((k[N1+2][0]-k[N1][0])*(k[N1+2][0]-k[N1+1][0]));}
       if(strcmp(spacing,"log") == 0){w1=(log(kinput)-log(k[N1][0]))*(log(kinput)-log(k[N1+1][0]))/((log(k[N1+2][0])-log(k[N1][0]))*(log(k[N1+2][0])-log(k[N1+1][0])));}
       if(strcmp(spacing,"log10") == 0){w1=(log10(kinput)-log10(k[N1][0]))*(log10(kinput)-log10(k[N1+1][0]))/((log10(k[N1+2][0])-log10(k[N1][0]))*(log10(k[N1+2][0])-log10(k[N1+1][0])));}


       return w1;
   }


double P_interpol_fast_doublearray(double k0, double **P, int index, int N, char *spacing,int interpolation_order, int Ninterpol,double w0,double w1,double w2 )
{
int shiftN;
double Pout;
interpolation_order=1;//1 or 2
shiftN=interpolation_order;
if(Ninterpol>=N-shiftN || Ninterpol<0 || k0<=0){Pout=0;}
else{
if(interpolation_order==1){Pout=P_interpol_w1_doublearray(P,index,Ninterpol,w1,spacing);}
if(interpolation_order==2){Pout=P_interpol_w012_doublearray(P,index,Ninterpol,w0,w1,w2,spacing);}
}

return Pout;
}



 int determine_N_doublearray(double **Theory,double kinput,int Nlin, char *spacing)
{
       int N1;
       int i;

  double check1;
       double check1log;
       double check1log10;
     
       check1=(Theory[Nlin-1][0]-Theory[0][0])/(Nlin*1.-1.);
       check1log=(log(Theory[Nlin-1][0])-log(Theory[0][0]))/(Nlin*1.-1.);
       check1log10=(log10(Theory[Nlin-1][0])-log10(Theory[0][0]))/(Nlin*1.-1.);

if(kinput<Theory[0][0] || kinput>Theory[Nlin-1][0] || kinput<=0){N1=-1;}
else{
       if(strcmp(spacing,"linear") == 0){N1=(int)((kinput-Theory[0][0])/check1);}
       if(strcmp(spacing,"log") == 0){N1=(int)((log(kinput)-log(Theory[0][0]))/check1log);}
       if(strcmp(spacing,"log10") == 0){N1=(int)((log10(kinput)-log10(Theory[0][0]))/check1log10);}
       if(strcmp(spacing,"irregular") == 0){
       i=-1;
            do{

               i++;
               }while(kinput>Theory[i][0]);
               N1=i-1;

       }

}
       return N1;

}

void freeTokens(double **tokens, int N)
{
int i;
for(i=0;i<N;++i)
{    
     free(tokens[i]);
}
free(tokens);

}
