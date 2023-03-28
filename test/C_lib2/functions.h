

double P_interpol_fast_doublearray(double k0, double **P, int index, int N, char *spacing,int interpolation_order, int Ninterpol,double w0,double w1,double w2 );

int determine_N_doublearray(double **Theory,double kinput,int Nlin, char *spacing);

///
double determine_w1_doublearray(double **Theory,double kinput,int N1,char *spacing);

double P_interpol_w1_singlearray(double *Theory, int N1,double w1, char *spacing);

double P_interpol_w1_doublearray(double **Theory,int index, int N1,double w1,char *spacing);

//


double determine_w0_2ndorder_doublearray(double **k,double kinput,int N1,char *spacing);

double determine_w1_2ndorder_doublearray(double **k,double kinput,int N1, char *spacing);

double determine_w2_2ndorder_doublearray(double **k,double kinput,int N1, char *spacing);

double P_interpol_w012_singlearray(double *Theory, int N1,double w0, double w1, double w2, char *spacing);

double P_interpol_w012_doublearray(double **Theory,int index, int N1,double w0, double w1, double w2,char *spacing);

void freeTokens(double** tokens, int N);
