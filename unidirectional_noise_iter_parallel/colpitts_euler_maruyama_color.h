//#include <random>
//#include <stdio.h>
//#include <string.h>
//#include <math.h>
//#include <ctime>    // For time()
//#include <cstdlib>  // For srand() and rand()
//#include <iostream>
//#include <fstream>
#include "nr3.h"
#include "ran.h"
#include "deviates.h"
#include "gamma.h"

//----------------------------Coupled-Colpitts Oscillator Equations---------------------

void rhs_ode(int Ncells, double Vc[], double Ve[], double IL[],double par[12]){

double Vcc = par[0];
double Rl  = par[1];
double L   = par[2];
double C1  = par[3];
double C2  = par[4];
double Ree = par[5];
double Vee = par[6];

double beta = par[7];
double Ron  = par[8];
double Vth  = par[9];
double Vb   = par[10];
double lambda = par[11];

double Vcnew[Ncells];
double Venew[Ncells];
double ILnew[Ncells];


double Ib;
double Ic;
double term1;

int jp = 0;

for (int jj = 1; jj < Ncells+1; jj++){

        jp = jj % Ncells;
        if ((Vb - Ve[jj-1]) <= Vth){
            Ib = 0;
        }
        else if ((Vb - Ve[jj-1]) > Vth){
            Ib = (Vb-Ve[jj-1]-Vth)/Ron;
        }
        Ic = beta*Ib;

        term1 = ((Vee + Vb - Ve[jj-1])/Ree + IL[jj-1] + Ib)/C2;

        Vcnew[jj-1] = (IL[jj-1]-Ic)/C1 + term1 + 1/C1*(lambda*(Vc[jp]-Vc[jj-1])/Ree);
        Venew[jj-1] = term1;
        ILnew[jj-1] = (Vcc-Vc[jj-1]+Vb-IL[jj-1]*Rl)/L;

}

for (int i1 = 0; i1 < Ncells; i1++){
    Vc[i1] = Vcnew[i1];
    Ve[i1] = Venew[i1];
    IL[i1] = ILnew[i1];
}
}


int intRand(const int & min, const int & max) {
    static thread_local std::mt19937 generator;
    std::uniform_int_distribution<int> distribution(min,max);
    return distribution(generator);
}


//---------------------------Euler-Maruyama Integration with Noise----------------------------------


void euler_maruyama(int N, int Ncells, int trans, double solution_data[], double solution_noise[], double eta[], double Vc[], double Ve[], double IL[], double par[12], double h, int seed, double np[]){

double Vcsol[Ncells];  //preallocating for noisy integration
for (int i1=0;i1<Ncells;i1++){
    Vcsol[i1] = Vc[i1];
}
double Vesol[Ncells];
for (int i2=0;i2<Ncells;i2++){
    Vesol[i2] = Ve[i2];
}
double ILsol[Ncells];
for (int i3=0;i3<Ncells;i3++){
    ILsol[i3] = IL[i3];
}



double Vcsol2[Ncells];  //preallocating for noiseless integration
for (int i1=0;i1<Ncells;i1++){
    Vcsol2[i1] = Vc[i1];
}
double Vesol2[Ncells];
for (int i2=0;i2<Ncells;i2++){
    Vesol2[i2] = Ve[i2];
}
double ILsol2[Ncells];
for (int i3=0;i3<Ncells;i3++){
    ILsol2[i3] = IL[i3];
}

double* randn;
randn = new double[Ncells];

double mu = 0, sigma = 1;
//int seed = 2;//37
std::random_device rd;                               //only used when randn[rr] = dist1(generator)
//std::mt19937 generator(time(0));
std::mt19937 generator(rd());
std::normal_distribution<double> dist1(mu,sigma);

Normaldev_BM GaussRNG(mu,sigma,seed);

//double D = 5.0e05;            double tau_c = 1.0e-03;

double D = np[0];
double tau_c = np[1];
//double sqrtdht = sqrt(h);
double htauc = h/tau_c,     sqrtdht = sqrt(2*D*h)/tau_c;



//int trans = 20; //how many times to integrate time array for transient (transient = 3*N)

//-------------------------Integration with Noise------------------
rhs_ode(Ncells,Vc,Ve,IL,par);

for (int iter = 0; iter < trans*N; iter++){
        for (int rr = 0; rr < Ncells; rr++){
            //Generate normal deviates and store in array randn
//            randn[rr] = GaussRNG.dev();
            randn[rr] = dist1(generator);
        }
    for (int jk = 0; jk < Ncells; jk++){
        //Generate normal deviates and store in array randn
//        randn[jk] = dist1(generator);

        Vcsol[jk] = Vcsol[jk] + h*Vc[jk] + eta[jk]; //D*sqrtdht*randn[jk];
        Vesol[jk] = Vesol[jk] + h*Ve[jk];
        ILsol[jk] = ILsol[jk] + h*IL[jk];


        Vc[jk] = Vcsol[jk];
        Ve[jk] = Vesol[jk];
        IL[jk] = ILsol[jk];

        eta[jk] =eta[jk] - htauc*eta[jk] + sqrtdht*randn[jk];

    }

    rhs_ode(Ncells,Vc,Ve,IL,par);

    //Make sure we do not save the transient
    if (iter >= (trans-1)*N){
        for (int aa = 0; aa < Ncells; aa++){
            int index = (iter-(trans-1)*N) + aa*N;
            solution_noise[index] = Vesol[aa];
        }
    }

}

//-----------------------Noiseless Integration---------------------------------

for (int jk1 = 0; jk1 < Ncells; jk1++){
        Vc[jk1] = Vcsol2[jk1];
        Ve[jk1] = Vesol2[jk1];
        IL[jk1] = ILsol2[jk1];
}

rhs_ode(Ncells,Vc,Ve,IL,par);

for (int iter = 0; iter < trans*N; iter++){
    for (int jk = 0; jk < Ncells; jk++){

        Vcsol2[jk] = Vcsol2[jk] + h*Vc[jk];
        Vesol2[jk] = Vesol2[jk] + h*Ve[jk];
        ILsol2[jk] = ILsol2[jk] + h*IL[jk];

        Vc[jk] = Vcsol2[jk];
        Ve[jk] = Vesol2[jk];
        IL[jk] = ILsol2[jk];

    }

    rhs_ode(Ncells,Vc,Ve,IL,par);

    //Make sure we do not save the transient
    if (iter >= (trans-1)*N){
        for (int aa = 0; aa < Ncells; aa++){
            int index = (iter-(trans-1)*N) + aa*N;
            solution_data[index] = Vesol2[aa];
        }
    }

}

delete[] randn;

}

////-----------------------absolute value of num1---------------------
//double dabs(double num1){
//
//    if(num1 >= 0.0){
//        return num1;
//    }
//    else{ return -num1;}
//}


//------------------------Average of an array-----------------------

double sol_avg(int N, double u1[]){

double avg_sum = 0;
double avg = 0;
int counter = 0;

for (int aa = 0; aa < N; aa++){
    if((std::isnan(u1[aa])==0) && (u1[aa]!=500)){//&&(u1[aa]>1e-15)){
            avg_sum = avg_sum+ u1[aa];
            counter++;
    }
//avg_sum = avg_sum + u1[aa];
}
if (counter != N){
        avg = avg_sum/counter;
}
else{
    avg = avg_sum/N;
}

return avg;

}


//------------------------------Absolute Deviations-------------------------

double abs_dev(int S, double u1[]){
    double avg1 = sol_avg(S,u1);
    double abs_diff[S];
    double temp2 = 0;
    for (int kk = 0; kk < S; kk++){
        temp2 = u1[kk] - avg1;
        abs_diff[kk] = abs(temp2);
    }

    double deviations = sol_avg(S,abs_diff);
    return deviations;

}

//------------------------Shifts solution by average------------------


void sol_shift(int N, double u1[]){
    double shift = sol_avg(N,u1);
    for (int ss = 0; ss < N; ss++){
        u1[ss] = u1[ss] - shift;
    }
}

//------------------------removes Standing Wave solutions-------------

int sol_class(int Ncells, int N, double solution_noise[]){

int beak = 5000;

double* u1;
u1 = new double[beak];
double* utemp;
utemp = new double[beak]();
//double* usum;
//usum = new double[beak]();

double classvec1[Ncells-1];

for (int zz = 0; zz< Ncells-1; zz++){
    classvec1[zz] = 0;
}

for (int ii = 0; ii < beak; ii++){                                //first oscillator solution
    u1[ii] = solution_noise[ii];
}

for (int jj = 1; jj < Ncells; jj++){
//    cout << "oops" << endl;
    for (int kk = 0; kk < beak; kk++){
        utemp[kk] = abs(u1[kk] - solution_noise[jj*N + kk]);   //diff between oscillators
    }

//    for (int tt = 0; tt < Ncells-1; tt++){                     //average difference
        classvec1[jj-1] = sol_avg(beak,utemp);
//    }

    for(int zz = 0; zz < beak; zz++){                              //resets utemp to 0
        utemp[zz] = 0;
    }
}

int bcount = 0;

for (int pp = 0; pp < Ncells-1; pp++){
    if (classvec1[pp] < 1e-2) bcount++;
}

//delete[] usum;
delete[] utemp;
delete[] u1;

if (bcount > 0) return 1;
else return 0;

}

//-------------------------Calculates Phase---------------------------

void phase_calc(int N,double t[], double u1[], double phase_data[], int cell_ind){


    double* p;                                      //initialize zero crossing vector
    p= new double[N];
        for(int jj=0;jj<N;jj++){                    //sets to 0
            p[jj]=0;
        }

    double slope=0;
    double a0=0;
    double a1=0;
    double a2=0;
    double A=0;
    double B=0;
    double C=0;
    double temp1=0;
    double temp2=0;

    double tol = 1.0e-14;
    int num1=0;

    int iter  = 0;
    int liter = 0;
    //-------------locate zero crossings using quadratic interpolation--------------
		for(int ii=1; ii<N-1;ii++){
            if((u1[ii]<0) && (0<u1[ii+1])){

                //Quadratic Interpolation
                a0 = u1[ii-1]/((t[ii-1]-t[ii])*(t[ii-1]-t[ii+1]));
                a1 = u1[ii]/((t[ii]-t[ii-1])*(t[ii]-t[ii+1]));
                a2 = u1[ii+1]/((t[ii+1]-t[ii])*(t[ii+1]-t[ii-1]));
 				 A = a0 + a1 + a2;
				 B = -((t[ii]+t[ii+1])*a0 + (t[ii-1]+t[ii+1])*a1 + (t[ii-1]+t[ii])*a2);
				 C = a0*t[ii]*t[ii+1] + a1*t[ii-1]*t[ii+1] + a2*t[ii-1]*t[ii];
                temp1 = (-B - sqrt(pow(B,2)-4*A*C))/(2*A);
                temp2 = (-B + sqrt(pow(B,2)-4*A*C))/(2*A);
                if(temp1 < t[ii+1] && temp1 > t[ii]){
                	p[num1] = temp1;
                	iter++;
				}
				else if(temp2 < t[ii+1] && temp2 > t[ii]){
					p[num1] = temp2;
					iter++;
				}
				else{
						slope = (u1[ii+1]-u1[ii-1])/(t[ii+1]-t[ii-1]);
                		p[num1]=t[ii]-u1[ii]/slope;     //linear interpolation
                        liter++;
				}
                num1++;
            }
        }


        double zero_crossing[num1];
        for (int zz = 0; zz < num1 ; zz++){
            zero_crossing[zz] = 0;
        }

        int num2 = 0;


        for (int kk = 0; kk < num1; kk++){   //p array with only crossing calculations
                zero_crossing[kk] = p[kk];
        }




        double period_array[num1];
        for (int zz = 0; zz < num1 ; zz++){
            period_array[zz] = 0;
        }

        double temp = 0;
        int num=0;

        for (int mm = 0; mm < num1-1; mm++){              //calculates period
            temp = zero_crossing[mm+1] - zero_crossing[mm];
			period_array[mm]= temp;
        }


        double phasetemp[600];
        for (int zz = 0; zz < 600 ; zz++){
            phasetemp[zz] = 0;
        }

        int badcounter = 0;
        for(int ii=0;ii<600;ii++){                         //saves first 600 periods
                phasetemp[ii]=period_array[ii];
                if (period_array[ii] < 1e-10) badcounter++;
			}

			if ((cell_ind >= 6)&&(badcounter!=0)){
                cout << "bad count: " << badcounter << endl;
			}

        double avg = 0;
        avg = sol_avg(600, phasetemp);

        for (int jj=0; jj<600; jj++){                    //saves data for phase calc
//            phase_data[jj] = abs(phasetemp[jj]-avg) /avg;  //check
              phase_data[jj] = phasetemp[jj];

        }



        delete[] p;

}




