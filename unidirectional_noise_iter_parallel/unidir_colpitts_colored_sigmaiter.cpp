//Iterates through D (noise signal strength) and calculates phase drift of Networks of Coupled Colpitts Oscillator written by Horacio Lopez, July 2019
//USE: unidir_txt_reader_colpitts.m to plot data.
//RUN AS: g++ --std=c++11 -fopenmp unidir_colpitts_colored_sigmaiter.cpp -o coloredSIGMA.exe
//        then run coloredSIGMA.exe
#include <random>
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctime>    // For time()
#include <sys/time.h>
#include <cstdlib>  // For srand() and rand()
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include "colpitts_euler_maruyama_color.h"
//#include "nr3.h"

//#include "deviates.h"
//#include "ran.h"
using namespace std;


void rhs_ode(int Ncells, double Vc[], double Ve[], double IL[], double par[12]);
void euler_maruyama(int N, int Ncells, int trans, double solution_data[], double solution_noise[], double eta[], double Vc[], double Ve[], double IL[], double par[12], double h, int seed, double np[]);
void sol_shift(int N,double u1[]);
void phase_calc(int N, double t[], double u1[], double phase_data[], int cell_ind);
double sol_avg(int N, double u1[]);
double abs_dev(int S, double u1[]);
int sol_class(int Ncells, int N, double solution_noise[]);



int main(){

//    clock_t time_s,time_e;
//    time_s = clock();
//    srand ( time(0) );

ofstream myfile_tsN;
myfile_tsN.open("unidir_colpitts_sigmaiter.txt");

struct timeval start, end;
gettimeofday(&start, NULL);

omp_set_nested(0);
int workers = 20;
//int osc = 2;

const double Vcc = 5.0;
const double Rl  = 200.0;
const double L   = 21.8e-06;
const double C1  = 1.0e-09;
const double C2  = 1.0e-09;
const double Ree = 1000.0;
const double Vee = -5.0;

const double beta   = 200;
const double Ron    = 100;
const double Vth    = 0.75;
const double Vb     = 0;
//double lambda = 0;

double h    = 1.00e-10;          // time-step size
int N       = 10*pow(2,20);        // steps to take
int trans   = 10;               //transient  = (N-1)*trans; samples at last iteration
long long N_iter;               //pre-allocation value for system integration




const int cellmax    = 5;           // number of oscillator cases you want to look at
                                    //i.e. cellmax = 2 => Ncells = [3, 5]
const int sample_max = 20;           // number of samples to take

const int lambda_size = 10;         // grid size for lambda sweeping

double lam_min = -0.05;              //min lambda

double lam_max = 0;                 //max lambda

double lam_step = abs(lam_max-lam_min)/(lambda_size-1);

double lambda_space[lambda_size];

int sigma_size = 2;                 //grid size for D (or sigma) sweeping

double sig_min = 1e-13;             //min D value

double sig_max = 1e-12;             //max D value

double sig_step = abs(sig_max-sig_min)/(sigma_size-1);

double sigma_space[sigma_size];

double tau_c = 1.0e-3;              //noise correlation time





for (int ii=0; ii<lambda_size; ii++){                //fills out vector with lambdas
        lambda_space[ii] = lam_min + lam_step*ii;
}

for (int ii=0; ii<sigma_size; ii++){                 //fills out vector with D
        sigma_space[ii]  = sig_min + sig_step*ii;
}

double print_sigma[cellmax+1];                 //saves a vector with D values to print
for (int zz = 0; zz < cellmax + 1; zz++){
    print_sigma[zz] = 0;
}

for (int ii=0; ii<sigma_size; ii++){                 //in .txt file
        print_sigma[ii]  = sigma_space[ii];
}

for (int kk = 0; kk < cellmax+1; kk++){              //prints sigma space to .txt file
        myfile_tsN<<print_sigma[kk]<<'\t';
}
myfile_tsN<<endl;

double mu = 0.0,sigma = 1.0;

int seed = 91; //rand()%100;//11//4
Normaldev_BM GaussRNG(mu,sigma,seed);



double* t; //Time array
t = new double[N];

for (int bb = 1; bb < N+1; bb++){
    t[bb-1] = h*(bb-1);
    }

for (int sigiter = 0; sigiter < sigma_size; sigiter ++){ //iterates through D values

    double D = sigma_space[sigiter];                         //D value

    double np[2] = {D,tau_c};

    double phase_space[lambda_size][cellmax];                //saves results for printing

    for (int ii=0; ii<lambda_size; ii++){
        for (int jj=0; jj<cellmax; jj++){
            phase_space[ii][jj] = 0.0;
        }
    }

    int lamiter = 0;                                           //iterates through lambda values
        do {                                                   //same as a for-loop. no specific reason why it is not.
        //for (int lamiter = 0; lamiter < lambda_size; lamiter++){

        double lambda = lambda_space[lamiter];
        double par[] = {Vcc,Rl,L,C1,C2,Ree,Vee,beta,Ron,Vth,Vb,lambda};


        for (int celliter = 0; celliter < cellmax; celliter++){ //iterates through N oscillator cases

            int Ncells = 2*(celliter+1) + 1; //# of oscillators


            N_iter = Ncells*N;

            double phase_drift[sample_max];             //array to save sample results

            for (int zz = 0; zz < sample_max ; zz++){   //setting to 0
                phase_drift[zz] = 0;
            }

            int NotTWcount = 0;

            #pragma omp parallel for num_threads(workers)
            for (int samples=0; samples < sample_max; samples++){

                static thread_local std::mt19937 generator;
                std::uniform_int_distribution<int> distribution(0,1000000);

                //Initializing random initial conditions
                double* Vc;
                Vc = new double[Ncells]();
                double* Ve;
                Ve = new double[Ncells]();
                double* IL;
                IL = new double[Ncells]();

                for (int r1=0;r1<Ncells;r1++){
                    Vc[r1] = (distribution(generator))*1.0e-07;
                }

                for (int r2=0;r2<Ncells;r2++){
                    Ve[r2] = (distribution(generator))*1.0e-07;
                }

                for (int r3=0;r3<Ncells;r3++){
                    IL[r3] = (distribution(generator))*1.0e-07;
                }


                //----------------------Euler-Maruyama Integration (excludes transient)---------------
                double* solution_data;
                solution_data = new double[N_iter];

                double* solution_noise;
                solution_noise = new double[N_iter];


                double* eta; //pre-allocation for random variables
                eta = new double[N]();
                //for (int zz = 0; zz < N ; zz++){
                //    eta[zz] = 0;
                //}

                euler_maruyama(N,Ncells,trans,solution_data,solution_noise,eta,Vc,Ve,IL,par,h,seed,np);

                int type = sol_class(Ncells,N,solution_noise);

                if (type != 0) NotTWcount = NotTWcount + 1;

                //----------------------Phase-Drift calculations--------------------
                int hh     = 600;                           //number of periods to capture per oscillator

                double* phase;                              //initialize dynamic array for phase average
                phase       = new double[hh]();

                if (type != 0) phase_drift[samples] = 500;

                else{

                    for (int dd=0; dd<Ncells; dd++){
                        double* phase_data;                     //initialize dynamic array for phase calc
                        phase_data  = new double[hh]();

                        double* u1;
                        u1 = new double[N]();

                        double* phase1;                         //initialized vectors for computing difference
                        phase1      = new double[hh]();         //between solution with noise and without
                        double* phase2;
                        phase2      = new double[hh]();
                        double* phase3;
                        phase3      = new double[hh]();

                        for (int d1=0; d1 < N; d1++){           //look at individual, noiseless oscillator
                            u1[d1] = solution_data[d1 + dd*N];

                        }

                        //sol_shift(N,u1);                        //remove average trend
                        phase_calc(N,t,u1,phase_data,dd);       //calculate phase vector

                        for (int d4 = 0; d4 < 600; d4++){
                            phase1[d4] = phase_data[d4];
                        }

                        for (int zz = 0; zz < N ; zz++){        //reset values manually
                            u1[zz] = 0;
                        }

                        for (int zz = 0; zz < hh ; zz++){       //reset values manually
                            phase_data[zz] = 0;
                        }

                        for (int d5=0; d5 < N; d5++){           //look at individual, noisy oscillator
                            u1[d5] = solution_noise[d5 + dd*N];

                        }

                        //sol_shift(N,u1);
                        phase_calc(N,t,u1,phase_data,dd);

                        for (int d6 = 0; d6 < 600; d6++){
                            phase2[d6] = phase_data[d6];
                        }

                        for (int d7 = 0; d7 < 600; d7++){       //calculate difference between noise and noiseless
                            phase3[d7] = abs(phase1[d7]-phase2[d7]);
                        }

                        for (int d2=0; d2<600;d2++){            //sum phase difference between noise and noiseless
                            phase[d2] = phase[d2] + phase3[d2];
                        }
                        //delete[] phase;
                        delete[] phase1;                            //deallocate dynamic memory
                        delete[] phase2;
                        delete[] phase3;
                        delete[] u1;
                        delete[] phase_data;
                    }

                    for (int d3=0; d3 < 600; d3++){             //calculate average
                        phase[d3] = phase[d3]/Ncells;
                    }
                    phase_drift[samples] = abs_dev(600,phase);  //calculate absolute deviations

                }

                delete[] Vc;          //deletes dynamic variables
                delete[] Ve;
                delete[] IL;
                delete[] phase;
                //delete[] phase1;
                //delete[] phase2;
                //delete[] phase3;
                //delete[] phase_data;
                delete[] eta;
                delete[] solution_data;
                delete[] solution_noise;
                //delete[] u1;
            }


            double slope = sol_avg(sample_max,phase_drift);   //average of samples

            phase_space[lamiter][celliter] = slope;           //saves average of samples

        }

        cout << "Sigma Iteration: " << sigiter+1 << " of " << sigma_size << "; Lambda Iteration: " << lamiter+1 << " of " << lambda_size << endl;
        lamiter++;

        } while(lamiter < lambda_size);

    //-------------------writes results to .txt file---------------------------------

    for (int kk = 0; kk < lambda_size; kk++){
            myfile_tsN<<lambda_space[kk]<<'\t';
        for(int i = 0;i<cellmax;i++){
                myfile_tsN<<phase_space[kk][i]<<'\t';
            }
             myfile_tsN<<endl;
    }


}
myfile_tsN.close();

gettimeofday(&end, NULL);

double delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
         end.tv_usec - start.tv_usec) / 1.e6;

cout << "TIME: " << delta << endl;

delete[] t;

}


