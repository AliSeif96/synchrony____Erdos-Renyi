/************************************************************************************************/
/*** Topic: Hodgkin-Huxley model with Runge-Kutta 4th Order Method in Erdos-Renyi network     ***/
/***	    for N neurons and Possibility of connection P                            Ali-Seif ***/
/*** Version Release 17.12 rev 11256                                                          ***/
/*** Date: 12/5/2020                                                                          ***/
/*** Code implemented in CodeBlocks C++ compiler (v. 17.12),                                  ***/
/*** MSI: PX60 6QD/ DDR4                                                                      ***/
/*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
/************************************************************************************************/
/*
The following an article have been used to write this program.
[1] Dynamic range in small-world networks of Hodgkin–Huxley neurons with chemical synapses
C.A.S. Batista a , R.L. Viana a, S.R. Lopes a, A.M. Batista
https://doi.org/10.1016/j.physa.2014.05.069
0378-4371/© 2014 Published by Elsevier B.V.
*/
#include <iostream>                                             //for cout
#include <math.h>                                               //for pow()
#include <fstream>                                              //for ofstream
#include <stdlib.h>                                             //srand, rand
#include <time.h>                                               //time
#include <sstream>   //for stringstream
using namespace std;                                            //for Standard program
//##############################################################
//####                                                      ####
//####                       Variables                      ####
//####                                                      ####
//##############################################################
#define size1 100                                                 //Number of neurons
#define probability 0.3                                         //Possibility of connection
#define endtime 2500                                              //finaltime for run
#define current 2                                             //External Current
#define lengthsteps 0.01                                        //length steps
//##############################################################
//####                                                      ####
//####                 Create class Neuron                  ####
//####                                                      ####
//##############################################################
class Neuron
{
	private:                                                    //Fixed and variable values that we do not have access to outside the class
		double Gna = 120.0;                                     //sodiom_constant
		double Gk = 36.0;                                       //Potasiom_constant
		double Gl = 0.3;                                        //Leak_constant
		double Gsyn = 0.06;                                        //Leak_constant
        double Ena = 55.17;                                     //Votage_sodiom
		double Ek = -72.14;                                     //Votage_Potasiom
		double El = -49.42;                                     //Votage_Leak
		double Esyn = 20.0;                                     //Votage_Leak
        double Cm = 1.0;                                        //Capacitor_capacity
        double Iapp= current;                                   //External Current
        double  k1, k2, k3, k4;                                 //define 4 point for calculate Runge-Kutta
	public:                                                     //Fixed and variable values that we have access to from outside the class
	    ofstream temp;                                          //create file for save data
	    void start(void);                                       //call start and Initial values
        void onedt(void);                                       //Run for one step dt
        void runprint(void);                                    //run for all time and print in file
        void runshow(void);                                     //run for all time and show in exe
        int num=size1;                                          //Number of neurons
		double alpha_n(double);                                 //calculate alpha n
		double beta_n(double);                                  //calculate beta n
	  	double alpha_m(double);                                 //calculate alpha m
	  	double beta_m(double);                                  //calculate beta m
	  	double alpha_h(double);                                 //calculate alpha h
	  	double beta_h(double);                                  //calculate beta h
	  	double alpha_r(double);                                 //calculate alpha r
	  	double beta_r(double);                                  //calculate beta r
	  	double n_inf(double);                                   //calculate infinitude n
	  	double m_inf(double);                                   //calculate infinitude m
	  	double h_inf(double);                                   //calculate infinitude h
        double r_inf(double);                                   //calculate infinitude r
	  	double INa1(double,double,double);                      //sodiom Current
	  	double IK1(double,double);                              //potasiom Current
	  	double Il1(double);                                     //leak Current
	  	double Isyn(double,int);                                     //leak Current
	  	double dvdt(double,double,double,double,double,int);        //Differential equation for voltage
	  	double dndt(double,double,double);                      //Differential equation for n
        double dmdt(double,double,double);                      //Differential equation for m
        double dhdt(double,double,double);                      //Differential equation for h
        double drdt(double,double,double);                      //Differential equation for r
	  	double rk4thOrder_v(double,double,double,double,double,double,int);//Runge-Kutta for voltage
	  	double rk4thOrder_n(double,double,double,double);       //Runge-Kutta for n
        double rk4thOrder_m(double,double,double,double);       //Runge-Kutta for m
	  	double rk4thOrder_h(double,double,double,double);       //Runge-Kutta for h
        double rk4thOrder_r(double,double,double,double);       //Runge-Kutta for r

        float *fillarr(float);
        float *readfile(int);



        float arr[];
        double t_final=endtime;                                 //final time
        double dt = lengthsteps;                                //length steps
        double V[size1+1][2+1];                                 //define matrix for pre and post voltage of any neuron
        double N[size1+1][2+1];                                 //define matrix for pre and post n of any neuron
        double H[size1+1][2+1];                                 //define matrix for pre and post h of any neuron
        double M[size1+1][2+1];                                 //define matrix for pre and post m of any neuron
        double R[size1+1][2+1];                                 //define matrix for pre and post r of any neuron
        double A[size1+1][size1+1];                              //Because matrix start from 0
	  	double t0;                                              //Pedometer
	  	int conter;
  		double v,n,m,h,r;};                                       //variable values for each step

//_________________________________Calculate alpha and betas_________________________________//

double Neuron::alpha_n(double v){return   0.01*(v+50)/(1-exp(-(v+50)/10));}
double Neuron::beta_n(double v){return   0.125*exp(-(v+60)/80);}
double Neuron::alpha_m(double v){return   0.1*(v+35)/(1-exp(-(v+35)/10));}
double Neuron::beta_m(double v){return   4.0*exp(-0.0556*(v+60));}
double Neuron::alpha_h(double v){return   0.07*exp(-0.05*(v+60));}
double Neuron::beta_h(double v){return   1/(1+exp(-(0.1)*(v+30)));}
double Neuron::alpha_r(double v){return   1.875/(1+exp(-(v+20)));}
double Neuron::beta_r(double v){return   0.125;}

//__________________________Calculate infinite activation variables__________________________//

double Neuron::n_inf(double v){return alpha_n(v)/(alpha_n(v)+beta_n(v));}
double Neuron::h_inf(double v){return alpha_h(v)/(alpha_h(v)+beta_h(v));}
double Neuron::m_inf(double v){return alpha_m(v)/(alpha_m(v)+beta_m(v));}
double Neuron::r_inf(double v){return alpha_r(v)/(alpha_r(v)+beta_r(v));}

//__________________________________Calculation of currents__________________________________//

double Neuron::INa1(double v,double h,double m) {return Gna*h*pow(m,3)*(v-Ena);}
double Neuron::IK1(double v,double n) {return Gk*pow(n,4)*(v-Ek);}
double Neuron::Il1(double v) {return Gl*(v-El);}

double Neuron::Isyn(double v,int conter) {

    double sigma=0;
    for(int i=1; i<=num ; i++){

            sigma=sigma+A[conter][i]*R[i][1]*(Esyn-V[i][1]);
            //cout<<"Isyn"<<i<<'\t'<<A[conter][i]<<endl;
    }
    return Gsyn*sigma;
    }

//___________________________________Differential Equations__________________________________//

double Neuron::dvdt(double t, double v,double n,double h,double m,int conter){return  (1/Cm)*(Iapp +Isyn(v,conter)-(INa1(v,h,m)+IK1(v,n)+Il1(v)));}
double Neuron::dndt(double t,double n, double v){return  ((alpha_n(v)*(1-n))-beta_n(v)*n);}
double Neuron::dhdt(double t, double h, double v){return   ((alpha_h(v)*(1-h))-beta_h(v)*h);}
double Neuron::dmdt(double t, double m, double v){return   ((alpha_m(v)*(1-m))-beta_m(v)*m);}
double Neuron::drdt(double t, double r, double v){return   ((alpha_r(v)*(1-r))-beta_r(v)*r);}

//__________________________________Runge-Kutta calculations_________________________________//

double Neuron::rk4thOrder_v(double t0, double v, double dt,double n,double h,double m,int conter) {
            k1=     dt*dvdt(t0, v,n,h,m,conter);
            k2=     dt*dvdt((t0+dt/2), (v+k1/2),n,h,m,conter);
            k3=     dt*dvdt((t0+dt/2), (v+k2/2),n,h,m,conter);
            k4=     dt*dvdt((t0+dt), (v+k3),n,h,m,conter);
            v=      v+double((1.0/6.0)*(k1+2*k2+2*k3+k4));
   return   v;}
double Neuron::rk4thOrder_n(double t0, double v, double dt, double n) {
            k1=     dt*dndt(t0, n,v);
            k2=     dt*dndt((t0+dt/2), (n+k1/2),v);
            k3=     dt*dndt((t0+dt/2), (n+k2/2),v);
            k4=     dt*dndt((t0+dt), (n+k3),v);
            n=      n+double((1.0/6.0)*(k1+2*k2+2*k3+k4));
   return   n;}
double Neuron::rk4thOrder_h(double t0, double v, double dt,double h) {
            k1=     dt*dhdt(t0, h,v);
            k2=     dt*dhdt((t0+dt/2), (h+k1/2),v);
            k3=     dt*dhdt((t0+dt/2), (h+k2/2),v);
            k4=     dt*dhdt((t0+dt), (h+k3),v);
            h=      h+double((1.0/6.0)*(k1+2*k2+2*k3+k4));
    return   h;}
double Neuron::rk4thOrder_m(double t0, double v, double dt,double m) {
            k1=     dt*dmdt(t0, m,v);
            k2=     dt*dmdt((t0+dt/2), (m+k1/2),v);
            k3=     dt*dmdt((t0+dt/2), (m+k2/2),v);
            k4=     dt*dmdt((t0+dt), (m+k3),v);
            m=      m+double((1.0/6.0)*(k1+2*k2+2*k3+k4));
    return   m;}

double Neuron::rk4thOrder_r(double t0, double v, double dt,double r) {
            k1=     dt*dmdt(t0, r,v);
            k2=     dt*dmdt((t0+dt/2), (r+k1/2),v);
            k3=     dt*dmdt((t0+dt/2), (r+k2/2),v);
            k4=     dt*dmdt((t0+dt), (r+k3),v);
            r=      r+double((1.0/6.0)*(k1+2*k2+2*k3+k4));
    return   r;}
//_________________________________One step run calculations_________________________________//

void Neuron::start(){

    for (int i=1 ; i<=num;i++){
       for (int j=1 ; j<3;j++){
            V[i][j]=0;
            N[i][j]=0;
            H[i][j]=0;
            M[i][j]=0;
            R[i][j]=0;}}
    v=-20.0;
    for(int i=1;i<=num;i++){
        V[i][1]=v;
        N[i][1]=n_inf(v);
        H[i][1]=h_inf(v);
        M[i][1]=m_inf(v);
        R[i][1]=r_inf(v);}
}

void Neuron::onedt(){
    v=rk4thOrder_v(t0, v, dt,n,h,m,conter);
    n=rk4thOrder_n(t0,v, dt ,n);
    h=rk4thOrder_h(t0, v, dt ,h);
    m=rk4thOrder_m(t0, v, dt ,m);
    r=rk4thOrder_r(t0, v, dt ,r);
    }

//__________________________run for all steps and print in file _____________________________//
void Neuron::runprint(){
    //ofstream temp;
    ofstream temp("temp.txt");
    for (t0=dt ; t0<=t_final ;t0=t0 + dt){
    conter=0;
        for ( int i=1 ; i<=num ; i++){
            conter=conter+1;
            v=V[i][1];
            n=N[i][1];
            h=H[i][1];
            m=M[i][1];
            r=R[i][1];
            onedt();

            V[i][2]=v;
            N[i][2]=n;
            H[i][2]=h;
            M[i][2]=m;
            R[i][2]=r;

            V[i][1]=V[i][2];
            N[i][1]=N[i][2];
            H[i][1]=H[i][2];
            M[i][1]=M[i][2];
            R[i][1]=R[i][2];
        }
        if(t0>=2000.0){
        double Rsum=0;
        temp<<t0;
        for(int k=1;k<=num;k++){
            Rsum= Rsum+ R[k][2];
            temp<<'\t'<<V[k][2];
        }
        temp<<'\t'<<Rsum<<endl;
        }
    }
    temp.close();
}
//___________________________ run for all steps and show in exe _____________________________//
void Neuron::runshow(){

    for (t0=dt ; t0<=t_final ;t0=t0 + dt){
        conter=0;
        for ( int i=1 ; i<=num ; i++){
            //cout<<endl<<i<<endl;
            conter=conter+1;
            v=V[i][1];
            n=N[i][1];
            h=H[i][1];
            m=M[i][1];
            r=R[i][1];

            onedt();

            V[i][2]=v;
            N[i][2]=n;
            H[i][2]=h;
            M[i][2]=m;
            R[i][2]=r;

            V[i][1]=V[i][2];
            N[i][1]=N[i][2];
            H[i][1]=H[i][2];
            M[i][1]=M[i][2];
            R[i][1]=R[i][2];
        }
        if(t0>=2000.0){
        cout<<t0;
        double Rsum=0;
        for(int k=1;k<=num;k++){
             Rsum= Rsum+ R[k][2];
            cout<<'\t'<<V[k][2];
        }
        cout<<'\t'<<Rsum<<endl;
        }
    }
}

//##############################################################
//####                                                      ####
//####            Create class Neighborhood                 ####
//####                                                      ####
//##############################################################

class Neighborhood
{
	private:                                                    //Fixed and variable values that we do not have access to outside the class
        float p=probability;                                    //Possibility of connection
        int iSecret;                                            //Random probability between zero and one
	public:                                                     //Fixed and variable values that we have access to from outside the class
	    ofstream temp1;                                         //create file for save data
	    ofstream temp2;                                         //create file for save data
        void oneA(void);                                        //create Adjacency matrix
        void show(void);                                        //show Adjacency matrix
        void print(void);                                       //print Adjacency matrix
        void bineryprint(void);                                       //print Adjacency matrix
        int n=size1;                                            //number of node
        double A[size1+1][size1+1];                              //Because matrix start from 0
        int numb;};                                             //Counting links
//________________________________ create Adjacency matrix __________________________________//

void Neighborhood::oneA(){
    srand (time(NULL));                                         //initialize random seed
    for (int i=1 ; i<=n;i++){
       for (int j=1 ; j<=n;j++){A[i][j]=0;}}
    numb=0;
    for (int i=1 ; i<=n;i++){
        for (int j=i ; j<=n;j++){
            iSecret = rand() % 10 + 1;
            int pp=p*10;
            if (iSecret<=pp && i!=j ){
                A[i][j]=1;
                numb=numb+1;
            }
            else{
               A[i][j]=0;
            }
            A[j][i]=A[i][j];}}}
//________________________________ show Adjacency matrix ____________________________________//

void Neighborhood::show(){
    for (int i=1 ; i<=n;i++){cout<<'\t'<<i;}
    cout<<endl<<endl;
    for (int i=1 ; i<=n;i++){
        cout<<i<<'\t';
        for (int j=1 ; j<=n;j++){cout<<A[i][j]<<'\t';}
        cout<<endl<<endl;}}

//_________________________________ print Adjacency matrix __________________________________//

 void Neighborhood::print(){
    ofstream temp1("datas.txt");
    temp1<<"{\"nodes\":[{\"name\":\"0\",\"group\":1}";
    for (int i=1 ; i<n;i++){temp1<<",{\"name\":\""<<i<<"\",\"group\":"<<i<<"}";}
    temp1<<"],\"links\":[";
    int number=0;
    for (int i=0 ; i<n;i++){
        for (int j=i ; j<n;j++){
            if (A[i+1][j+1]==1){
                number=number+1;
                temp1<<"{\"source\":"<<i<<",\"target\":"<<j<<",\"value\":1}";
                if(number<numb){temp1<<",";}}}}
    temp1<<"]}";
    temp1.close();
    cout << "\nFinish" << endl;}



 void Neighborhood::bineryprint(){
    ofstream temp2("bineryprint.txt");


    for (int i=1 ; i<=n;i++){
        for (int j=1 ; j<=n;j++){temp2<<A[i][j]<<'\t';}
        temp2<<endl;}


    temp2.close();
    cout << "\nFinish" << endl;}


//##############################################################
//####                                                      ####
//####                     fill array                       ####
//####                                                      ####
//##############################################################
/*float Neuron::*fillarr( float arr[] ) {
    return arr;
}
//##############################################################
//####                                                      ####
//####              Convert file to array                   ####
//####                                                      ####
//##############################################################
float Neuron::*readfile(int neurun ) {


    float Array[6];

    for(int i = 0; i < 6; ++i){
        Array[i]=i;

    }


    float *nn = fillarr(Array);
    return   nn;
}
*/


//_______________________________________________________________________________________\\
//_____________                                                             _____________\\
//_____________                                      @                      _____________\\
//_____________           @@       @@       @            @@     @           _____________\\
//_____________           @ @     @ @      @ @       @   @ @    @           _____________\\
//_____________           @  @   @  @     @   @      @   @  @   @           _____________\\
//_____________           @   @@@   @    @@@@@@@     @   @   @  @           _____________\\
//_____________           @    @    @   @       @    @   @    @ @           _____________\\
//_____________           @         @  @         @   @   @     @@           _____________\\
//_______________________________________________________________________________________

int main() {
    Neighborhood neig;                                          //call class Adjacency matrix
    neig.oneA();                                                //create Adjacency matrix
    //neig.show();                                                //show Adjacency matrix
    neig.print();                                               //print Adjacency matrix
    neig.bineryprint();



//______________________________________________________________

    Neuron neuron;                                              //call class Hodgkin-Huxley
    for (int i=1 ; i<=size1;i++){
        for (int j=1 ; j<=size1;j++){
            neuron.A[i][j]=neig.A[i][j];
                //cout<<neig.A[i][j]<<'\t'<<neuron.A[i][j]<<endl;
        }
    }





    neuron.start();                                             //start and create variable
    //neuron.runshow();                                           //run Hodgkin-Huxley and show
    neuron.runprint();                                          //run Hodgkin-Huxley and print

    cout << "\nFinish" << endl;
    return 0;
}
