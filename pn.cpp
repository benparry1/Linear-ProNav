#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>


int main() {

    //initial engagement condition
    long double aT = 0;
    double HE = 10* M_PI / 180;
    int N = 3;                          //in general 2<n<5 generates lead behavior
    double h = 1;

    long double beta_rad = 180 * M_PI / 180;
    long double RTX = 30000;    //m
    long double RTZ = 10000;    //m

    long double RMX = 0;    //m
    long double RMZ = 0;    //m

    long double VM = 15000;   //m/s
    long double VT = 10020;        //m/s

    //relative position components and magnitude:
    long double R_TMX = RTX - RMX;
    long double R_TMZ = RTZ - RMZ;
    long double RTM = sqrt(pow(R_TMX, 2) + pow(R_TMZ, 2));      //relative position magnitude

    //line of sight angle and lead angle:
    long double lambda_r = atan(R_TMZ / R_TMX);
    long double L = (asin(VT * sin(beta_rad + lambda_r) / VM)) - HE;

    //missile and target velocity components:
    long double VMX = VM * cos(lambda_r + L + HE);
    long double VMZ = VM * sin(lambda_r + L + HE);

    long double VTX = -VT * cos(beta_rad);
    long double VTZ = VT * sin(beta_rad);

    //relative velocity components:
    long double V_TMX = VTX - VMX;
    long double V_TMZ = VTZ - VMZ;

    //change of LOS angle in time
    long double lambda_dot = ((R_TMX * V_TMZ) - (R_TMZ * V_TMX)) / pow(RTM, 2);

    //closing velocity:
    long double nc;
    long double VC = -((R_TMX * V_TMX) - (R_TMZ * V_TMZ)) / RTM;

    //intercept variables and loop conditions initialization
    double t = 0;
    double exit_time = 0;
    long double Txint, Tzint, Mxint, Mzint;
    bool intercept = false;
    bool ground = false;
    int count=0;

    long double y0[9] ={aT, HE, beta_rad, RTX, RTZ, RMX, RMZ, VT, VM };


    std::ofstream os;
    os.open("missile1.dat");

    //runtime start
    auto start = std::chrono::high_resolution_clock::now();

    while (!intercept && !ground) {

        os<<t<<"\t\t"<<RTX<<'\t'<<RTZ<<"\t\t"<<RMX<<'\t'<<RMZ<<'\n';

        t = t + h;

        //forward euler method using each functions time derivative (ie beta_dot, RTX_dot ...)
        //y[n+1]=y[n]+h*f_dot
        beta_rad += h * (aT / VT);
        RTX += h * (-VT * cos(beta_rad));
        RTZ += h * (VT * sin(beta_rad));
        RMX += h * (VMX);
        RMZ += h * (VMZ);
        VTX += h * (aT * sin(beta_rad));
        VTZ += h * (aT * cos(beta_rad));

        //closing velocity adjustments
        nc = N * VC * lambda_dot;
        VMX += h * (-nc * sin(lambda_r));
        VMZ += h * (nc * cos(lambda_r));

        std::cout << std::fixed;
        std::cout.precision(6);

        std::cout << t << '\t' << beta_rad << '\t' << RTX << '\t' << RTZ << '\t' << RMX << '\t' << RMZ << '\t' << VTX << '\t'
                  << VTZ << '\t' << VMX << '\t' << VMZ << '\n';

        //recalculate values using n+1 variables
        VT = sqrt(pow(VTX, 2) + pow(VTZ, 2));

        R_TMX = RTX - RMX;
        R_TMZ = RTZ - RMZ;
        RTM = sqrt(pow(R_TMX, 2) + pow(R_TMZ, 2));

        V_TMX = VTX - VMX;
        V_TMZ = VTZ - VMZ;

        lambda_r = atan(R_TMZ / R_TMX);
        lambda_dot = ((R_TMX * V_TMZ) - (R_TMZ * V_TMX)) / pow(RTM, 2);

        VC = -((R_TMX * V_TMX) - (R_TMZ * V_TMZ)) / RTM;

        //to reduce the number of computations and speedup runtime, the time step decreases as the target and pursuer
        //grow closer. This means that less time steps are used to reach an engagement and then once they are close
        // enough, a smaller time step is used to ensure the target and pursuer intercept.
        if(RTM>=100000)
            h=.5;
        else if (RTM < 100000 && RTM>= 50000) {
            h=.1;
        }
        else if(RTM<50000 && RTM>=10000)
            h=.01;
        else if(RTM<10000 && RTM>=5000)
            h=.001;
        else if(RTM<5000 && RTM>=1000)
            h=.0001;
        else if(RTM<1000 && RTM>=100)
            h=.00001;
        else if (RTM<=.1)           //if the distance between the two are within .1 m, intercept condition is
        {                           //set to true and the loop ends.
            exit_time = t;
            Txint = RTX;
            Tzint = RTZ;
            Mxint = RMX;
            Mzint = RMZ;
            intercept = true;
        }
        if (RTZ < 0 || RMZ<0) {         //if the target or pursuer hit the ground, the loop ends as well
            exit_time = t;
            ground = true;
        }
        count++;                        //records how many times the while loop runs

    }
    auto stop = std::chrono::high_resolution_clock::now();          //runtime stop
    os.close();

    //to get duration subtract time points.
    //to get proper units use duration cast method
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "\nTime taken by function: "
              << duration.count() << " milliseconds\nloops completed: " <<count<<'\n';

    std::cout << "\nrun conditions: {" << y0[0] << ", " << y0[1] << ", " << N << ", " << y0[2]
              << ", " << y0[3] << ", " << y0[4] << ", " << y0[5] << ", " << y0[6] << ", " << y0[7] << ", " << y0[8] << "}\n\n";

    if(ground)
        std::cout<<"missile hit the ground at time: " << exit_time << '\n';
    else
    {
        std::cout<<"target was intercepted at time: "<<exit_time<<" target pos: ("<<Txint<<", "<<Tzint<<") missile pos: ("<<Mxint<<", "<<Mzint<<")\n";
    }

    return 0;
}


//formulas and variables needed:

//kinematics of the target:

//target evasive acceleration: a_T=beta_dot * x-comp of target velocity vector

//beta dot is the time rate of change of the target's heading angle

// x comp of target velocity in inertial coord: V^i_Tx = -V^t_Tx * cos(beta)
// z comp of the target velocity in the inertial cood: V^i_Tz = V^t_Tz * sin(beta)

//x compt of time rate of chane of the range in inertial coord system: Rdot^i_Tx = V^i_Tx
//z comp (upward) of the time rate of change of the range in inertial cood system: Rdot^i_Tz = V^i_Tz

//x comp of the time rate of change of the velocity in inertial coord: Vdot^i_Tx = a^t_Tz * sin(beta)
//z comp of the time rate of change of the velocity in inertial coord: Vdot^i_Tz = a^t_Tz * cos(beta)

//kinematics of the pursuer:

//x comp of pursuer acc: a^i_Px = -a^p_Pz * sin(lambda+L+HE)
//z comp of pursuer acc: a^i_Pz= a^p_Pz * cos(lambda+L+HE)

//x compt of time rate of change of the pursuer range in inertial coord system: Rdot^i_Px = V^i_Px
//z compt of time rate of change of the pursuer range in inertial coord system: Rdot^i_Pz = V^i_Pz

//x comp of the time rate of change of the pursuer velocity: Vdot^i_Px = a^i_Px
//z comp of the time rate of change of the pursuer velocity: Vdot^i_Pz = a^i_Pz


//relative position and velocity:

//relative position (vector): R^i_T/P = R^i_P - R^i_T
//relative velocity (vecotr): V^i_T/P = V^i_P - R^i_T

//closing velocity: V_c = - (R^i_TPx * V^i_TPx) - (R^i_TPz * V^i_TPz) / |R^i_T/P|

//line of sight angle and rate:

//line of sight (lambda): lambda = tan^-1(R^i_TPz/R^i_TPx)
//time rate of change of line of sight: lambda^dot = (R^i_TPx * V^i_TPz) - (R^i_TPz * V^i_TPx) / |R^i_T/P|^2

//lead angle: L = sin^-1[sin(beta + lambda) * V_T /V_P] - HE


//time rate of change of beta: beta^dot = a^t_T / V^t_Tx

//engagement conditions/parameters:

//pursuer: L, HE, R^i_P, V_P
//target: beta, a_T, R^i_T, V_T
//navigation constant: N

//initial conditions:

//inertial positions: beta_0, R_Tx0, R_Tz0, R_Px0, R_Pz0,
//inertial velocities: V_Tx0, V_Tz0, V_Px0, V_Pz0

//ProNav: lambda_dot(t_i), V_c(t_i)

//a^Pure_P = N * V_c * lambda_dot
