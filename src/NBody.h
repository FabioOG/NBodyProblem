#ifndef NCORPS_H
#define NCORPS_H

#include <vector>
#include <utility>
#include "body.h"
#include "vector3.h"
#include <string>
#include <cmath>

enum class Integrator {
    Euler,
    Heun,
    Taylor,
    RungeKutta,
    Symplectic4,
    Symplectic6
};

enum class StepType {
    FixedStep,
    AdaptativeStep,
};

struct ComputationMethod {
    Integrator integrator;
    StepType stepType;
    bool energie;
    unsigned int order;
    std::vector<std::vector<long double>> RKCoeffDerivatives;
    std::vector<long double> RKCoeffAverage;
};

const ComputationMethod Euler = {Integrator::Euler,StepType::FixedStep,false};
const ComputationMethod Heun = {Integrator::Heun,StepType::FixedStep,false};
const ComputationMethod Bricolage = {Integrator::Euler,StepType::FixedStep,true};
const ComputationMethod Taylor4 = {Integrator::Taylor,StepType::FixedStep,false,4};
const ComputationMethod Taylor16 = {Integrator::Taylor,StepType::FixedStep,false,16};
const ComputationMethod RungeKutta4 = {Integrator::RungeKutta,StepType::FixedStep,false,4,\
        {{0,0,0,0},\
        {0.5,0,0,0},\
        {0,0.5,0,0},\
        {0,0,1.,0}},\
    {1./6.,1./3.,1./3.,1./6.}
    };
const ComputationMethod RungeKutta4_38 = {Integrator::RungeKutta,StepType::FixedStep,false,4,\
        {{0,0,0,0},\
        {1./3.,0,0,0},\
        {-1./3.,1.,0,0},\
        {1,-1.,1.,0}},\
    {1./8.,3./8.,3./8.,1./8.}
    };
const ComputationMethod RungeKutta5 = {Integrator::RungeKutta,StepType::FixedStep,false,6,\
        {{},\
        {0.2},\
        {3./40.,9./40.},\
        {3./10.,-9./10.,6./5.},\
        {-11./54.,5./2.,-70/27,35./27.},\
        {1631./55296.,175./512.,575./13824.,44275./110592.,253./4096.}},\
    {37./378.,0,250./621.,125./594.,0,512./1771.}
    };
const ComputationMethod Symplectic4 = {Integrator::Symplectic4,StepType::FixedStep,false};


class NBody
{
    public:

        NBody();
        std::vector<Body> body;

        void tick(long double delta);
        void tick();

        long double getStep();
        long double normalStep;
        long double stepFactor;
        long double stepPower;
        ComputationMethod computationMethod;
        long double positionPulseRatio;

        void Compute(long double tmax, unsigned int skip);
        std::vector<std::pair<vector3l,vector3l>> getResults();
        long double getRecordedStep(int i=1);
        unsigned int recordedTrajectoryIndex;
        bool getComputationDone();
        void save(std::string filename="trajectoire.csv");

        long double getEnergie();
        long double accelerationSum();
        void linkAll();

        unsigned int binomialCoeff(unsigned int n, unsigned int k);
        unsigned int factorial(unsigned int n);
        static void initialize();
        static const unsigned int maxn = 32;

    protected:
        std::vector< std::pair<long double,std::vector<  std::pair<vector3l,vector3l> >> > trajectories; // Choix du temps, puis du couple de vecteur, puis position/vitesse
        bool computationDone;
        static bool initialized;
        long double initEnergy;
        static unsigned int C[maxn + 1][maxn + 1];
        static unsigned int f[maxn+1];
        static long double fInv[maxn+1];
};

#endif // NCORPS_H
