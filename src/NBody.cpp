#include "NBody.h"

#include "vector3.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>

using namespace std;

const unsigned int NBody::maxn;
bool NBody::initialized = false;
unsigned int NBody::C[NBody::maxn+1][NBody::maxn+1] = {0};
unsigned int NBody::f[NBody::maxn+1] = {0};
long double NBody::fInv[NBody::maxn+1] = {0};

NBody::NBody()
{
    recordedTrajectoryIndex=0;
    stepPower = 0;
    stepFactor = 1;
    computationDone = false;
    computationMethod = RungeKutta4;
    positionPulseRatio = 1;
    if (!initialized) initialize();
}

void NBody::initialize()
{
    C[0][0] = 1;
    for (unsigned int n = 1; n < maxn; n++) {
        // Set every nCr = 1 where r = 0
        C[n][0] = 1;
        for (unsigned int k = 1; k < n + 1; k++) {
            // Value for the current cell of Pascal's triangle
            C[n][k] = (C[k - 1][k - 1] + C[n - 1][k]);
        }
    }

    f[0]=1;
    for (unsigned int i=1; i<=maxn; i++)
    {
        f[i] = i*f[i-1];
    }
    for (unsigned int i=1; i<=maxn; i++)
    {
        fInv[i] = 1/f[i];
    }

    initialized=true;
}

void NBody::Compute(long double tmax, unsigned int skip)
{
    initEnergy = getEnergie(); // Usefull for energy conservation
    stepFactor = normalStep * (long double)pow((long double)accelerationSum(),stepPower);
    computationDone=false;
    long double t = 0;
    trajectories = {};
    vector<pair<vector3l,vector3l>> resultatsActuels;

    for (unsigned int i=0; t<tmax; i++)
    {
        if (i%skip==0) // save this point
        {
            resultatsActuels = {};
            for (unsigned int j=0; j<=body.size(); j++) {
                resultatsActuels.push_back( pair<vector3l,vector3l>( body[j].getPosition(), body[j].getVelocity() ));
            }
            trajectories.push_back( pair<long double,vector<pair<vector3l,vector3l>>>(t,resultatsActuels)  );
        }
        long double currentStep = getStep();
        if (currentStep==0) {cout << "Error : Delay dropped to 0" << endl; exit(1);}
        t += currentStep;
        tick();
        cout << "\r" << int(100*t/tmax) << " %";
    }
    cout << "\r";
    computationDone=true;
}

long double NBody::getRecordedStep(int i)
{
    if (recordedTrajectoryIndex+i>=trajectories.size() || recordedTrajectoryIndex+i<0) return -1;
    return trajectories[recordedTrajectoryIndex+i].first-trajectories[recordedTrajectoryIndex].first; // Delai avant la prochaine frame
}

vector<pair<vector3l,vector3l>> NBody::getResults()
{
    return trajectories[recordedTrajectoryIndex].second;
}

unsigned int NBody::factorial(unsigned int n)
{
    if (n<=maxn)
    {
        return f[n];
    }
    return n*factorial(n-1);
}

unsigned int NBody::binomialCoeff(unsigned int n, unsigned int k)
{
    if (k>n) {cout << "k>n, impossible de faire k parmis n" << endl; exit(1);}
    if (n>maxn) return factorial(n)/(factorial(k)*factorial(n-k));
    return C[n][k];
}

void NBody::tick(long double delta)
{
    if (computationMethod.integrator == Integrator::Euler || computationMethod.integrator==Integrator::Heun)
    {
        vector<vector3l> accelerations,accelerations2;
        bool heun = (computationMethod.integrator==Integrator::Heun);
        for (unsigned int i=0; i<body.size(); i++)
        {
            accelerations.push_back(body[i].getAcc());
            if (heun)
            {
                body[i].velocity += delta*accelerations[i];
                accelerations2.push_back(body[i].getAcc());
                body[i].velocity -= delta*accelerations[i];
            }
        }
        for (unsigned int i=0; i<body.size(); i++)
        {
            body[i].position += delta*body[i].velocity+ heun*(long double)0.5*(long double)pow((long double)delta,2)*accelerations[i];
            if (heun) body[i].velocity += (delta/2)*(accelerations[i]+accelerations2[i]);
            else body[i].velocity += delta*accelerations[i];
        }
    }

    if (computationMethod.integrator==Integrator::Taylor)
    {
        unsigned int order = computationMethod.order;
        if (order<1) {cout << "Error : Taylor cannot be less than order 1 for working computations" << endl; exit(EXIT_FAILURE);}
        if (order==2) {computationMethod.integrator=Integrator::Euler; tick(delta);} // Order 1 Taylor is Euler

        // Initialisation des listes utilisées
        vector<vector<vector<long double>>> derivativesGamma(body.size()), derivativesDirSquare(body.size()); // i j n
        vector<vector<vector<vector3l>>> derivativesDir(body.size()), derivativesMu(body.size());
        vector<vector<vector3l>> derivativesPos(body.size());
        vector<vector<long double>> emptyScalar(body.size());
        vector<vector<vector3l>> emptyVector(body.size());
        for (unsigned int i=0; i<body.size(); i++) // fill lists by i j
        {
            derivativesDir[i] = emptyVector;
            derivativesGamma[i] = emptyScalar;
            derivativesDirSquare[i] = emptyScalar;
            derivativesMu[i] = emptyVector;
        }

        for (unsigned int n=0; n<=(unsigned)order-2; n++) // Derivatives computations
        {
            for (unsigned int i=0; i<body.size(); i++) // Computations gamma, dir square, mu
            {
                for (unsigned int j=0; j<body.size(); j++)
                {
                    if (i<j)
                    {
                        if (n==0) // Init
                        {
                            derivativesDir[i][j].push_back(body[j].getPosition()-body[i].getPosition());
                            derivativesDir[i][j].push_back(body[j].getVelocity()-body[i].getVelocity());
                            derivativesDirSquare[i][j].push_back(derivativesDir[i][j][0].square());
                            derivativesMu[i][j].push_back(derivativesDir[i][j][0]/(long double)pow(derivativesDirSquare[i][j][0],1.5));
                            derivativesGamma[i][j].push_back(-1/derivativesDir[i][j][0].norm());
                            continue;
                        }

                        long double newDirCarre=0;
                        for  (unsigned int k=0; k<=n; k++)
                        {
                            newDirCarre += (long double)binomialCoeff(n,k)*derivativesDir[i][j][k]*derivativesDir[i][j][n-k];
                        };
                        derivativesDirSquare[i][j].push_back(newDirCarre);

                        long double newGamma=0;
                        for  (unsigned int k=0; k<=n-1; k++)
                        {
                            newGamma += (long double)binomialCoeff(n-1,k)*derivativesDir[i][j][k+1]*derivativesMu[i][j][n-1-k];
                        }
                        derivativesGamma[i][j].push_back(newGamma);

                        vector3l newMu=VNULL;
                        for  (unsigned int k=0; k<=n-1; k++)
                        {
                            newMu += (long double)binomialCoeff(n,k)*( derivativesGamma[i][j][k]*derivativesDir[i][j][n-k] + derivativesDirSquare[i][j][n-k]*derivativesMu[i][j][k]);
                        }
                        newMu += derivativesGamma[i][j][n]*derivativesDir[i][j][0];
                        newMu *= -1/derivativesDirSquare[i][j][0];
                        derivativesMu[i][j].push_back(newMu);
                    }
                }
            }

            for (unsigned int i=0; i<body.size(); i++) // Position derivatives
            {
                if (n==0) // Init
                {
                    derivativesPos[i].push_back(body[i].getPosition());
                    derivativesPos[i].push_back(body[i].getVelocity());
                }
                vector3l newPos;
                for (unsigned int j=0; j<body.size(); j++)
                {
                    if (i<j) newPos += body[j].mass*derivativesMu[i][j][n];
                    else if (i>j) newPos -= body[j].mass*derivativesMu[j][i][n];
                }
                newPos *= Body::G;
                derivativesPos[i].push_back(newPos);
            }

            for (unsigned int i=0; i<body.size(); i++) // Derivatives of Delta position
            {
                for (unsigned int j=0; j<body.size(); j++)
                {
                    if (i<j) derivativesDir[i][j].push_back(derivativesPos[j][n+2]-derivativesPos[i][n+2]);
                }
            }
        }

        // Taylor series
        for (unsigned int i=0; i<body.size(); i++)
        {
            for (unsigned int n=1; n<=order-1; n++)
            {
                body[i].position += derivativesPos[i][n] * pow(delta,n)/(long double)factorial(n);
                body[i].velocity += derivativesPos[i][n+1] * pow(delta,n)/(long double)factorial(n);
                if (isinf( body[i].position.square())) exit(1);
            }
            body[i].position += derivativesPos[i][order] * pow(delta,order)/(long double)factorial(order);
        }
    }

    if (computationMethod.integrator == Integrator::RungeKutta)
    {
        unsigned int order = computationMethod.order;
        vector<vector<vector3l>> vitessesInter(order),accInter(order); // choose point then body
        vector3l incPos,incVit;
        for (unsigned int i=0; i<order; i++)
        {
            for (unsigned int k=0; k<body.size(); k++)
            {
                incPos = VNULL;
                incVit = VNULL;
                for (unsigned int j=0; j<i; j++)
                {
                    incPos += computationMethod.RKCoeffDerivatives[i][j]*vitessesInter[j][k];
                    incVit += computationMethod.RKCoeffDerivatives[i][j]*accInter[j][k];
                }
                body[k].position += incPos*delta;
                body[k].velocity += incVit*delta;
                vitessesInter[i].push_back(body[k].velocity);
                accInter[i].push_back(body[k].getAcc());
                body[k].position -= incPos*delta;
                body[k].velocity -= incVit*delta;
            }
        }
        for (unsigned int k=0; k<body.size(); k++) // Integral computation for each body
        {
            incPos = VNULL;
            incVit = VNULL;
            for (unsigned int i=0; i<order; i++) // for each point
            {
                incPos += computationMethod.RKCoeffAverage[i]*vitessesInter[i][k];
                incVit += computationMethod.RKCoeffAverage[i]*accInter[i][k];
            }
            body[k].position += incPos*delta;
            body[k].velocity += incVit*delta;
        }
    }

    if (computationMethod.integrator == Integrator::Symplectic4)
    {
        const long double a[] = {1/(2-pow(2,1./3.)),-pow(2,1./3.)/(2-pow(2,1./3.)),1/(2-pow(2,1./3.)),0}, \
            b[] = {1/(2*(2-pow(2,1./3.))), (1-pow(2,1./3.))/(2*(2-pow(2,1./3.))), (1-pow(2,1./3.))/(2*(2-pow(2,1./3.))), 1/(2*(2-pow(2,1./3.)))};
        vector<vector3l> accelerations;
        for (unsigned int j = 0; j<4; j++)
        {
            accelerations = {};
            for (unsigned int i=0; i<body.size(); i++)
            {
                accelerations.push_back(body[i].getAcc());
            }
            for (unsigned int i=0; i<body.size(); i++)
            {
                body[i].velocity += a[j]*delta*accelerations[i];
                body[i].position += b[j]*delta*body[i].velocity;
            }
        }
    }

    if (computationMethod.integrator == Integrator::Symplectic6)
    {
        const long double a[] = {0.0378593198406116, 0.102635633102435, -0.0258678882665587, 0.314241403071447, -0.130144459517415, 0.106417700369543, -0.00879424312851058, 0.207305069}, \
            b[] = {0.09171915262446165, 0.183983170005006, -0.05653436583288827, 0.004914688774712854, 0.143761127168358, 0.328567693746804, -0.196411466, 0};
        vector<vector3l> accelerations;
        long double stepb = delta/0.5, stepa = delta/0.603652534;
        for (unsigned int j = 0; j<8; j++)
        {
            accelerations = {};
            for (unsigned int i=0; i<body.size(); i++)
            {
                accelerations.push_back(body[i].getAcc());
            }
            for (unsigned int i=0; i<body.size(); i++)
            {
                body[i].velocity += a[j]*stepa*accelerations[i];
                body[i].position += b[j]*stepb*body[i].velocity;
            }
        }
    }


    if (computationMethod.energie)
    {
        long double energyDiff = getEnergie()-initEnergy;
        energyDiff = getEnergie()-initEnergy;
        vector<vector3l> forces,velocities;
        long double sForces(0),sVelocities(0);
        for (unsigned int i=0; i<body.size(); i++) // Computation of the required data
        {
            vector3l force = body[i].getForces();
            forces.push_back(force);
            sForces += force.square();
            velocities.push_back(body[i].getVelocity());
            sVelocities += body[i].getVelocity().square();
        }
        long double k = 1/(1+pow(positionPulseRatio,2)*(sVelocities/sForces));
        for (unsigned int i=0; i<body.size(); i++) // Energy correction
        {
            if (sForces!=0) body[i].position -= k*energyDiff*forces[i]/sForces;
            if (body[i].mass!=0) body[i].velocity += (1-k)*energyDiff*velocities[i]/(sForces*body[i].mass);
        }
    }
}

long double NBody::getStep()
{
    if (computationMethod.stepType==StepType::FixedStep) return normalStep;
    else return std::min( normalStep, stepFactor/pow(accelerationSum(),stepPower) );
}

bool NBody::getComputationDone()
{
    return computationDone;
}

void NBody::save(string filename)
{
    string resultsString;
    long double pmax(0),vmax(0);
    for (unsigned int i=0; i<trajectories.size(); i++)
    {
        resultsString.append(to_string(trajectories[i].first)); // time
        for (unsigned int j=0; j<body.size(); j++)
        {
            resultsString.append("," + to_string(trajectories[i].second[j].first*VX) + "," + to_string(trajectories[i].second[j].first*VY) + "," + to_string(trajectories[i].second[j].first*VZ)); // position
            resultsString.append("," + to_string(trajectories[i].second[j].second*VX) + "," + to_string(trajectories[i].second[j].second*VY) + "," + to_string(trajectories[i].second[j].second*VZ)); // vitesse
            pmax = max(pmax, trajectories[i].second[j].first.norm());
            vmax = max(vmax, trajectories[i].second[j].second.norm());
        }
        resultsString.append("\n");
    }

    ofstream fichier(filename);
    fichier << body.size() << "," << trajectories.size() << "," << pmax << "," << vmax << endl; // first line
    fichier << resultsString; // Enter the results

}

void NBody::tick()
{
    tick(getStep());
}

long double NBody::getEnergie()
{
    long double E=0;
    for (unsigned int i=0; i<body.size(); i++)
    {
        E+= body[i].getEnergy();
    }
    return E;
}

long double NBody::accelerationSum()
{
    long double sum = 0;
    for (unsigned int i=0; i<body.size(); i++)
    {
        sum += body[i].getAcc().norm();
    }
    return sum;
}

void NBody::linkAll()
{
    for (unsigned int i=0; i<this->body.size(); i++) {
        for (unsigned int j=0; j<this->body.size(); j++) {
            if (i!=j) this->body[i].attractors.push_back(&this->body[j]);
        }
    }
}
