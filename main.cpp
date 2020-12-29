#include "src/NBody.h"
#include <string.h>
#include "src/rapidcsv.h"

using namespace std;

int main(int argc, char *argv[])
{
    string initFilename,outFilename="out.csv";
    ComputationMethod computationMethod = RungeKutta4;
    bool firstArgGiven = false;

    for (int iarg=1; iarg<argc; iarg++)
    {
        if (argv[iarg][0] == '-')
        {
            if (strcmp(argv[iarg],"-h")==0 || strcmp(argv[iarg],"--help")==0 )
            {
                cout << "Usage : nBodyProblem [options] initFile [outputFile]" << endl << endl \
                << "Options : " << endl << "\t -euler \t uses the explicit Euler integrator" << endl << "\t -taylor \t Uses a 4-order Taylor serie" << endl \
                << "\t -symplectic4 \t Uses a 4-order symplectic integrator" << endl << "\t -symplectic6 \t Uses a 6-order symplectic integrator" << endl \
                 << "\t -energy \t corrects the energy" << endl \
                 << endl << "The default integrator is the classic RK4" << endl;
                exit(EXIT_SUCCESS);
            }
            else if (strcmp(argv[iarg],"-euler")==0)
            {
                computationMethod.integrator = Integrator::Euler;
            }
            else if (strcmp(argv[iarg],"-taylor")==0)
            {
                computationMethod.integrator = Integrator::Taylor;
                computationMethod.order = 4;
            }
            else if (strcmp(argv[iarg],"-symplectic4")==0)
            {
                computationMethod.integrator = Integrator::Symplectic4;
            }
            else if (strcmp(argv[iarg],"-symplectic6")==0)
            {
                computationMethod.integrator = Integrator::Symplectic6;
            }
            else if (strcmp(argv[iarg],"-energy")==0)
            {
                computationMethod.energie = true;
            }
            else
            {
                cout << "Unkown option : '" << argv[iarg] <<"'" << endl << "Usage : nBodyProblem [options] initFile [outputFile]" << endl \
                    << "Use -h or --help for more informations" << endl;
                exit(EXIT_SUCCESS);
            }
        }
        else
        {
            if (!firstArgGiven)
            {
                initFilename = argv[iarg];
                firstArgGiven = true;
            }
            else
            {
                outFilename = argv[iarg];
            }
        }
    }
    if (!firstArgGiven)
    {
        cout << "Initialisation file not given" << endl << "Usage : nBodyProblem [options] initFile [outputFile]" << endl \
        << "Use -h or --help for more informations" << endl;
        exit(EXIT_SUCCESS);
    }


    NBody system;
    system.computationMethod = computationMethod;

    // Init file read
    rapidcsv::Document fichier(initFilename);
    vector<long double> line = fichier.GetRow<long double>(0); // First line
    unsigned int lineSize = line.size();
    if (lineSize!=3 && lineSize!=4 && lineSize!=5) {cout << "Error : the first line size is " << lineSize << ", should be between 3 and 5" << endl; exit(EXIT_FAILURE);}
    unsigned int n = line[0];
    system.normalStep = line[1];
    long double computationTime = line[2];
    unsigned int skip=100;
    if (lineSize>=4) skip = line[3];
    if (lineSize==5) Body::G = line[4];

    for (unsigned int i=0; i<n; i++)
    {
        line = fichier.GetRow<long double>(1+i);
        if (line.size()!=7) {cout << "Error : the size body number " << i << "'s line does not have the correct size (7 numbers expected)" << endl; exit(EXIT_FAILURE);}
        system.body.push_back(Body( vector3l(line[0],line[1],line[2]), vector3l(line[3],line[4],line[5]), line[6] ));
    }

    system.linkAll();
    system.Compute(computationTime,skip);
    cout << "Computation completed" << endl;

    system.save(outFilename);
    cout << "File saved" << endl;

    return EXIT_SUCCESS;
}
