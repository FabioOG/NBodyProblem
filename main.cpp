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

    system.computePreset(initFilename);
    cout << "Computation completed" << endl;

    system.save(outFilename);
    cout << "File saved" << endl;

    return EXIT_SUCCESS;
}
