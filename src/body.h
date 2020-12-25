#ifndef CORPS_H
#define CORPS_H

#include <vector>
#include "vector3.h"

class NBody;

class Body
{
    public:
        Body(vector3l p, vector3l v, long double m);
        vector3l getPosition();
        vector3l getVelocity();
        vector3l getForces();

        std::vector<Body*> attractors;

        vector3l getAcc();

        long double getKineticEnergy();
        long double getPotentialEnergy();
        long double getEnergy();

        friend class NBody;
        static long double G;

    protected:
        vector3l position,velocity;
        long double mass;

};

#endif // CORPS_H
