#include "body.h"
#include "vector3.h"
#include <cmath>
#include <iostream>

long double Body::G=1;

using namespace std;

Body::Body(vector3l p, vector3l v, long double m)
{
    position = p;
    velocity = v;
    mass = (m>=0)*m - (m<0)*m;
}

vector3l Body::getPosition()
{
    return position;
}

vector3l Body::getVelocity()
{
    return velocity;
}

vector3l Body::getForces()
{
    return mass*getAcc();
}

vector3l Body::getAcc()
{
    vector3l acc = VNULL;
    for (unsigned int i=0; i<attractors.size(); i++)
    {
        vector3l dir = attractors[i]->position - position;
        if (dir.norm()!=0) acc += G * attractors[i]->mass * dir / (dir.square()*dir.norm());
    }
    return acc;
}

long double Body::getKineticEnergy()
{
    return 0.5*mass*velocity.square();
}

long double Body::getPotentialEnergy()
{
    long double E=0.f;
    for (unsigned int i=0; i<attractors.size(); i++)
    {
        vector3l dir = attractors[i]->position - position;
        E += attractors[i]->mass * G / dir.norm();
    }
    return E*mass;
}

long double Body::getEnergy()
{
    return getKineticEnergy()+getPotentialEnergy();
}

