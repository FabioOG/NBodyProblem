#ifndef VECTEUR_H_INCLUDED
#define VECTEUR_H_INCLUDED

#include <math.h>
#include <iostream>

#define VNULL vector3l(0,0,0)
#define VX vector3l(1,0,0)
#define VY vector3l(0,1,0)
#define VZ vector3l(0,0,1)

template <typename K> class vector3;

template <typename K> vector3<K> operator+(vector3<K> u, vector3<K> v);
template <typename K> vector3<K> operator-(vector3<K> u, vector3<K> v);
template <typename K> K operator*(vector3<K> u, vector3<K> v);
template <typename K> vector3<K> operator*(K a, vector3<K> v);
template <typename K> vector3<K> operator*(vector3<K> v, K a);
template <typename K> vector3<K> operator/(vector3<K> v, K a);
template <typename K> vector3<K> operator/(vector3<K> v, vector3<K> theta);
template <typename K> vector3<K> operator%(vector3<K> u, vector3<K> v);
template <typename K> bool operator==(vector3<K> u, vector3<K> v);
template <typename K> std::ostream &operator<<(std::ostream &flux, vector3<K> u);
template <typename K> std::istream &operator>>(std::istream &flux, vector3<K>& u);

template <typename K>
class vector3
{

    protected :

        float m_x,m_y,m_z;


    public :

        vector3(K x=0, K y=0, K z=0);
        vector3(const vector3<K>& u);

        void setSpherical(K m=0, K a_x_y=0, K a_xy_z=0);

        K norm() const;
        K square() const;
        vector3 normaliser() const;

        vector3 operator+=(vector3<K> u);
        vector3 operator-=(vector3<K> u);
        vector3 operator*=(K a);
        vector3 operator/=(K a);
        vector3 operator/=(vector3<K> theta);
        vector3 operator%=(vector3<K> u);

        friend vector3 operator+ <>(vector3<K> u, vector3<K> v);
        friend vector3 operator- <>(vector3<K> u, vector3<K> v);
        friend K operator* <>(vector3<K> u, vector3<K> v);
        friend vector3 operator* <>(K a, vector3<K> v);
        friend vector3 operator* <>(vector3<K> v, K a);
        friend vector3 operator/ <>(vector3<K> v, K a);
        friend vector3 operator/ <>(vector3<K> v, vector3<K> theta);
        friend vector3 operator% <>(vector3<K> u, vector3<K> v);
        friend bool operator== <>(vector3<K> u, vector3<K> v);
        friend std::ostream &operator<< <>(std::ostream &flux, vector3<K> u);
        friend std::istream &operator>> <>(std::istream &flux, vector3<K>& u);

};

typedef vector3<long double> vector3l;

template <typename K>
vector3<K>::vector3 (K x, K y, K z) {m_x=x; m_y=y; m_z=z;}
template <typename K>
vector3<K>::vector3 (const vector3<K>& u) {m_x=u.m_x; m_y=u.m_y; m_z=u.m_z;}

template <typename K>
void vector3<K>::setSpherical (K n, K a_x_y, K a_xy_z) {
    m_z = sin(a_xy_z)*n;
    K normxy = cos(a_xy_z)*n;
    m_x = cos(a_x_y)*normxy;
    m_y = sin(a_x_y)*normxy;
}

template <typename K>
K vector3<K>::norm() const {return sqrt((long double)square());}

template <typename K>
K vector3<K>::square() const {return *this * *this;}

template <typename K>
vector3<K> vector3<K>::normaliser() const {return *this/(this->getNorm());}

template <typename K>
vector3<K> operator+ (vector3<K> u, vector3<K> v) {
    return vector3<K>(u.m_x+v.m_x, u.m_y+v.m_y, u.m_z+v.m_z);
}

template <typename K>
vector3<K> operator- (vector3<K> u, vector3<K> v) {
    return vector3<K>(u.m_x-v.m_x, u.m_y-v.m_y, u.m_z-v.m_z);
}

template <typename K>
K operator* (vector3<K> u, vector3<K> v) {
    return u.m_x*v.m_x + u.m_y*v.m_y + u.m_z*v.m_z;
}

template <typename K>
vector3<K> operator* (K a, vector3<K> u) {
    return vector3<K>(a*u.m_x, a*u.m_y, a*u.m_z);
}

template <typename K>
vector3<K> operator* (vector3<K> u, K a) {
    return a*u;
}

template <typename K>
vector3<K> operator/ (vector3<K> u, K a) {
    return u*(1/a);
}

template <typename K>
vector3<K> operator/ (vector3<K> v, vector3<K> theta) { /** Rotation by the angle vector, Rodrigues formula **/
    return v*cos(theta.getNorm())+(1-cos(theta.getNorm()))*(v*(theta/theta.getNorm()))*(theta/theta.getNorm())+sin(theta.getNorm())*((theta/theta.getNorm())%v);
}

template <typename K>
vector3<K> operator% (vector3<K> u, vector3<K> v) { /** Vectorial product **/
    return vector3<K>(u.m_y*v.m_z-u.m_z*v.m_y, u.m_z*v.m_x-u.m_x*v.m_z, u.m_x*v.m_y-u.m_y*v.m_x);
}

template <typename K>
bool operator== (vector3<K> u, vector3<K> v) {
    return bool(u.m_x==v.m_x && u.m_y==v.m_y && u.m_z==v.m_z);
}

template <typename K>
bool operator!= (vector3<K> u, vector3<K> v) {
    return !(u==v);
}

template <typename K>
std::ostream &operator<< (std::ostream &flux, vector3<K> u) {
    return flux<<'('<<u.m_x<<';'<<u.m_y<<';'<<u.m_z<<')';
}

template <typename K>
std::istream &operator>> (std::istream &flux, vector3<K>& u) {
    char c1,c2,c3,c4;
    K x,y,z;
    flux>>c1>>x>>c2>>y>>c3>>z>>c4;
    if (c1=='(' && c2==';' && c3==';' && c4==')') {
        u.m_x = x;
        u.m_y = y;
        u.m_z = z;
    }
    return flux;
}

template <typename K>
vector3<K> vector3<K>::operator+=(vector3<K> u) {
    *this = *this + u;
    return *this;
}

template <typename K>
vector3<K> vector3<K>::operator-=(vector3<K> u) {
    *this = *this - u;
    return *this;
}

template <typename K>
vector3<K> vector3<K>::operator*=(K a) {
    *this = *this*a;
    return *this;
}

template <typename K>
vector3<K> vector3<K>::operator/=(K a) {
    *this = *this/a;
    return *this;
}

template <typename K>
vector3<K> vector3<K>::operator/=(vector3<K> theta) {
    *this = *this/theta;
    return *this;
}

template <typename K>
vector3<K> vector3<K>::operator%=(vector3<K> u) {
    *this = *this%u;
    return *this;
}

#endif // VECTEUR_H_INCLUDED
