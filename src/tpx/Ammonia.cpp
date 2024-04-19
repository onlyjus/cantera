/**
 * @file Ammonia.cpp representation of substance Ammonia.
 *
 * Values and functions are from "Thermodynamic Properties in SI" by W.C.
 * Reynolds AUTHOR: me@rebeccahhunt.com: GCEP, Stanford University
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "Ammonia.h"
#include "cantera/base/stringUtils.h"

using namespace Cantera;

namespace tpx
{

/*
 * Ammonia constants
 */
static const double Tmn = 50; // [K] minimum temperature for which calculations are valid
static const double Tmx = 1500.0; // [K] maximum temperature for which calculations are valid
static const double Tc=406.80; // [K] critical temperature
static const double Roc=237.64; // [kg/m^3] critical density
static const double To=200; // [K] reference Temperature
static const double R=488.20981; // [] gas constant for NH3 J/kg/K
static const double Ta=500; // [K] ??
static const double tauc=1.2333498;
static const double Gamma=5.0E-6; // [??]
static const double u0=1.3814023E6; // [] internal energy at To
static const double s0=6.2092055E3; // [] entropy at To
static const double Tp=300; // [K] ??
static const double Pc=11.627E6; // [Pa] critical pressure
static const double M=17.03; // [kg/kmol] molar density

// array A is used by the function named Pp
static const double A[9][6]= {
    {-6.453022304053E-3,  -1.371992677050E-2,  -8.100620315713E-3,  -4.880096421085E-3,  -1.202877562682E-2,   6.806345929616E-3,  },
    { 8.080094367688E-6,   1.434692000561E-5,  -4.505297669943E-5,  -1.661889985705E-4,   3.790895022982E-5   -4.073020833373E-5,  },
    { 1.032994880724E-9,   5.584395580933E-8,   4.920166508177E-7,   1.737835999473E-6,  -3.087491526377E-8,   7.148353041627E-8,  },
    {-8.948264632008E-12, -1.697777441391E-10, -1.236532371672E-9,  -7.812161168317E-9,   1.779548269140E-12, -3.897461095850E-11, },
    {-6.692285882015E-14, -1.753943775320E-15,  2.085533713355E-13,  2.134894661440E-11,  0.0,                 0.0,                },
    { 2.473417459954E-16,  2.999839155475E-16,  4.509080578790E-15, -3.798084988179E-14,  0.0,                 0.0,                },
    {-3.065578854310E-19,  2.411655109855E-20, -9.323356799989E-18,  4.272409853059E-17,  0.0,                 0.0,                },
    { 1.617910033375E-22, -5.074780704643E-22,  8.139470397409E-21, -2.745871062656E-20,  0.0,                 0.0,                },
    {-2.782168879368E-26,  2.988129173133E-25, -2.772597352058E-24,  7.668928677925E-24,  0.0,                 0.0,                },
};

// array F is used by the function named Psat
static const double F[]= {
    -6.7232038,
    -1.4928492E-3,
    -2.1966350,
     1.8152441E-1,
     3.4255443E-1,
    -1.2772013E1,
    -5.8344087E1,
    -6.5163169E1,
};

// array D is used by the function ldens
static const double D[]= {
     2.3763863E2,
     2.2030340E2,
     1.1999997E3,
    -1.9145612E3,
     1.7358862E3,
    -5.5587491E2,
};

// array G is used by the function sp
static const double G[]= {
     1.469259288E3,
     2.411085448E-1,
    -7.038236532E-3,
     5.157906857E-5,
    -1.209815448E-7,
     1.440829341E-10,
    -9.429402197E-14,
     3.229595395E-17,
    -4.528318341E-21,
};

static const double taua[] = {1.544912, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5};


// equation Q-1 in Reynolds
double Ammonia::Q()
{
    double tau = Ta/T;
    double sum = 0.0;
    for (int i=1; i<=9; i++) {
        for (int j=1; j<=6; j++) {
            sum += A[i-1][j-1]*pow(Rho,i-1)*pow(tau-tauc,j-1);
        }
    }
    return sum;
}

double Ammonia::Qprime()
{
    double tau = Ta/T;
    double sum = 0.0;
    for (int i=1; i<=9; i++) {
        for (int j=1; j<=6; j++) {
            sum += A[i-1][j-1]*pow(Rho,i-1)*double(j-1)*pow(tau-tauc,j-2);
        }
    }
    return sum;
}

double Ammonia::C(int i)
{
    double tau = Ta/T;
    return (i == 0 ? R*T : R*T*(tau - tauc)*pow(tau - taua[i],i-1));
}

double Ammonia::Cprime(int i)
{
    double tau = Ta/T;
    return (i == 0 ? R : (i == 1 ? -R*tauc :
                          -R*pow(tau - taua[i],i-2)*(tauc*(tau - taua[i])
                                  + (i-1)*tau*(tau - tauc))));
}


// internal energy. See Reynolds eqn (15) section 2
double Ammonia::up()
{
    double sum = 0.0;
    int i;
    for (i=0; i<7; i++) {
        sum += (C(i) - T*Cprime(i));//*I(i);
    }
    // C-2 in Reynolds, integrated from To to T
    for (i=1; i<6; i++) {
        sum += G[i]*(pow(T,i) - pow(To,i))/double(i);
    }
    sum += G[0]*log(T/To) + u0;
    return sum + m_energy_offset;
}

// entropy. See Reynolds eqn (16) section 2
double Ammonia::sp()
{
    double sum = 0.0;
    int i;
    for (i=2; i<6; i++) {
        sum += G[i]*(pow(T,i-1) - pow(To,i-1))/double(i-1);
    }
    sum += G[1]*log(T/To);
    sum -= G[0]*(1.0/T - 1.0/To);
    sum += s0 - R*log(Rho);
    for (i=0; i<7; i++) {
        sum -= Cprime(i);//*I(i);
    }
    return sum + m_entropy_offset;
}

// Pressure. Equation P-3 in Reynolds.
double Ammonia::Pp()
{
    double P = Rho*R*T;

    P *= (1 + Rho*Q() + pow(Rho,2)*Qprime());
    return P;
}

// Saturation pressure. Equation S-2 in Reynolds.
double Ammonia::Psat()
{
    double log, sum=0,P;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("Ammonia::Psat",
                           "Temperature out of range. T = {}", T);
    }
    for (int i=1; i<=8; i++) {
        sum += F[i-1] * pow((T/Tp -1),double(i-1));
    }

    log = ((Tc/T)-1)*sum;
    P=exp(log)*Pc;
    return P;
}

// Liquid density. Equation D-2 in Reynolds.
double Ammonia::ldens()
{
    double xx=1-(T/Tc), sum=0;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("Ammonia::ldens",
                           "Temperature out of range. T = {}", T);
    }
    for (int i=1; i<=6; i++) {
        sum+=D[i-1]*pow(xx,double(i-1)/3.0);
    }
    return sum;
}

// The following functions allow users to get the properties of Ammonia
// that are not dependent on the state

double Ammonia::Tcrit()
{
    return Tc;
}
double Ammonia::Pcrit()
{
    return Pc;
}
double Ammonia::Vcrit()
{
    return 1.0/Roc;
}
double Ammonia::Tmin()
{
    return Tmn;
}
double Ammonia::Tmax()
{
    return Tmx;
}
double Ammonia::MolWt()
{
    return M;
}

}
