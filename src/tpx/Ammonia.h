//! @file Ammonia.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_AMMONIA_H
#define TPX_AMMONIA_H

#include "cantera/tpx/Sub.h"

namespace tpx
{

//! Pure species representation of carbon dioxide. Values and functions are
//! from Reynolds @cite reynolds1979.
class Ammonia : public Substance
{
public:
    Ammonia() {
        m_name="ammonia";
        m_formula="NH3";
    }

    double MolWt() override;
    double Tcrit() override;
    double Pcrit() override;
    double Vcrit() override;
    double Tmin() override;
    double Tmax() override;

    //! Pressure. Equation P-6 in Reynolds. P(rho, T).
    double Pp() override;

    /**
     * internal energy. See Reynolds eqn (15) section 2
     *
     *  u = (the integral from T to To of co(T)dT) +
     *         sum from i to N ([C(i) - T*Cprime(i)] + uo
     */
    double up() override;

    //! entropy. See Reynolds eqn (16) section 2
    double sp() override;

    //! Pressure at Saturation. Equation S-2 in Reynolds.
    double Psat() override;

private:
    //! Liquid density. Equation D2 in Reynolds.
    double ldens() override;

    //! Equation Q-1 in Reynolds.
    double Q();
    double Qprime();

    double C(int i);
    double Cprime(int i);

};

}

#endif // ! TPX_AMMONIA_H
