// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: andreas.lauser _at_ iws.uni-stuttgart.de                         *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief This file tests the properties system.
 *
 * We define a few type tags and property tags, then we attach values
 * to (TypeTag, PropertyTag) tuples and finally we use them in the
 * main function and print some diagnostic messages.
 */
#include <dumux/common/propertysystem.hh>

#include <iostream>

namespace Dumux {
namespace Properties {

///////////////////
// Define some hierarchy of type tags:
//
//  CompactCar --- Sedan -_
//                         \.
//                          +- Pickup ---_
//                         /              \.
//  Truck ----------------^                \.
//  Tank ---------------------------------- HummerH1
///////////////////
NEW_TYPE_TAG(CompactCar);
NEW_TYPE_TAG(Truck);
NEW_TYPE_TAG(Tank);

NEW_TYPE_TAG(Sedan, INHERITS_FROM(CompactCar));
NEW_TYPE_TAG(Pickup, INHERITS_FROM(Sedan, Truck));

NEW_TYPE_TAG(HummerH1, INHERITS_FROM(Sedan, Pickup, Tank));

///////////////////
// Define the property tags:
// TopSpeed, NumSeats, CanonCaliber, GasUsage, AutomaticTransmission, Payload
///////////////////
NEW_PROP_TAG(TopSpeed); // [km/h]
NEW_PROP_TAG(NumSeats); // []
NEW_PROP_TAG(CanonCaliber); // [mm]
NEW_PROP_TAG(GasUsage); // [l/100km]
NEW_PROP_TAG(AutomaticTransmission); // true/false
NEW_PROP_TAG(Payload); // [t]

///////////////////
// Make the AutomaticTransmission default to false
SET_PROP_DEFAULT(AutomaticTransmission)
{ static const bool value = false; };

///////////////////
// Define some values for the properties on the type tags:
//
// (CompactCar, TopSpeed) = GasUsage*35
// (CompactCar, NumSeats) = 5
// (CompactCar, GasUsage) = 4
//
// (Truck, TopSpeed) = 100
// (Truck, NumSeats) = 2
// (Truck, GasUsage) = 12
// (Truck, Payload) = 35
//
// (Tank, TopSpeed) = 60
// (Tank, GasUsage) = 65
// (Tank, CanonCaliber) = 120
//
// (Sedan, GasUsage) = 7
// (Sedan, AutomaticTransmission) = true
//
// (Pickup, TopSpeed) = 120
// (Pickup, Payload) = 5
//
// (HummmerH1, TopSpeed) = (Pickup, TopSpeed)
///////////////////

SET_INT_PROP(CompactCar, TopSpeed, GET_PROP_VALUE(TypeTag, PTAG(GasUsage)) * 30);
SET_INT_PROP(CompactCar, NumSeats, 5);
SET_INT_PROP(CompactCar, GasUsage, 4);

SET_INT_PROP(Truck, TopSpeed, 100);
SET_INT_PROP(Truck, NumSeats, 2);
SET_INT_PROP(Truck, GasUsage, 12);
SET_INT_PROP(Truck, Payload, 35);

SET_INT_PROP(Tank, TopSpeed, 60);
SET_INT_PROP(Tank, GasUsage, 65);
SET_INT_PROP(Tank, CanonCaliber, 120);

SET_INT_PROP(Sedan, GasUsage, 7);
SET_BOOL_PROP(Sedan, AutomaticTransmission, true);

SET_INT_PROP(Pickup, TopSpeed, 120);
SET_INT_PROP(Pickup, Payload, 5);

SET_INT_PROP(HummerH1, TopSpeed, GET_PROP_VALUE(TTAG(Pickup), PTAG(TopSpeed)));

///////////////////
// Unmount the canon from the Hummer
UNSET_PROP(HummerH1, CanonCaliber);

}; // namespace Properties
}; // namespace Dumux


int main()
{
    // print all properties for all type tags
    std::cout << "---------------------------------------\n";
    std::cout << "-- Property values\n";
    std::cout << "---------------------------------------\n";

    std::cout << "---------- Values for CompactCar ----------\n";

    std::cout << "(CompactCar, TopSpeed) = " << GET_PROP_VALUE(TTAG(CompactCar), PTAG(TopSpeed)) << "\n";
    std::cout << "(CompactCar, NumSeats) = " << GET_PROP_VALUE(TTAG(CompactCar), PTAG(NumSeats)) << "\n";
    std::cout << "(CompactCar, GasUsage) = " << GET_PROP_VALUE(TTAG(CompactCar), PTAG(GasUsage)) << "\n";
    std::cout << "(CompactCar, AutomaticTransmission) = " << GET_PROP_VALUE(TTAG(CompactCar), PTAG(AutomaticTransmission)) << "\n";

    std::cout << "---------- Values for Truck ----------\n";

    std::cout << "(Truck, TopSpeed) = " << GET_PROP_VALUE(TTAG(Truck), PTAG(TopSpeed)) << "\n";
    std::cout << "(Truck, NumSeats) = " << GET_PROP_VALUE(TTAG(Truck), PTAG(NumSeats)) << "\n";
    std::cout << "(Truck, GasUsage) = " << GET_PROP_VALUE(TTAG(Truck), PTAG(GasUsage)) << "\n";
    std::cout << "(Truck, Payload) = " << GET_PROP_VALUE(TTAG(Truck), PTAG(Payload)) << "\n";
    std::cout << "(Truck, AutomaticTransmission) = " << GET_PROP_VALUE(TTAG(Truck), PTAG(AutomaticTransmission)) << "\n";

    std::cout << "---------- Values for Tank ----------\n";

    std::cout << "(Tank, TopSpeed) = " << GET_PROP_VALUE(TTAG(Tank), PTAG(TopSpeed)) << "\n";
    std::cout << "(Tank, GasUsage) = " << GET_PROP_VALUE(TTAG(Tank), PTAG(GasUsage)) << "\n";
    std::cout << "(Tank, AutomaticTransmission) = " << GET_PROP_VALUE(TTAG(Tank), PTAG(AutomaticTransmission)) << "\n";
    std::cout << "(Tank, CanonCaliber) = " << GET_PROP_VALUE(TTAG(Tank), PTAG(CanonCaliber)) << "\n";

    std::cout << "---------- Values for Sedan ----------\n";

    std::cout << "(Sedan, TopSpeed) = " << GET_PROP_VALUE(TTAG(Sedan), PTAG(TopSpeed)) << "\n";
    std::cout << "(Sedan, NumSeats) = " << GET_PROP_VALUE(TTAG(Sedan), PTAG(NumSeats)) << "\n";
    std::cout << "(Sedan, GasUsage) = " << GET_PROP_VALUE(TTAG(Sedan), PTAG(GasUsage)) << "\n";
    std::cout << "(Sedan, AutomaticTransmission) = " << GET_PROP_VALUE(TTAG(Sedan), PTAG(AutomaticTransmission)) << "\n";

    std::cout << "---------- Values for Pickup ----------\n";
    std::cout << "(Pickup, TopSpeed) = " << GET_PROP_VALUE(TTAG(Pickup), PTAG(TopSpeed)) << "\n";
    std::cout << "(Pickup, NumSeats) = " << GET_PROP_VALUE(TTAG(Pickup), PTAG(NumSeats)) << "\n";
    std::cout << "(Pickup, GasUsage) = " << GET_PROP_VALUE(TTAG(Pickup), PTAG(GasUsage)) << "\n";
    std::cout << "(Pickup, Payload) = " << GET_PROP_VALUE(TTAG(Pickup), PTAG(Payload)) << "\n";
    std::cout << "(Pickup, AutomaticTransmission) = " << GET_PROP_VALUE(TTAG(Pickup), PTAG(AutomaticTransmission)) << "\n";

    std::cout << "---------- Values for HummerH1 ----------\n";
    std::cout << "(HummerH1, TopSpeed) = " << GET_PROP_VALUE(TTAG(HummerH1), PTAG(TopSpeed)) << "\n";
    std::cout << "(HummerH1, NumSeats) = " << GET_PROP_VALUE(TTAG(HummerH1), PTAG(NumSeats)) << "\n";
    std::cout << "(HummerH1, GasUsage) = " << GET_PROP_VALUE(TTAG(HummerH1), PTAG(GasUsage)) << "\n";
    std::cout << "(HummerH1, Payload) = " << GET_PROP_VALUE(TTAG(HummerH1), PTAG(Payload)) << "\n";
    std::cout << "(HummerH1, AutomaticTransmission) = " << GET_PROP_VALUE(TTAG(HummerH1), PTAG(AutomaticTransmission)) << "\n";
    // CanonCaliber is explcitly unset for the Hummer -> this would not compile:
    // std::cout << "(HummerH1, CanonCaliber) = " << GET_PROP_VALUE(TTAG(HummerH1), PTAG(CanonCaliber)) << "\n";

    std::cout << "\n";
    std::cout << "---------------------------------------\n";
    std::cout << "-- Diagnostic messages\n";
    std::cout << "---------------------------------------\n";
    std::cout << "---- Message for (HummerH1, CanonCaliber) ---\n"
              << PROP_DIAGNOSTIC(TTAG(HummerH1), PTAG(CanonCaliber));
    std::cout << "---- Message for (HummerH1, GasUsage) ---\n"
              << PROP_DIAGNOSTIC(TTAG(HummerH1), PTAG(GasUsage));
    std::cout << "---- Message for (HummerH1, AutomaticTransmission) ---\n"
              << PROP_DIAGNOSTIC(TTAG(HummerH1), PTAG(AutomaticTransmission));

    return 0;
};
