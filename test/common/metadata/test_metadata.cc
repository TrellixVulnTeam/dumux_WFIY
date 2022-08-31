// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief This file tests the metadata functionality
 */
#include <config.h>

#include <dune/common/float_cmp.hh>

#include <dumux/common/metadata.hh>

struct MyDataClass
{
    double val;
    std::string name;
};

int main()
{
    using namespace Dumux::MetaData;

    Collector collector;

    // Test brackets operator and Json tree
    collector["MyData"]["values"] = {1,2};
    collector["MyData"]["myBool"] = true;

    if(!collector["MyData"]["myBool"])
        DUNE_THROW(Dune::Exception, "Wrong bool in JsonTree");

    auto& vals = collector["MyData"]["values"];
    vals.push_back(3);

    for(int i=0; i<vals.size(); i++)
        if (vals[i] != i+1)
            DUNE_THROW(Dune::Exception, "Wrong vector in JsonTree");

    // Test className function
    MyDataClass data({3.1,"mydata"});
    std::string className("MyDataClass");
    collector[className]["className"] = Collector::className(data, true);
    collector[className]["value"] = data.val;
    collector[className]["name"] = data.name;

    if (!(collector[className]["className"].get<std::string>().find(className) != std::string::npos))
        DUNE_THROW(Dune::Exception, "Class name is wrong");

    if (!Dune::FloatCmp::eq(3.1, collector[className]["value"].get<double>() , 1e-6))
        DUNE_THROW(Dune::Exception, "Wrong double value in JsonTree");

    std::string s(collector[className]["name"].get<std::string>());
    if(!(data.name.compare(s)==0))
        DUNE_THROW(Dune::Exception, "Wrong string in JsonTree");

    // Test writing and reading Json files
    writeJsonFile(collector, "output");
    Collector collector2;
    readJsonFile(collector2, "output");

    if(!(collector.getTree() == collector2.getTree()))
        DUNE_THROW(Dune::Exception, "Json trees differ");

    collector2[className]["value"] = 3.2;
    if((collector.getTree() == collector2.getTree()))
        DUNE_THROW(Dune::Exception, "Json trees do not differ");

    return 0;
}
