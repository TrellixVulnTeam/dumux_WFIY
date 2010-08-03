// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief Generic class to read/write restart files.
 */
#ifndef DUMUX_RESTART_HH
#define DUMUX_RESTART_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <boost/format.hpp>

#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

namespace Dumux {
/*!
 * \brief Load or save a state of a model to/from the harddisk.
 */
class Restart {
    //! \brief Create a magic cookie for restart files, so that it is
    //!        unlikely to load a restart file for an incorrectly.
    template <class GridView>
    static const std::string magicRestartCookie_(const GridView &gridView)
    {
        const std::string gridName = "blubb"; //gridView.grid().name();
        const int dim = GridView::dimension;

        int numVertices = gridView.template size(dim);
        int numElements = gridView.template size(0);
        int numEdges = gridView.template size(dim-1);
        int numCPUs = gridView.comm().size();
        int rank = gridView.comm().rank();

        return
            (boost::format("DuMux restart file: "
                           "gridName='%s' "
                           "numCPUs=%d "
                           "myRank=%d "
                           "numElements=%d "
                           "numEdges=%d "
                           "numVertices=%d")
             %gridName
             %numCPUs
             %rank
             %numElements
             %numEdges
             %numVertices).str();
    }

    //! \brief Return the restart file name.
    template <class GridView>
    static const std::string restartFileName_(const GridView &gridView,
                                              const std::string &simName,
                                              double t)
    {
        int rank = gridView.comm().rank();
        return (boost::format("%s_time=%.3lf_rank=%05d.drs")
                %simName
                %t
                %rank).str();
    };


public:
    /*!
     * \brief Returns the name of the file which is (de-)serialized.
     */
    const std::string &fileName() const
    { return fileName_; }

    /*!
     * \brief Write the current state of the model to disk.
     */
    template <class Problem>
    void serializeBegin(Problem &problem)
    {
        const std::string magicCookie = magicRestartCookie_(problem.gridView());
        fileName_ = restartFileName_(problem.gridView(), 
                                     problem.name(),
                                     problem.timeManager().time());

        // open output file and write magic cookie
        outStream_.open(fileName_.c_str());
        outStream_.precision(20);

        serializeSectionBegin(magicCookie);
        serializeSectionEnd();
    }

    /*!
     * \brief The output stream to write the serialized data.
     */
    std::ostream &serializeStream()
    {
        return outStream_;
    }

    /*!
     * \brief Start a new section in the serialized output.
     */
    void serializeSectionBegin(const std::string &cookie)
    {
        outStream_ << cookie << "\n";
    }

    /*!
     * \brief End of a section in the serialized output.
     */
    void serializeSectionEnd()
    { outStream_ << "\n"; }

    /*!
     * \brief Serialize all leaf entities of a codim in a gridView.
     *
     * The actual work is done by Serializer::serialize(Entity)
     */
    template <int codim, class Serializer, class GridView>
    void serializeEntities(Serializer &serializer,
                           const GridView &gridView)
    {
        std::string cookie = (boost::format("Entities: Codim %d")%codim).str();
        serializeSectionBegin(cookie);

        // write element data
        typedef typename GridView::template Codim<codim>::Iterator Iterator;

        Iterator it = gridView.template begin<codim>();
        const Iterator &endIt = gridView.template end<codim>();
        for (; it != endIt; ++it) {
            serializer.serializeEntity(outStream_, *it);
            outStream_ << "\n";
        };

        serializeSectionEnd();
    }

    /*!
     * \brief Finish the restart file.
     */
    void serializeEnd()
    {
        outStream_.close();
    };


    /*!
     * \brief Start reading a restart file at a certain simulated
     *        time.
     */
    template <class Problem>
    void deserializeBegin(Problem &problem, double t)
    {
        fileName_ = restartFileName_(problem.gridView(),
                                     problem.name(),
                                     t);

        // open input file and read magic cookie
        inStream_.open(fileName_.c_str());
        if (!inStream_.good()) {
            DUNE_THROW(Dune::IOError,
                       "Restart file '"
                       << fileName_
                       << "' could not be opened properly");
        }

        // make sure that we don't open an empty file
        inStream_.seekg(0, std::ios::end);
        int pos = inStream_.tellg();
        if (pos == 0) {
            DUNE_THROW(Dune::IOError,
                       "Restart file '"
                       << fileName_
                       << "' is empty");
        }
        inStream_.seekg(0, std::ios::beg);

        const std::string magicCookie = 
            magicRestartCookie_(problem.gridView());

        deserializeSectionBegin(magicCookie);
        deserializeSectionEnd();
    }

    /*!
     * \brief The input stream to read the data which ought to be
     *        deserialized.
     */
    std::istream &deserializeStream()
    {
        return inStream_;
    }

    /*!
     * \brief Start reading a new section of the restart file.
     */
    void deserializeSectionBegin(const std::string &cookie)
    {
        if (!inStream_.good())
            DUNE_THROW(Dune::IOError,
                       "Encountered unexpected EOF in restart file.");
        std::string buf;
        std::getline(inStream_, buf);
        if (buf != cookie)
            DUNE_THROW(Dune::IOError,
                       "Could not start section '" << cookie << "'");
    };

    /*!
     * \brief End of a section in the serialized output.
     */
    void deserializeSectionEnd()
    {
        std::string dummy;
        std::getline(inStream_, dummy);
        for (int i = 0; i < dummy.length(); ++i) {
            if (!std::isspace(dummy[i])) {
                DUNE_THROW(Dune::InvalidStateException,
                           "Encountered unread values while deserializing");
            };
        }
    }


    /*!
     * \brief Deserialize all leaf entities of a codim in a grid.
     *
     * The actual work is done by Deserializer::deserialize(Entity)
     */
    template <int codim, class Deserializer, class GridView>
    void deserializeEntities(Deserializer &deserializer,
                             const GridView &gridView)
    {
        std::string cookie = (boost::format("Entities: Codim %d")%codim).str();
        deserializeSectionBegin(cookie);

        std::string curLine;

        // read entity data
        typedef typename GridView::template Codim<codim>::Iterator Iterator;
        Iterator it = gridView.template begin<codim>();
        const Iterator &endIt = gridView.template end<codim>();
        for (; it != endIt; ++it) {
            if (!inStream_.good()) {
                DUNE_THROW(Dune::IOError,
                           "Restart file is corrupted");
            }

            std::getline(inStream_, curLine);
            std::istringstream curLineStream(curLine);
            deserializer.deserializeEntity(curLineStream, *it);
        };

        deserializeSectionEnd();
    }

    /*!
     * \brief Stop reading the restart file.
     */
    void deserializeEnd()
    {
        inStream_.close();
    }

    /*!
     * \brief Returns the list of restart files in the current directory.
     */
    static void restartFileList(std::list<std::string> &fileList,
                                const std::string directory=".")
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Dumux::Restart::restartFileList()");
    };


private:
    std::string fileName_;
    std::ifstream inStream_;
    std::ofstream outStream_;
};
}

#endif
