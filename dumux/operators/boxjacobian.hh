// $Id$

#ifndef DUNE_BOXJACOBIAN_HH
#define DUNE_BOXJACOBIAN_HH

#include<map>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<sstream>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>
#include <dune/grid/utility/intersectiongetter.hh>

#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>

#include"dumux/functions/p1functionextended.hh"
#include"localjacobian.hh"

/**
 * @file
 * @brief  compute local jacobian matrix for conforming finite elements for diffusion equation
 * @author Peter Bastian
 */

namespace Dune {
/** @addtogroup DISC_Disc
 *
 * @{
 */
/**
 * @brief compute local jacobian matrix for conforming finite elements for diffusion equation
 *
 */

//! A class for computing local jacobian matrices
/*! A class for computing local jacobian matrix for the
 diffusion equation

 div j = q; j = -K grad u; in Omega

 u = g on Gamma1; j*n = J on Gamma2.

 Uses conforming finite elements with the Lagrange shape functions.
 It should work for all dimensions and element types.
 All the numbering is with respect to the reference element and the
 Lagrange shape functions

 Template parameters are:

 - Grid  a DUNE grid type
 - RT    type used for return values
 */
template<class Imp, class G, class RT, int m,
	class BoxFunction = LeafP1FunctionExtended<G, RT, m> >
class BoxJacobian :
public LocalJacobian<Imp,G,RT,m> {
	// mapper: one data element per vertex
	template<int dim> struct P1Layout {
		bool contains(Dune::GeometryType gt) {
			return gt.dim() == 0;
		}
	};

	typedef typename G::ctype DT;
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename Entity::Geometry Geometry;
	typedef typename LocalJacobian<Imp,G,RT,m>::VBlockType VBlockType;
	typedef typename LocalJacobian<Imp,G,RT,m>::MBlockType MBlockType;
	typedef MultipleCodimMultipleGeomTypeMapper<G, typename G::Traits::LeafIndexSet, P1Layout>
	VertexMapper;

public:
	// define the number of components of your system, this is used outside
	// to allocate the correct size of (dense) blocks with a FieldMatrix
	enum {dim=G::dimension};
	enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,dim>::maxsize};

	//! Constructor
	BoxJacobian(bool levelBoundaryAsDirichlet_, const G& grid,
			BoxFunction& sol, bool procBoundaryAsDirichlet_=true) :
				vertexMapper(grid, grid.leafIndexSet()),
				levelBoundaryAsDirichlet(levelBoundaryAsDirichlet_),
				procBoundaryAsDirichlet(procBoundaryAsDirichlet_),
				currentSolution(sol), oldSolution(grid), dt(1) {
	}

	//**********************************************************
	//*																			*
	//*	Computation of the local defect								*
	//*																			*
	//**********************************************************

	template<class TypeTag>
	void localDefect(const Entity& e, const VBlockType* sol, bool withBC = true) {
    	for (int i=0; i < this->fvGeom.nNodes; i++) // begin loop over vertices / sub control volumes
			{
				// implicit Euler
				bool old = true;
				VBlockType massContrib = computeM(e, this->uold, i, old);
				massContrib *= -1.0;
				this->def[i] = massContrib;
			}

			//updateVariableData(e, sol);
			for (int i=0; i < this->fvGeom.nNodes; i++) // begin loop over vertices / sub control volumes
			{
				VBlockType massContrib = computeM(e, sol, i);
				this->def[i] += massContrib;
				this->def[i] *= this->fvGeom.subContVol[i].volume/dt;

				// get source term
				VBlockType q = computeQ(e, sol, i);
				q *= this->fvGeom.subContVol[i].volume;
				this->def[i] -= q;
			} // end loop over vertices / sub control volumes

			for (int k = 0; k < this->fvGeom.nEdges; k++) // begin loop over edges / sub control volume faces
			{
				int i = this->fvGeom.subContVolFace[k].i;
				int j = this->fvGeom.subContVolFace[k].j;

				VBlockType flux = computeA(e, sol, k);

				// add to defect
				this->def[i] -= flux;
				this->def[j] += flux;
//				std::cout << "i = " << i << ", j = " << j << ", flux = " << flux << std::endl;
			} // end loop over edges / sub control volume faces

			if (withBC) {
				// assemble boundary conditions
				assembleBC<TypeTag> (e);

				// add to defect
				for (int i=0; i < this->fvGeom.nNodes; i++) {
					this->def[i] += this->b[i];
				}
			}

			return;
	}

	void setLocalSolution(const Entity& e)
	{
		Dune::GeometryType gt = e.geometry().type();
		const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
		&sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt, 1);
		int size = sfs.size();
		this->setcurrentsize(size);

		for (int i = 0; i < size; i++)
			for (int comp = 0; comp < m; comp++) {
				this->u[i][comp] = currentSolution.evallocal(comp, e, sfs[i].position());
				this->uold[i][comp] = oldSolution.evallocal(comp, e, sfs[i].position());
			}

		return;
	}

	void localToGlobal(const Entity& e, const VBlockType* sol)
	{
		int size = e.template count<dim>();
		for (int i = 0; i < size; i++) {
			int globalIdx = vertexMapper.template map<dim>(e, i);
			(*currentSolution)[globalIdx] = this->u[i];
		}
	}


	void setDt(double d) {
		dt = d;

		return;
	}

	double getDt() {
		return dt;
	}

	void setOldSolution(BoxFunction& uOld) {
		*oldSolution = *uOld;
	}

	VBlockType computeM(const Entity& e, const VBlockType* sol, int node, bool old = false) {
		return this->getImp().computeM(e, sol, node, old);
	}

	VBlockType computeQ(const Entity& e, const VBlockType* sol, int node) {
		return this->getImp().computeQ(e, sol, node);
	}

	VBlockType computeA(const Entity& e, const VBlockType* sol, int face) {
		return this->getImp().computeA(e, sol, face);
	}

	// analog to EvalStaticData in MUFTE
	virtual void updateStaticData(const Entity& e, VBlockType* sol) {
		return this->getImp().updateStaticData(e, sol);
	}

	template<class TypeTag> void assembleBC(const Entity& e) {
		Dune::GeometryType gt = e.geometry().type();
		const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
		&sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt, 1);
		setcurrentsize(sfs.size());
		this->fvGeom.update(e);

		const typename ReferenceElementContainer<DT,dim>::value_type
		&referenceElement = ReferenceElements<DT, dim>::general(gt);

		for (int i = 0; i < sfs.size(); i++) {
			this->bctype[i].assign(BoundaryConditions::neumann);
			this->b[i] = 0;
			this->dirichletIndex[i] = 0;
		}

		// evaluate boundary conditions via intersection iterator
		typedef typename IntersectionIteratorGetter<G,TypeTag>::IntersectionIterator IntersectionIterator;

		IntersectionIterator endit = IntersectionIteratorGetter<G, TypeTag>::end(e);
		for (IntersectionIterator it = IntersectionIteratorGetter<G, TypeTag>::begin(e); it!=endit; ++it)
		{
			// if we have a neighbor then we assume there is no boundary (forget interior boundaries)
			// in level assemble treat non-level neighbors as boundary
			if (it->neighbor()) {
				if (levelBoundaryAsDirichlet && it->outside()->level()==e.level())
					continue;
				if (!levelBoundaryAsDirichlet)
					continue;
			}

			// determine boundary condition type for this face, initialize with processor boundary
			FieldVector<typename BoundaryConditions::Flags, m> bctypeface(BoundaryConditions::process);
			FieldVector<int,m> dirichletIdx(0);

			// handle face on exterior boundary, this assumes there are no interior boundaries
			if (it->boundary()) {
				int faceIdx = it->numberInSelf();
				// 				std::cout << "faceIdx = " << faceIdx << ", beginning: " << std::endl;
				// 				for (int i = 0; i < 4; i++)
				// 				  std::cout << "bctype[" << i << "] = " << this->bctype[i] << std::endl;

				int nNodesOfFace = referenceElement.size(faceIdx, 1, dim);
				for (int nodeInFace = 0; nodeInFace < nNodesOfFace; nodeInFace++) {
					int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);
					for (int equationNumber = 0; equationNumber < m; equationNumber++) {
						if (this->bctype[nodeInElement][equationNumber] == BoundaryConditions::neumann) {
							int bfIdx = this->fvGeom.boundaryFaceIndex(faceIdx,	nodeInFace);
							FieldVector<DT,dim> local = this->fvGeom.boundaryFace[bfIdx].ipLocal;
							FieldVector<DT,dim> global = this->fvGeom.boundaryFace[bfIdx].ipGlobal;
							bctypeface = this->getImp().problem.bctype(global, e, it, local); // eval bctype
							this->getImp().problem.dirichletIndex(global, e, it, local, dirichletIdx); // eval bctype
							//							 						std::cout << "faceIdx = " << faceIdx << ", nodeInElement = " << nodeInElement
							//							 							  << ", bfIdx = " << bfIdx << ", local = " << local << ", global = " << global
							//							 							  << ", bctypeface = " << bctypeface << std::endl;
							if (bctypeface[equationNumber]!=BoundaryConditions::neumann)
								break;
							VBlockType J = this->getImp().problem.J(global, e, it, local);
							J[equationNumber] *= this->fvGeom.boundaryFace[bfIdx].area;
							this->b[nodeInElement][equationNumber] += J[equationNumber];
						}
					}
				}

				bool nface(true); // check if face is a neumann face
				for(int i=0; i<m; i++)
				{
					if(bctypeface[i] != BoundaryConditions::neumann)
						nface = false; // was not a neumann face
				}
				if(nface == true)
					continue; // was a neumann face, go to next face
			}

			// If we are here, then it is
			// (i)   an exterior boundary face with Dirichlet condition, or
			// (ii)  a processor boundary (i.e. neither boundary() nor neighbor() was true), or
			// (iii) a level boundary in case of level-wise assemble
			// How processor boundaries are handled depends on the processor boundary mode

			bool pface(false);  // check if face is a process boundary
			for(int i=0; i<m; i++)
			{
				if (bctypeface[i]==BoundaryConditions::process
						&& procBoundaryAsDirichlet==false
						&& levelBoundaryAsDirichlet==false)
				{
					pface = true;
					break;
				}
			}
			if(pface == true)
				continue;   // if face is a process boundary it acts like homogeneous Neumann


			for (int equationNumber=0; equationNumber<m; equationNumber++) {
				for (int i=0; i<sfs.size(); i++) // loop over test function number
				{
					//this->dirichletIndex[i][equationNumber] = equationNumber;

					//std::cout<<"i = "<<i<<std::endl;
					if (sfs[i].codim()==0)
						continue; // skip interior dof
					if (sfs[i].codim()==1) // handle face dofs
					{
						if (sfs[i].entity()==it->numberInSelf()) {
							if (this->bctype[i][equationNumber] < bctypeface[equationNumber]) {
								this->bctype[i][equationNumber] = bctypeface[equationNumber];
								this->dirichletIndex[i][equationNumber] = dirichletIdx[equationNumber];

								if (bctypeface[equationNumber] == BoundaryConditions::process)
									this->b[i][equationNumber] = 0;
								if (bctypeface[equationNumber] == BoundaryConditions::dirichlet) {
									this->b[i][equationNumber] = 0;
								}
							}
						}
						continue;
					}
					// handle subentities of this face
					for (int j=0; j<ReferenceElements<DT,dim>::general(gt).size(it->numberInSelf(), 1, sfs[i].codim()); j++)
						if (sfs[i].entity()==ReferenceElements<DT,dim>::general(gt).subEntity(it->numberInSelf(), 1, j, sfs[i].codim()))
						{
							if (this->bctype[i][equationNumber] < bctypeface[equationNumber]) {
								this->bctype[i][equationNumber] = bctypeface[equationNumber];
								this->dirichletIndex[i][equationNumber] = dirichletIdx[equationNumber];
								if (bctypeface[equationNumber] == BoundaryConditions::process)
									this->b[i][equationNumber] = 0;
								if (bctypeface[equationNumber] == BoundaryConditions::dirichlet) {
									this->b[i][equationNumber] = 0;
								}
							}
						}
				}
			}
		}

	}

	// parameters given in constructor
	VertexMapper vertexMapper;
	bool levelBoundaryAsDirichlet;
	bool procBoundaryAsDirichlet;
	BoxFunction& currentSolution;
	BoxFunction oldSolution;

public:
	double dt;
};
}
#endif
