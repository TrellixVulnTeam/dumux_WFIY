// $Id$

#ifndef DUNE_ONEPHASEMODEL_HH
#define DUNE_ONEPHASEMODEL_HH

#include <dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include "dumux/operators/p1operatorextended.hh"
#include "dumux/nonlinear/nonlinearmodel.hh"
#include "dumux/nonlinear/newtonmethod.hh"
#include "dumux/fvgeometry/fvelementgeometry.hh"


namespace Dune {

/** \todo Please doc me! */

template<class Grid, class Scalar, class ProblemType, class LocalJacobian,
		class FunctionType, class OperatorAssembler> class OnePhaseModel :
	public NonlinearModel<Grid, Scalar, ProblemType, LocalJacobian, FunctionType, OperatorAssembler> {
public:
	typedef NonlinearModel<Grid, Scalar, ProblemType, LocalJacobian,
	FunctionType, OperatorAssembler> ThisNonlinearModel;

	OnePhaseModel(const Grid& grid, ProblemType& prob) :
		ThisNonlinearModel(grid, prob), uOldTimeStep(grid, grid.overlapSize(0)==0) {
	}

	OnePhaseModel(const Grid& grid, ProblemType& prob, int level) :
		ThisNonlinearModel(grid, prob, level), uOldTimeStep(grid, level, grid.overlapSize(0)==0) {
	}

	virtual void initial() = 0;

	virtual void update(double& dt) = 0;

	virtual void solve() = 0;

	FunctionType uOldTimeStep;
};

/** \todo Please doc me! */

template<class Grid, class Scalar, class ProblemType, class LocalJac, int numEq=2>
class LeafP1OnePhaseModel
: public OnePhaseModel<Grid, Scalar, ProblemType, LocalJac,
		LeafP1Function<Grid, Scalar, numEq>, LeafP1OperatorAssembler<Grid, Scalar, numEq> >
{
public:
	// define the function type:
	typedef LeafP1Function<Grid, Scalar, numEq> FunctionType;

	// define the operator assembler type:
	typedef LeafP1OperatorAssembler<Grid, Scalar, numEq> OperatorAssembler;

	typedef OnePhaseModel<Grid, Scalar, ProblemType, LocalJac,
	FunctionType, OperatorAssembler> ThisOnePhaseModel;

	typedef LeafP1OnePhaseModel<Grid, Scalar, ProblemType, LocalJac, numEq> ThisType;

	typedef LocalJac LocalJacobian;

	// mapper: one data element per vertex
	template<int dim> struct P1Layout {
		bool contains(Dune::GeometryType gt) {
			return gt.dim() == 0;
		}
	};

	typedef typename Grid::LeafGridView GV;
    typedef typename GV::IndexSet IS;
	typedef MultipleCodimMultipleGeomTypeMapper<Grid,IS,P1Layout> VertexMapper;
	typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator
			IntersectionIterator;

	LeafP1OnePhaseModel(const Grid& grid, ProblemType& prob) :
		ThisOnePhaseModel(grid, prob), problem(prob), grid_(grid), vertexmapper(grid,
				grid.leafIndexSet()), size((*(this->u)).size()), p(size), x(size) { }

	virtual void initial() {
		typedef typename Grid::Traits::template Codim<0>::Entity Element;
		typedef typename Grid::ctype CoordScalar;
		typedef typename GV::template Codim<0>::Iterator Iterator;
		enum {dim = Grid::dimension};
		enum {dimworld = Grid::dimensionworld};

		const GV& gridview(this->grid_.leafView());

		// iterate through leaf grid an evaluate c0 at cell center
		Iterator eendit = gridview.template end<0>();
		for (Iterator it = gridview.template begin<0>(); it
				!= eendit; ++it) {
			// get geometry type
			Dune::GeometryType gt = it->geometry().type();

			// get element
			const Element& element = *it;

			const typename Dune::LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,dim>::value_type
					&sfs=Dune::LagrangeShapeFunctions<CoordScalar, Scalar, dim>::general(gt,
							1);
			int size = sfs.size();

			for (int i = 0; i < size; i++) {
				// get cell center in reference element
				const Dune::FieldVector<CoordScalar,dim>&local = sfs[i].position();

				// get global coordinate of cell center
				Dune::FieldVector<CoordScalar,dimworld> global = it->geometry().global(local);

				int globalId = vertexmapper.template map<dim>(element,
						sfs[i].entity());

				// initialize cell concentration
				(*(this->u))[globalId] = this->problem.initial(
						global, element, local);
			}
		}

		// set Dirichlet boundary conditions
		for (Iterator it = gridview.template begin<0>(); it
				!= eendit; ++it) {
			// get geometry type
			Dune::GeometryType gt = it->geometry().type();

			// get element
			const Element& element = *it;

			const typename Dune::LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,dim>::value_type
					&sfs=Dune::LagrangeShapeFunctions<CoordScalar, Scalar, dim>::general(gt,
							1);
			int size = sfs.size();

			// set type of boundary conditions
			this->localJacobian().template assembleBC<LeafTag>(element);

			IntersectionIterator
					endit = IntersectionIteratorGetter<Grid, LeafTag>::end(element);
			for (IntersectionIterator is = IntersectionIteratorGetter<Grid,
					LeafTag>::begin(element); is!=endit; ++is)
				if (is->boundary()) {
					for (int i = 0; i < size; i++)
						// handle subentities of this face
						for (int j = 0; j < ReferenceElements<CoordScalar,dim>::general(gt).size(is->numberInSelf(), 1, sfs[i].codim()); j++)
							if (sfs[i].entity()
									== ReferenceElements<CoordScalar,dim>::general(gt).subEntity(is->numberInSelf(), 1,
											j, sfs[i].codim())) {
								for (int equationNumber = 0; equationNumber<numEq; equationNumber++) {
									if (this->localJacobian().bc(i)[equationNumber]
											== BoundaryConditions::dirichlet) {
										// get cell center in reference element
										Dune::FieldVector<CoordScalar,dim>
												local = sfs[i].position();

										// get global coordinate of cell center
										Dune::FieldVector<CoordScalar,dimworld>
												global = it->geometry().global(local);

										int
												globalId = vertexmapper.template map<dim>(
														element, sfs[i].entity());
										FieldVector<int,numEq> dirichletIndex;
										FieldVector<BoundaryConditions::Flags, numEq>
												bctype = this->problem.bctype(
														global, element, is,
														local);
												this->problem.dirichletIndex(global, element, is,
														local, dirichletIndex);

										if (bctype[equationNumber]
												== BoundaryConditions::dirichlet) {
											FieldVector<Scalar,numEq>
													ghelp = this->problem.g(
															global, element, is,
															local);
											(*(this->u))[globalId][dirichletIndex[equationNumber]]
													= ghelp[dirichletIndex[equationNumber]];
										}
									}
								}
							}
				}
		}

		*(this->uOldTimeStep) = *(this->u);
		return;
	}

	virtual void update(double& dt) {
		this->localJacobian().setDt(dt);
		this->localJacobian().setOldSolution(this->uOldTimeStep);
		NewtonMethod<Grid, ThisType> newtonMethod(this->grid_, *this);
		newtonMethod.execute();
		dt = this->localJacobian().getDt();
		*(this->uOldTimeStep) = *(this->u);

		if (this->problem.exsolution)
			this->problem.updateExSol(dt, *(this->u));

		return;
	}

	virtual void globalDefect(FunctionType& defectGlobal) {
		typedef typename Grid::Traits::template Codim<0>::Entity Element;
		typedef typename Grid::ctype CoordScalar;
		typedef typename GV::template Codim<0>::Iterator Iterator;
		enum {dim = Grid::dimension};
		typedef array<BoundaryConditions::Flags, numEq> BCBlockType;

		const GV& gridview(this->grid_.leafView());
		(*defectGlobal)=0;

		// allocate flag vector to hold flags for essential boundary conditions
		std::vector<BCBlockType> essential(this->vertexmapper.size());
		for (typename std::vector<BCBlockType>::size_type i=0; i
				<essential.size(); i++)
			essential[i].assign(BoundaryConditions::neumann);

		// iterate through leaf grid
		Iterator eendit = gridview.template end<0>();
		for (Iterator it = gridview.template begin<0>(); it
				!= eendit; ++it) {
			// get geometry type
			Dune::GeometryType gt = it->geometry().type();

			// get element
			const Element& element = *it;
			this->localJacobian().fvGeom.update(element);
			int size = this->localJacobian().fvGeom.numVertices;

			this->localJacobian().setLocalSolution(element);
			this->localJacobian().computeElementData(element);
			bool old = true;
			this->localJacobian().updateVariableData(element, this->localJacobian().uold, old);
			this->localJacobian().updateVariableData(element, this->localJacobian().u);
			this->localJacobian().template localDefect<LeafTag>(element, this->localJacobian().u);

			// begin loop over vertices
			for (int i=0; i < size; i++) {
				int globalId = this->vertexmapper.template map<dim>(element,i);
				for (int equationnumber = 0; equationnumber < numEq; equationnumber++) {
					if (this->localJacobian().bc(i)[equationnumber] == BoundaryConditions::neumann)
						(*defectGlobal)[globalId][equationnumber]
								+= this->localJacobian().def[i][equationnumber];
					else
						essential[globalId].assign(BoundaryConditions::dirichlet);
				}
			}
		}

		for (typename std::vector<BCBlockType>::size_type i=0; i
				<essential.size(); i++)
			for (int equationnumber = 0; equationnumber < numEq; equationnumber++) {
			if (essential[i][equationnumber] == BoundaryConditions::dirichlet)
				(*defectGlobal)[i][equationnumber] = 0;
			}
	}

//	virtual double injected(double& upperMass, double& oldUpperMass) {
//		typedef typename Grid::Traits::template Codim<0>::Entity Element;
//		typedef typename Grid::ctype CoordScalar;
//		typedef typename GV::template Codim<0>::Iterator Iterator;
//		enum {dim = Grid::dimension};
//		enum {dimworld = Grid::dimensionworld};
//
//		const GV& gridview(this->grid_.leafView());
//		double totalMass = 0;
//		upperMass = 0;
//		oldUpperMass = 0;
//		// iterate through leaf grid an evaluate c0 at cell center
//		Iterator eendit = gridview.template end<0>();
//		for (Iterator it = gridview.template begin<0>(); it
//				!= eendit; ++it) {
//			// get geometry type
//			Dune::GeometryType gt = it->geometry().type();
//
//			// get element
//			const Element& element = *it;
//
//			FVElementGeometry<G> fvGeom;
//			fvGeom.update(element);
//
//			const typename Dune::LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,dim>::value_type
//					&sfs=Dune::LagrangeShapeFunctions<CoordScalar, Scalar, dim>::general(gt,
//							1);
//			int size = sfs.size();
//
//			for (int i = 0; i < size; i++) {
//				// get cell center in reference element
//				const Dune::FieldVector<CoordScalar,dim>&local = sfs[i].position();
//
//				// get global coordinate of cell center
//				Dune::FieldVector<CoordScalar,dimworld> global = it->geometry().global(local);
//
//				int globalId = vertexmapper.template map<dim>(element,
//						sfs[i].entity());
//
//				double volume = fvGeom.subContVol[i].volume;
//
//				double porosity = this->problem.soil().porosity(global, element, local);
//
//				double density = this->problem.materialLaw().nonwettingPhase.density();
//
//				double mass = volume*porosity*density*((*(this->u))[globalId][1]);
//
//				totalMass += mass;
//
//				if (global[2] > 80.0) {
//					upperMass += mass;
//					oldUpperMass += volume*porosity*density*((*(this->uOldTimeStep))[globalId][1]);
//				}
//			}
//		}
//
//		return totalMass;
//	}

	virtual void vtkout(const char* name, int k) {
		VTKWriter<typename Grid::LeafGridView> vtkwriter(this->grid_.leafView());
		char fname[128];
		sprintf(fname, "%s-%05d", name, k);
//		double minSat = 1e100;
//		double maxSat = -1e100;
//		if (problem.exsolution) {
//			satEx.resize(size);
//			satError.resize(size);
//		}
		for (int i = 0; i < size; i++) {
			p[i] = (*(this->u))[i][0];
			x[i] = (*(this->u))[i][1];
//			satW[i] = 1 - satN[i];
//			double satNI = satN[i];
//			minSat = std::min(minSat, satNI);
//			maxSat = std::max(maxSat, satNI);
//			pN[i] = (*(this->u))[i][1];
//			pC[i] = pN[i] - pW[i];
//			satW[i] = this->problem.materialLaw().saturationW(pC[i]);
//			satN[i] = 1 - satW[i];
//			if (problem.exsolution) {
//				satEx[i]=problem.uExOutVertex(i, 1);
//				satError[i]=problem.uExOutVertex(i, 2);
//			}
		}
		vtkwriter.addVertexData(p, "pressure");
		vtkwriter.addVertexData(x, "mole fraction");
//		if (problem.exsolution) {
//			vtkwriter.addVertexData(satEx, "saturation, exact solution");
//			vtkwriter.addVertexData(satError, "saturation error");
//		}
		vtkwriter.write(fname, VTKOptions::ascii);
//		std::cout << "nonwetting phase saturation: min = "<< minSat
//				<< ", max = "<< maxSat << std::endl;
//		if (minSat< -0.5 || maxSat > 1.5)DUNE_THROW(MathError, "Saturation exceeds range.");
	}

    const Grid &grid() const
        { return grid_; }

protected:
	ProblemType& problem;
	const Grid& grid_;
	VertexMapper vertexmapper;
	int size;
	BlockVector<FieldVector<Scalar, 1> > p;
	BlockVector<FieldVector<Scalar, 1> > x;
//	BlockVector<FieldVector<Scalar, 1> > pC;
//	BlockVector<FieldVector<Scalar, 1> > satW;
//	BlockVector<FieldVector<Scalar, 1> > satN;
//	BlockVector<FieldVector<Scalar, 1> > satEx;
//	BlockVector<FieldVector<Scalar, 1> > pEx;
//	BlockVector<FieldVector<Scalar, 1> > satError;
};

}
#endif
