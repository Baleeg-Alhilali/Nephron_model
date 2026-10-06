/*
 * Copyright (c) 2014:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/*
 * http://crd-legacy.lbl.gov/~xiaoye/SuperLU
 */

#ifndef __H__UG__SUPER_LU_SOLVER__
#define __H__UG__SUPER_LU_SOLVER__


#define SUPERLU_6_EXPERIMENTAL 1

#include "lib_algebra/operator/linear_solver/external_solvers/external_solvers.h"
#include <vector>


namespace ug{

struct SuperLUConfiguration
{
	bool bPrintStat;
	bool equil;
	enum
	{
	CPT_NATURAL, CPT_MMD_ATA, CPT_MMD_AT_PLUS_A, CPT_COLAMD
	} colPerm;
};


IExternalSolverImplementation *CreateSuperLUImplementation(SuperLUConfiguration &config);


template<typename TAlgebra>
class SuperLUSolver : public IExternalSolver<TAlgebra>
{
	SuperLUConfiguration config;
	IExternalSolverImplementation *impl;
public:

	using IExternalSolver<TAlgebra>::init;

	SuperLUSolver()
	{
		config.bPrintStat = false;
		config.equil = true;
		config.colPerm = SuperLUConfiguration::CPT_COLAMD;

		impl = CreateSuperLUImplementation(config);
	}

	void print_stat(bool b)
	{
		config.bPrintStat = b;
	}

	void equil(bool b)
	{
		config.equil = b;
	}

	void col_perm_natural()
	{
		config.colPerm = SuperLUConfiguration::CPT_NATURAL;
	}

	void col_perm_mdo_ATA()
	{
		config.colPerm = SuperLUConfiguration::CPT_MMD_ATA;
	}

	void col_perm_mdo_AT_plus_A()
	{
		config.colPerm = SuperLUConfiguration::CPT_MMD_AT_PLUS_A;
	}

	void col_perm_approximate()
	{
		config.colPerm = SuperLUConfiguration::CPT_COLAMD;
	}


	virtual bool double_apply(CPUAlgebra::vector_type &c, const CPUAlgebra::vector_type &d)
	{
		return impl->apply(c, d);
	}

	virtual const char *double_name() const
	{
		return impl->name();
	}

	virtual void double_init(const CPUAlgebra::matrix_type &mat)
	{
		impl->init(mat);
	}

};

} // end namespace ug

#endif

