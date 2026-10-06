/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
 * Authors: Martin Rupp, Christian Wehner
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


#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FVCR__PCRILUT__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FVCR__PCRILUT__

#include "lib_algebra/operator/interface/operator.h"
#include "lib_algebra/operator/interface/preconditioner.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{
	
/// PCRILUTPreconditioner class
/*	Threshold ILU method
*
*   parallelization equivalent to ILU parallelization (see ilu.h)
*
*	The matrix entries can be devided into four sectors: 
*	velocity-velocity, velocity-pressure, pressure-velocity, pressure-pressure
*	Threshold tolerance depends on sector of new entry
*	(given by parameters m_eps_vv, m_eps_vp, m_eps_pv, m_eps_pp)		
*	Type of entry a_ij is determined by diagonal entry of rows i and j
*   if a_ii=0 i is a pressure index, else a velocity index
*
*/
template <typename TAlgebra>
class PCRILUTPreconditioner : public IPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Matrix Operator type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	private:
		typedef typename matrix_type::value_type block_type;
		using IPreconditioner<TAlgebra>::debug_writer;
		using IPreconditioner<TAlgebra>::set_debug;

	public:
	//	Constructor
		PCRILUTPreconditioner(double eps=1e-6) :
			m_eps(eps), m_eps_vv(eps), m_eps_vp(eps), m_eps_pv(eps), m_eps_pp(eps)
		{
			m_info = false;
		};

		PCRILUTPreconditioner(number threshvv,number thresh_vp_pv_pp)
		{
			set_threshold(threshvv,thresh_vp_pv_pp);
			m_info = false;
		}

		PCRILUTPreconditioner(number threshvv,number threshvp,number threshpv,number threshpp)
		{
			set_threshold(threshvv,threshvp,threshpv,threshpp);
			m_info = false;
		};

		PCRILUTPreconditioner(double eps,bool info) :
		m_eps(eps), m_eps_vv(eps), m_eps_vp(eps), m_eps_pv(eps), m_eps_pp(eps), m_info(info)
		{
		}

		PCRILUTPreconditioner(number threshvv,number thresh_vp_pv_pp,bool info)
		{
			set_threshold(threshvv,thresh_vp_pv_pp);
			m_info = info;
		}

		PCRILUTPreconditioner(number threshvv,number threshvp,number threshpv,number threshpp,bool info)
		{
			set_threshold(threshvv,threshvp,threshpv,threshpp);
			m_info = info;
		}

	// 	Clone
	
	SmartPtr<ILinearIterator<vector_type> > clone()
	{
		SmartPtr<PCRILUTPreconditioner<algebra_type> > newInst(new PCRILUTPreconditioner<algebra_type>(m_eps));
		newInst->set_debug(debug_writer());
		newInst->set_damp(this->damping());
		newInst->set_threshold(this->m_eps_vv,this->m_eps_vp,this->m_eps_pv,this->m_eps_pp);
		newInst->set_info(this->m_info);
		return newInst;
	}

	// 	Destructor
		virtual ~PCRILUTPreconditioner()
		{
		};

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const {return true;}

	///	sets threshold for incomplete LU factorisation
	void set_threshold(number threshvv,number threshvp,number threshpv,number threshpp)
	{
		m_eps_vv = threshvv;
		m_eps_vp = threshvp;
		m_eps_pv = threshpv;
		m_eps_pp = threshpp;
		m_eps = std::min(std::min(m_eps_vv,m_eps_vp),std::min(m_eps_pv,m_eps_pp));
	}
	
	///	sets threshold for incomplete LU factorisation
	void set_threshold(number threshvv,number thresh_vp_pv_pp)
	{
		m_eps_vv = threshvv;
		m_eps_vp = thresh_vp_pv_pp;
		m_eps_pv = thresh_vp_pv_pp;
		m_eps_pp = thresh_vp_pv_pp;
		m_eps = std::min(m_eps_vv,m_eps_vp);
	}

	void set_threshold(number thresh)
	{
		m_eps_vv = thresh;
		m_eps_vp = thresh;
		m_eps_pv = thresh;
		m_eps_pp = thresh;
		m_eps = thresh;
	}
	
	///	sets storage information output to true or false
	void set_info(bool info)
	{
		m_info = info;
	}
		
	
protected:
	//	Name of preconditioner
	virtual const char* name() const {return "PCRILUTPreconditioner";}
	
	//	Preprocess routine
	virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
	{
		STATIC_ASSERT(matrix_type::rows_sorted, Matrix_has_to_have_sorted_rows);
		matrix_type& mat = *pOp;

#ifdef 	UG_PARALLEL
		//	copy original matrix
		// commented-out function seems to be broken (A.Vogel)
		//			MakeConsistent(mat, m_ILUT);

		// we use new instead (A.Vogel)
		m_A = mat;
		MatAddSlaveRowsToMasterRowOverlap0(m_A);

		//	set zero on slaves
		std::vector<IndexLayout::Element> vIndex;
		CollectUniqueElements(vIndex,  m_A.layouts()->slave());
		SetDirichletRow(m_A, vIndex);
#else
		//	copy original matrix
		m_A = mat;
#endif
		typedef typename matrix_type::connection connection;
		m_L.resize_and_clear(m_A.num_rows(), m_A.num_cols());
		m_U.resize_and_clear(m_A.num_rows(), m_A.num_cols());
		
		// con is the current line of L/U
		std::vector<typename matrix_type::connection> con;
		con.reserve(300);
		con.resize(0);
		
		static const size_t velocity=0;
		static const size_t pressure=1; 
		
		// init row 0 of U
		for(typename matrix_type::row_iterator i_it = m_A.begin_row(0); i_it != m_A.end_row(0); ++i_it)
			con.push_back(connection(i_it.index(), i_it.value()));
		m_U.set_matrix_row(0, &con[0], con.size());
		
		size_t totalentries=0;
		size_t maxentries=0;
		
		for(size_t i=1; i<m_A.num_rows(); i++)
		{				
			size_t itype;
			if (m_A(i,i)!=0) itype=velocity; else itype=pressure;
			
			con.resize(0);
			size_t u_part=0;
			
			// get the row A(i, .) into con
			double dmax=0;
			m_remove_zeros = false;
			if (m_remove_zeros==false){
				for(typename matrix_type::row_iterator i_it = m_A.begin_row(i); i_it != m_A.end_row(i); ++i_it)
				{
					con.push_back(connection(i_it.index(), i_it.value()));
					if(dmax < BlockNorm(i_it.value()))
						dmax = BlockNorm(i_it.value());
				}
			} else {
				typename matrix_type::row_iterator i_it = m_A.begin_row(i);
				con.push_back(connection(i_it.index(), i_it.value()));
				++i_it;
				for(; i_it != m_A.end_row(i); ++i_it)
				{
					if (i_it.value()==0) continue;
					con.push_back(connection(i_it.index(), i_it.value()));
					if(dmax < BlockNorm(i_it.value()))
						dmax = BlockNorm(i_it.value());
				}
			}
			// eliminate all entries A(i, k) with k<i with rows U(k, .) and k<i
			for(size_t i_it = 0; i_it < con.size(); ++i_it)
			{
				size_t k = con[i_it].iIndex;
				if(k >= i)
				{
					// safe where U begins / L ends in con
					u_part = i_it;
					break;
				}
				if(con[i_it].dValue == 0.0) continue;
				UG_ASSERT(m_U.num_connections(k) != 0 && m_U.begin_row(k).index() == k, "");
				block_type &ukk = m_U.begin_row(k).value();
				
				// add row k to row i by A(i, .) -= U(k,.)  A(i,k) / U(k,k)
				// so that A(i,k) is zero.
				// safe A(i,k)/U(k,k) in con, (later L(i,k) )
				block_type &d = con[i_it].dValue = con[i_it].dValue / ukk;
				
				typename matrix_type::row_iterator k_it = m_U.begin_row(k); // upper row iterator
				++k_it; // skip diag
				size_t j = i_it+1;
				while(k_it != m_U.end_row(k) && j < con.size())
				{
					
					// (since con and U[k] is sorted, we can do sth like a merge on the two lists)
					if(k_it.index() == con[j].iIndex)
					{
						// match
						con[j].dValue -= k_it.value() * d;
						++k_it;	++j;
					}
					else if(k_it.index() < con[j].iIndex)
					{
						// we have a value in U(k, (*k_it).iIndex), but not in A.
						// check tolerance criteria
						
						typename matrix_type::connection 
						c(k_it.index(), k_it.value() * d * -1.0);
						
						if(BlockNorm(c.dValue) > dmax * m_eps){
							
							size_t kitType;
							if (m_A(k_it.index(),k_it.index())!=0) kitType=velocity; else kitType=pressure;
						
						if ((itype==velocity)&&(kitType==velocity)){
							if(BlockNorm(c.dValue) > dmax * m_eps_vv)
							{
								// insert sorted
								con.insert(con.begin()+j, c);
								++j;
							}
						}
						if ((itype==velocity)&&(kitType==pressure)){
							if(BlockNorm(c.dValue) > dmax * m_eps_vp)
							{
								// insert sorted
								con.insert(con.begin()+j, c);
								++j;
							}
						}
						if ((itype==pressure)&&(kitType==velocity)){
							if(BlockNorm(c.dValue) > dmax * m_eps_pv)
							{
								// insert sorted
								con.insert(con.begin()+j, c);
								++j;
							}
						}
						if ((itype==pressure)&&(kitType==pressure)){
							if(BlockNorm(c.dValue) > dmax * m_eps_pp)
							{
								// insert sorted
								con.insert(con.begin()+j, c);
								++j;
							}
						} 
						}
						// else do some lumping
						++k_it;
					}
					else
					{
						// we have a value in A(k, con[j].iIndex), but not in U.
						++j;
					}
				}
				// insert new connections after last connection of row i
				if (k_it!=m_U.end_row(k)){
					for (;k_it!=m_U.end_row(k);++k_it){
						typename matrix_type::connection c(k_it.index(),-k_it.value()*d);
						if(BlockNorm(c.dValue) > dmax * m_eps){
							size_t kitType;
							if (m_A(k_it.index(),k_it.index())!=0) kitType=velocity; else kitType=pressure;
							if ((itype==velocity)&&(kitType==velocity)){
								if(BlockNorm(c.dValue) > dmax * m_eps_vv)
								{
									con.push_back(c);
								}
							} else {
								if ((itype==velocity)&&(kitType==pressure)){
									if(BlockNorm(c.dValue) > dmax * m_eps_vp)
									{
										con.push_back(c);
									}
								} else {
									if ((itype==pressure)&&(kitType==velocity)){
										if(BlockNorm(c.dValue) > dmax * m_eps_pv)
										{
											con.push_back(c);
										}
									} else {
								//		if ((itype==pressure)&&(kitType==pressure)){
											if(BlockNorm(c.dValue) > dmax * m_eps_pp)
											{
												con.push_back(c); 
											}
								//		} 
									}
								}
							}

						};
					}
				};
			}
			
			totalentries+=con.size();
			if(maxentries < con.size()) maxentries = con.size();
			
			// safe L and U
			m_L.set_matrix_row(i, &con[0], u_part);
			m_U.set_matrix_row(i, &con[u_part], con.size()-u_part);
			
		}
		
		m_L.defragment();
		m_U.defragment();
		m_A.defragment();
		
		if (m_info==true){
			UG_LOG("\n	PCRILUT storage information:\n");
			UG_LOG("	A nr of connections: " << m_A.total_num_connections()  << "\n");
			UG_LOG("	L+U nr of connections: " << m_L.total_num_connections()+m_U.total_num_connections() << "\n");
			UG_LOG("	Increase factor: " << (float)(m_L.total_num_connections() + m_U.total_num_connections() )/m_A.total_num_connections() << "\n");
		}
		
		return true;
	}
	
	//	Stepping routine
	virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
	{
#ifdef 	UG_PARALLEL
		SmartPtr<vector_type> spDtmp = d.clone();
		spDtmp->change_storage_type(PST_UNIQUE);
		// apply iterator: c = LU^{-1}*d (damp is not used)
		// L
		for(size_t i=0; i < m_L.num_rows(); i++)
		{
			// c[i] = d[i] - m_L[i]*c;
			c[i] = (*spDtmp)[i];
			for(typename matrix_type::row_iterator it = m_L.begin_row(i); it != m_L.end_row(i); ++it)
				MatMultAdd(c[i], 1.0, c[i], -1.0, it.value(), c[it.index()] );
			// lii = 1.0.
		}
		// U
		//
		// last row diagonal U entry might be close to zero with corresponding zero rhs 
		// when solving Navier Stokes system, therefore handle separately
		{
			size_t i=m_U.num_rows()-1;
			typename matrix_type::row_iterator it = m_U.begin_row(i);
			UG_ASSERT(it.index() == i, "");
			block_type &uii = it.value();
			typename vector_type::value_type s = c[i];
			// check if close to zero
			if (BlockNorm(uii)<m_small_lower){
				// set correction to zero
				c[i] = 0;
				if (BlockNorm(s)>m_small_upper){
					UG_LOG("Warning: zero entry in last row of U with corresponding non-zero rhs entry (" << BlockNorm(s) << ")\n");
				}
			} else {
				// c[i] = s/uii;
				InverseMatMult(c[i], 1.0, uii, s);
			}
		}
		// handle all other rows
		for(size_t i=m_U.num_rows()-2; ; i--)
		{
			typename matrix_type::row_iterator it = m_U.begin_row(i);
			UG_ASSERT(it.index() == i, "");
			block_type &uii = it.value();
			
			typename vector_type::value_type s = c[i];
			++it; // skip diag
			for(; it != m_U.end_row(i); ++it)
				// s -= it.value() * c[it.index()];
				MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()] );
			
			// c[i] = s/uii;
			InverseMatMult(c[i], 1.0, uii, s);
			
			if(i==0) break;
		}
		//	Correction is always consistent
		c.set_storage_type(PST_ADDITIVE);
		c.change_storage_type(PST_CONSISTENT);
#else
		// apply iterator: c = LU^{-1}*d (damp is not used)
		// L
		for(size_t i=0; i < m_L.num_rows(); i++)
		{
			// c[i] = d[i] - m_L[i]*c;
			c[i] = d[i];
			for(typename matrix_type::row_iterator it = m_L.begin_row(i); it != m_L.end_row(i); ++it)
				MatMultAdd(c[i], 1.0, c[i], -1.0, it.value(), c[it.index()] );
			// lii = 1.0.
		}
		// U
		//
		// last row diagonal U entry might be close to zero with corresponding zero rhs
		// when solving Navier Stokes system, therefore handle separately
		{
			size_t i=m_U.num_rows()-1;
			typename matrix_type::row_iterator it = m_U.begin_row(i);
			UG_ASSERT(it.index() == i, "");
			block_type &uii = it.value();
			typename vector_type::value_type s = c[i];
			// check if close to zero
			if (BlockNorm(uii)<m_small_lower){
				// set correction to zero
				c[i] = 0;
				if (BlockNorm(s)>m_small_upper){
					UG_LOG("Warning: zero entry in last row of U with corresponding non-zero rhs entry (" << BlockNorm(s) << ")\n");
				}
			} else {
				// c[i] = s/uii;
				InverseMatMult(c[i], 1.0, uii, s);
			}
		}
		// handle all other rows
		for(size_t i=m_U.num_rows()-2; ; i--)
		{
			typename matrix_type::row_iterator it = m_U.begin_row(i);
			UG_ASSERT(it.index() == i, "");
			block_type &uii = it.value();

			typename vector_type::value_type s = c[i];
			++it; // skip diag
			for(; it != m_U.end_row(i); ++it)
				// s -= it.value() * c[it.index()];
				MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()] );

			// c[i] = s/uii;
			InverseMatMult(c[i], 1.0, uii, s);

			if(i==0) break;
		}
#endif
		return true;
	}
	
	//	Postprocess routine
	virtual bool postprocess() {return true;}
	
protected:
	matrix_type m_L;
	matrix_type m_U;
	matrix_type m_A;
	number m_eps;
	number m_eps_vv;
	number m_eps_vp;
	number m_eps_pv;
	number m_eps_pp;
	bool m_info;
	bool m_remove_zeros;
	static const number m_small_lower;
	static const number m_small_upper;
};
	
template <typename TAlgebra>
const number PCRILUTPreconditioner<TAlgebra>::m_small_lower = 1e-9;
template <typename TAlgebra>
const number PCRILUTPreconditioner<TAlgebra>::m_small_upper = 1e-6;

} // end namespace ug

#endif
