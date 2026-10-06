/*
 * Copyright (c) 2014:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG_lua_shell
#define __H__UG_lua_shell

#include "bindings/lua/lua_util.h"
#include "bindings/lua/lua_parsing.h"

namespace ug{
namespace luashell{

/** 
 *  \defgroup lua_shell Lua Shell
 *  \ingroup plugins
 *  This plugin provides the 'LuaShell' class through which one can execute
 *	Lua-Scripts with custom variables.
 */

/// \addtogroup lua_shell
/// \{
class LuaShell{
	public:
		LuaShell();

		void reset();

		void parse_file(const char* filename);
		void run(const char* buffer);

		void abort_run(const char* message);

		template <class TVal>
		void set(const char* name, TVal value)
		{
			bridge::lua::LuaParsing<TVal>::push(m_luaState, value);
			lua_setglobal(m_luaState, name);
		}

		void set(const char* name, void* pval, const char* className);
		void set(const char* name, const void* pval, const char* className);
		void set(const char* name, SmartPtr<void> pval, const char* className);
		void set(const char* name, ConstSmartPtr<void> pval, const char* className);

		template <class TVal>
		TVal get_val(const char* name)
		{
			lua_getglobal(m_luaState, name);
			if(!bridge::lua::LuaParsing<TVal>::check(m_luaState, -1)){
				lua_pop(m_luaState, 1);
				UG_THROW("LuaShell error: Couldn't convert " << name << " to requested type.");
			}

			TVal val = bridge::lua::LuaParsing<TVal>::get(m_luaState, -1);
			lua_pop(m_luaState, 1);
			return val;
		}

	private:
		lua_State* m_luaState;

		void init_lua_state();
};

/// \}
}// namespace LuaShell
}//	end of namespace


#endif	//__H__lua_shell
