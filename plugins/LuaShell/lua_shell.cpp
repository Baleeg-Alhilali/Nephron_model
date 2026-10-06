/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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


#include "bridge/util.h"
#include "registry/registry.h"
#include "bindings/lua/lua_util.h"
#include "bindings/lua/lua_parsing.h"
#include "common/profiler/profiler.h"
#include "registry/class_name_provider.h"
#include "ug.h"
#include "lua_shell.h"
#include "common/util/file_util.h"

using namespace std;
using namespace ug::bridge;
using namespace ug::bridge::lua;
using namespace ug::script;

namespace ug{
namespace luashell{

/// error function to be used for lua_pcall
int LuaCallStackError( lua_State *L )
{
//	UG_LOG("LUA-ERROR! Call stack:\n");
//    ug::bridge::LuaStackTrace();
    return 1;
}

LuaShell::LuaShell()
{
	init_lua_state();
}

void LuaShell::reset()
{
	ReleaseDefaultLuaState();
	init_lua_state();
}

void LuaShell::parse_file(const char* filename)
{
	vector<char> fileBuf;
	if(ReadFile(filename, fileBuf, true)){
		if(!fileBuf.empty())
			run(&fileBuf.front());
	}
}

void LuaShell::run(const char* buffer)
{
	PROFILE_FUNC();
	lua_State* L = m_luaState;

	ClearAbortRunFlag();
	lua_pushcfunction(L, LuaCallStackError);

	int error = luaL_loadbuffer(L, buffer, strlen(buffer), "luashell-run-buffer");

	try{
		if(error == 0){
			error = lua_pcall(L, 0, 0, -2);
		}
	}
	catch(SoftAbort& err){
		UG_LOG("Execution of LuaShell::run aborted with the following message:\n")
		UG_LOG(err.get_msg() << std::endl);
		reset();
	}

	if(error){
		string msg = lua_tostring(L, -1);
		lua_pop(L, 1);
		if(msg.find("__UG__LUA__EMPTY__MSG__") == string::npos)
			throw(LuaError(msg.c_str()));
		else
			throw(LuaError());
	}
	lua_pop(m_luaState, 1);
}

void LuaShell::abort_run(const char* message)
{
	AbortRun();
}

void LuaShell::set(const char* name, void* pval, const char* className)
{
	LuaParsing<void*>::push(m_luaState, pval, className);
	lua_setglobal(m_luaState, name);
}

void LuaShell::set(const char* name, const void* pval, const char* className)
{
	LuaParsing<const void*>::push(m_luaState, pval, className);
	lua_setglobal(m_luaState, name);
}

void LuaShell::set(const char* name, SmartPtr<void> pval, const char* className)
{
	LuaParsing<SmartPtr<void> >::push(m_luaState, pval, className);
	lua_setglobal(m_luaState, name);
}

void LuaShell::set(const char* name, ConstSmartPtr<void> pval, const char* className)
{
	LuaParsing<ConstSmartPtr<void> >::push(m_luaState, pval, className);
	lua_setglobal(m_luaState, name);
}

void LuaShell::init_lua_state()
{
//	instead of the global lua state, one could also use a local one		
	m_luaState = GetDefaultLuaState();
//	replace LUAs print function with our own, to use UG_LOG
	lua_register(m_luaState, "print", UGLuaPrint );
	lua_register(m_luaState, "print_all", UGLuaPrintAllProcs );
	lua_register(m_luaState, "write", UGLuaWrite );

//	todo: allow for argv arguments
	char* argv = 0;
	SetLuaUGArgs(m_luaState, 0, &argv);
}

}// namespace LuaShell


extern "C" void
InitUGPlugin_LuaShell(Registry* reg, string grp)
{
	using namespace luashell;
	grp.append("LuaShell");

	RegisterDefaultLuaBridge(reg);
	
	typedef LuaShell	T;
	reg->add_class_<T>("LuaShell", grp)
		.add_constructor()
		.add_method("reset", &T::reset)
		.add_method("run", &T::run)
		.add_method("abort_run", &T::abort_run)
		.add_method("set", static_cast<void (T::*)(const char*, double)>(&T::set<double>))
		.add_method("set", static_cast<void (T::*)(const char*, bool)>(&T::set<bool>))
		.add_method("set", static_cast<void (T::*)(const char*, const char*)>(&T::set<const char*>))
		.add_method("set", static_cast<void (T::*)(const char*, std::string)>(&T::set<std::string>))
		.add_method("set", static_cast<void (T::*)(const char*, const std::string&)>(&T::set<const std::string&>))
		.add_method("set", static_cast<void (T::*)(const char*, void*, const char*)>(&T::set))
		.add_method("set", static_cast<void (T::*)(const char*, const void*, const char*)>(&T::set))
		.add_method("set", static_cast<void (T::*)(const char*, SmartPtr<void>, const char*)>(&T::set))
		.add_method("set", static_cast<void (T::*)(const char*, ConstSmartPtr<void>, const char*)>(&T::set))
		.add_method("get_number", &T::get_val<double>)
		.add_method("get_bool", &T::get_val<bool>)
		.add_method("get_string", &T::get_val<std::string>)
		.set_construct_as_smart_pointer(true);
}

}// namespace ug
