

#include "super_lu_bridge.h"


#ifdef UG_USE_PYBIND11



PYBIND11_MODULE(pysuperlu, m)
{
	m.doc() = "SuperLU module";
	m.attr("__name__") = "ug4py.superlu";

	ug::pybind::Registry registry(m);
	std::string name("SuperLU");

	ug::bridge::SuperLUBridge::InitUGPlugin(&registry, name);
}
#endif
