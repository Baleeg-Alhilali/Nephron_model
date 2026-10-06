#pragma once

#include "bridge/util.h"

extern "C" void InitUGPlugin_SuperLU6(ug::bridge::Registry* reg, std::string grp);

#ifdef UG_USE_PYBIND11

#include "bindings/pybind/ug_pybind.h"

namespace ug {
namespace bridge {
namespace SuperLUBridge{
	void InitUGPlugin(ug::pybind::Registry* reg, std::string grp);
}
}
}
#endif
