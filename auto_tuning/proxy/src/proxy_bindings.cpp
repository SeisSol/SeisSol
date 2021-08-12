#include "proxy_common.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


PYBIND11_MODULE(seissol_proxy_bindings, module) {
  py::enum_<Kernel>(module, "Kernel")
      .value("all", Kernel::all)
      .value("local", Kernel::local)
      .value("neigh", Kernel::neigh)
      .value("ader", Kernel::ader)
      .value("localwoader", Kernel::localwoader)
      .value("neigh_dr", Kernel::neigh_dr)
      .value("godunov_dr", Kernel::godunov_dr)
      .export_values();

  py::class_<ProxyConfig>(module, "ProxyConfig")
      .def(py::init<>())
      .def_readwrite("cells", &ProxyConfig::cells)
      .def_readwrite("timesteps", &ProxyConfig::timesteps)
      .def_readwrite("kernel", &ProxyConfig::kernel)
      .def_readwrite("verbose", &ProxyConfig::verbose);

  py::class_<ProxyOutput>(module, "ProxyOutput")
      .def(py::init<>())
      .def_readwrite("time", &ProxyOutput::time)
      .def_readwrite("cycles", &ProxyOutput::cycles)
      .def_readwrite("libxsmm_num_total_gflop", &ProxyOutput::libxsmmNumTotalGFlop)
      .def_readwrite("pspamm_num_total_gflop", &ProxyOutput::pspammNumTotalGFlop)
      .def_readwrite("libxsmm_and_pspamm_num_total_gflop", &ProxyOutput::libxsmmAndpspammNumTotalGFlop)
      .def_readwrite("actual_non_zero_gflop", &ProxyOutput::actualNonZeroGFlop)
      .def_readwrite("actual_hardware_gflop", &ProxyOutput::actualHardwareGFlop)
      .def_readwrite("gib", &ProxyOutput::gib)
      .def_readwrite("non_zero_flop_per_cycle", &ProxyOutput::nonZeroFlopPerCycle)
      .def_readwrite("hardware_flop_per_cycle", &ProxyOutput::hardwareFlopPerCycle)
      .def_readwrite("bytes_per_cycle", &ProxyOutput::bytesPerCycle)
      .def_readwrite("non_zero_gflops", &ProxyOutput::nonZeroGFlops)
      .def_readwrite("hardware_gflops", &ProxyOutput::hardwareGFlops)
      .def_readwrite("gib_per_second", &ProxyOutput::gibPerSecond);

  py::class_<Aux>(module, "Aux")
      .def(py::init<>())
      .def("str_to_kernel", &Aux::str2kernel)
      .def("kernel_to_string", &Aux::kernel2str)
      .def("get_allowed_kernels", &Aux::getAllowedKernels)
      .def("display_output", &Aux::displayOutput);

  module.def("run_proxy", &runProxy, "runs seissol proxy");
}
