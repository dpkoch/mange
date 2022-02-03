#include <string>

#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include "mange/mange.hh"

namespace py = pybind11;

namespace mange::python {

template <class G>
auto bind_class(auto handle, const std::string& name) {
    using DomainType = typename G::DomainType;

    return py::class_<G>(handle, name.c_str())
        .def_readonly_static("DOF", &G::DOF)
        .def_readonly_static("DIM", &G::DIM)
        .def(py::init<>())
        .def_static("Identity", &G::Identity)
        .def_static("Random", &G::Random)
        .def_static("Exp", &G::Exp)
        .def_static("Jl", &G::Jl)
        .def_static("Jr", &G::Jr)
        .def_static("JlInverse", &G::JlInverse)
        .def_static("JrInverse", &G::JrInverse)
        .def_static("hat", &G::hat)
        .def_static("vee", &G::vee)
        .def("Log", &G::Log)
        .def("Ad", &G::Ad)
        .def("inverse", &G::inverse)
        .def(py::self * py::self)
        .def(py::self * DomainType())
        .def("setIdentity", &G::setIdentity)
        .def("isApprox", &G::isApprox)
        .def("isIdentity", &G::isIdentity)
        .def("matrix", &G::matrix);
}

PYBIND11_MODULE(mange_python, m) {
    m.doc() = "A Lie group library";

    bind_class<SO2>(m, "SO2").def("rotation", &SO2::rotation);

    bind_class<SE2>(m, "SE2")
        .def(py::init<SO2>())
        .def(py::init<SE2::DomainType>())
        .def(py::init<SO2, SE2::DomainType>());

    bind_class<SO3>(m, "SO3");

    bind_class<SE3>(m, "SE3")
        .def(py::init<SO3>())
        .def(py::init<SE3::DomainType>())
        .def(py::init<SO3, SE3::DomainType>());
}

}  // namespace mange::python