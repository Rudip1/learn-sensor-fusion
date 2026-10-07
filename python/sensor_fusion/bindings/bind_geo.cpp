#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>

#include "bindings.hpp"
#include "sensor_fusion/geodesy.hpp"

namespace sensor_fusion::python {

namespace {

using geo::deg2rad;
using geo::rad2deg;
using Array = py::array_t<double, py::array::c_style | py::array::forcecast>;

// Inputs may be scalars or 1-D arrays; scalars give a single row back (shape (3,) instead of (1, 3)).
struct Batch {
    Eigen::VectorXd v;
    bool scalar;
};

Batch flat(const Array& a) {
    Batch b{Eigen::VectorXd(a.size()), a.ndim() == 0};
    const double* p = a.data();
    for (py::ssize_t i = 0; i < a.size(); ++i) b.v(i) = p[i];
    return b;
}

Eigen::VectorXd broadcast(const Batch& b, Eigen::Index n) {
    if (b.v.size() == n) return b.v;
    if (b.v.size() == 1) return Eigen::VectorXd::Constant(n, b.v(0));
    throw std::invalid_argument("inputs must have the same length or be scalars");
}

// (n, 3) or (3,) array -> n x 3 matrix
RowsX3 rows3(const Array& a) {
    if (a.size() % 3 != 0 || a.ndim() > 2 || (a.ndim() == 2 && a.shape(1) != 3))
        throw std::invalid_argument("expected an array of shape (n, 3) or (3,)");
    RowsX3 m(a.size() / 3, 3);
    std::copy(a.data(), a.data() + a.size(), m.data());
    return m;
}

py::object rows_out(const RowsX3& m, bool scalar) {
    if (scalar) return py::cast(Eigen::Vector3d(m.row(0).transpose()));
    return py::cast(m);
}

py::object vec_out(const Eigen::VectorXd& v, bool scalar) {
    if (scalar) return py::float_(v(0));
    return py::cast(v);
}

}  // namespace

void bind_geo(py::module_ m) {
    py::class_<geo::Ellipsoid>(m, "Ellipsoid")
        .def(py::init<double, double>(), py::arg("a"), py::arg("f"))
        .def_readonly("a", &geo::Ellipsoid::a)
        .def_readonly("f", &geo::Ellipsoid::f)
        .def_property_readonly("b", &geo::Ellipsoid::b)
        .def_property_readonly("e2", &geo::Ellipsoid::e2)
        .def_property_readonly("ep2", &geo::Ellipsoid::ep2);
    m.attr("WGS84") = geo::kWgs84;

    m.def(
        "prime_vertical_radius",
        [](const Array& lat_deg) {
            const Batch b = flat(lat_deg);
            Eigen::VectorXd out(b.v.size());
            for (Eigen::Index i = 0; i < out.size(); ++i) out(i) = geo::prime_vertical_radius(deg2rad(b.v(i)));
            return vec_out(out, b.scalar);
        },
        py::arg("lat_deg"), "N(phi), eq. (2.2)");
    m.def(
        "meridian_radius",
        [](const Array& lat_deg) {
            const Batch b = flat(lat_deg);
            Eigen::VectorXd out(b.v.size());
            for (Eigen::Index i = 0; i < out.size(); ++i) out(i) = geo::meridian_radius(deg2rad(b.v(i)));
            return vec_out(out, b.scalar);
        },
        py::arg("lat_deg"), "M(phi), eq. (2.6)");

    m.def(
        "geodetic_to_ecef",
        [](const Array& lat_deg, const Array& lon_deg, const Array& h) {
            const Batch la = flat(lat_deg), lo = flat(lon_deg), hh = flat(h);
            const Eigen::Index n = std::max({la.v.size(), lo.v.size(), hh.v.size()});
            const Eigen::VectorXd a = broadcast(la, n), b = broadcast(lo, n), c = broadcast(hh, n);
            RowsX3 out(n, 3);
            for (Eigen::Index i = 0; i < n; ++i)
                out.row(i) = geo::geodetic_to_ecef({deg2rad(a(i)), deg2rad(b(i)), c(i)}).transpose();
            return rows_out(out, la.scalar && lo.scalar && hh.scalar);
        },
        py::arg("lat_deg"), py::arg("lon_deg"), py::arg("h") = 0.0, "eq. (2.3); returns (n, 3) ECEF [m]");
    m.def(
        "ecef_to_geodetic",
        [](const Array& ecef) {
            const bool single = ecef.ndim() == 1;
            const auto pts = rows3(ecef);
            RowsX3 out(pts.rows(), 3);
            for (Eigen::Index i = 0; i < pts.rows(); ++i) {
                const geo::Geodetic g = geo::ecef_to_geodetic(pts.row(i).transpose());
                out.row(i) << rad2deg(g.lat), rad2deg(g.lon), g.h;
            }
            return rows_out(out, single);
        },
        py::arg("ecef"), "Algorithm 2.1; returns columns lat_deg, lon_deg, h");
    m.def(
        "rotation_ecef_enu",
        [](double lat_deg, double lon_deg) { return geo::rotation_ecef_enu(deg2rad(lat_deg), deg2rad(lon_deg)); },
        py::arg("lat_deg"), py::arg("lon_deg"), "R_EN: columns east, north, up in ECEF, eq. (2.4)");

    py::class_<geo::LocalTangentPlane>(m, "LocalTangentPlane", "local ENU frame at a geodetic origin, eq. (2.5)")
        .def(py::init([](double lat_deg, double lon_deg, double h) {
                 return geo::LocalTangentPlane({deg2rad(lat_deg), deg2rad(lon_deg), h});
             }),
             py::arg("lat_deg"), py::arg("lon_deg"), py::arg("h") = 0.0)
        .def(
            "geodetic_to_enu",
            [](const geo::LocalTangentPlane& self, const Array& lat_deg, const Array& lon_deg, const Array& h) {
                const Batch la = flat(lat_deg), lo = flat(lon_deg), hh = flat(h);
                const Eigen::Index n = std::max({la.v.size(), lo.v.size(), hh.v.size()});
                const Eigen::VectorXd a = broadcast(la, n), b = broadcast(lo, n), c = broadcast(hh, n);
                RowsX3 out(n, 3);
                for (Eigen::Index i = 0; i < n; ++i)
                    out.row(i) = self.geodetic_to_enu({deg2rad(a(i)), deg2rad(b(i)), c(i)}).transpose();
                return rows_out(out, la.scalar && lo.scalar && hh.scalar);
            },
            py::arg("lat_deg"), py::arg("lon_deg"), py::arg("h") = 0.0)
        .def(
            "enu_to_geodetic",
            [](const geo::LocalTangentPlane& self, const Array& enu) {
                const bool single = enu.ndim() == 1;
                const auto pts = rows3(enu);
                RowsX3 out(pts.rows(), 3);
                for (Eigen::Index i = 0; i < pts.rows(); ++i) {
                    const geo::Geodetic g = self.enu_to_geodetic(pts.row(i).transpose());
                    out.row(i) << rad2deg(g.lat), rad2deg(g.lon), g.h;
                }
                return rows_out(out, single);
            },
            py::arg("enu"))
        .def(
            "ecef_to_enu",
            [](const geo::LocalTangentPlane& self, const Array& ecef) {
                const bool single = ecef.ndim() == 1;
                const auto pts = rows3(ecef);
                RowsX3 out(pts.rows(), 3);
                for (Eigen::Index i = 0; i < pts.rows(); ++i) out.row(i) = self.ecef_to_enu(pts.row(i).transpose());
                return rows_out(out, single);
            },
            py::arg("ecef"))
        .def(
            "enu_to_ecef",
            [](const geo::LocalTangentPlane& self, const Array& enu) {
                const bool single = enu.ndim() == 1;
                const auto pts = rows3(enu);
                RowsX3 out(pts.rows(), 3);
                for (Eigen::Index i = 0; i < pts.rows(); ++i) out.row(i) = self.enu_to_ecef(pts.row(i).transpose());
                return rows_out(out, single);
            },
            py::arg("enu"))
        .def_property_readonly("rotation", &geo::LocalTangentPlane::rotation)
        .def_property_readonly("origin_ecef", &geo::LocalTangentPlane::origin_ecef);

    m.def(
        "geodetic_to_enu_flat",
        [](const Array& lat_deg, const Array& lon_deg, const Array& h, double lat0_deg, double lon0_deg, double h0) {
            const Batch la = flat(lat_deg), lo = flat(lon_deg), hh = flat(h);
            const Eigen::Index n = std::max({la.v.size(), lo.v.size(), hh.v.size()});
            const Eigen::VectorXd a = broadcast(la, n), b = broadcast(lo, n), c = broadcast(hh, n);
            const geo::Geodetic origin{deg2rad(lat0_deg), deg2rad(lon0_deg), h0};
            RowsX3 out(n, 3);
            for (Eigen::Index i = 0; i < n; ++i)
                out.row(i) = geo::geodetic_to_enu_flat({deg2rad(a(i)), deg2rad(b(i)), c(i)}, origin).transpose();
            return rows_out(out, la.scalar && lo.scalar && hh.scalar);
        },
        py::arg("lat_deg"), py::arg("lon_deg"), py::arg("h"), py::arg("lat0_deg"), py::arg("lon0_deg"),
        py::arg("h0") = 0.0, "first-order local coordinates, eq. (2.7)");

    m.def(
        "utm_zone", [](double lat_deg, double lon_deg) { return geo::utm_zone(deg2rad(lat_deg), deg2rad(lon_deg)); },
        py::arg("lat_deg"), py::arg("lon_deg"));
    m.def(
        "utm_central_meridian", [](int zone) { return rad2deg(geo::utm_central_meridian(zone)); }, py::arg("zone"),
        "degrees");
    m.def(
        "geodetic_to_utm",
        [](const Array& lat_deg, const Array& lon_deg, int zone) {
            const Batch la = flat(lat_deg), lo = flat(lon_deg);
            const Eigen::Index n = std::max(la.v.size(), lo.v.size());
            const Eigen::VectorXd a = broadcast(la, n), b = broadcast(lo, n);
            if (zone == 0 && n > 0) zone = geo::utm_zone(deg2rad(a(0)), deg2rad(b(0)));  // one zone per batch
            Eigen::VectorXd e(n), no(n), g(n), k(n);
            bool north = true;
            for (Eigen::Index i = 0; i < n; ++i) {
                const geo::Utm u = geo::geodetic_to_utm(deg2rad(a(i)), deg2rad(b(i)), zone);
                e(i) = u.easting;
                no(i) = u.northing;
                g(i) = rad2deg(u.convergence);
                k(i) = u.scale;
                north = u.north;
            }
            const bool s = la.scalar && lo.scalar;
            py::dict out;
            out["easting"] = vec_out(e, s);
            out["northing"] = vec_out(no, s);
            out["zone"] = zone;
            out["north"] = north;
            out["convergence_deg"] = vec_out(g, s);
            out["scale"] = vec_out(k, s);
            return out;
        },
        py::arg("lat_deg"), py::arg("lon_deg"), py::arg("zone") = 0,
        "eqs. (2.8)-(2.11); all points of a call are projected into one zone (the first point's by default)");
    m.def(
        "utm_to_geodetic",
        [](const Array& easting, const Array& northing, int zone, bool north) {
            const Batch eb = flat(easting), nb = flat(northing);
            const Eigen::Index n = std::max(eb.v.size(), nb.v.size());
            const Eigen::VectorXd e = broadcast(eb, n), no = broadcast(nb, n);
            Eigen::VectorXd lat(n), lon(n);
            for (Eigen::Index i = 0; i < n; ++i) {
                geo::Utm u;
                u.zone = zone;
                u.north = north;
                u.easting = e(i);
                u.northing = no(i);
                const geo::Geodetic g = geo::utm_to_geodetic(u);
                lat(i) = rad2deg(g.lat);
                lon(i) = rad2deg(g.lon);
            }
            const bool s = eb.scalar && nb.scalar;
            return py::make_tuple(vec_out(lat, s), vec_out(lon, s));
        },
        py::arg("easting"), py::arg("northing"), py::arg("zone"), py::arg("north") = true,
        "inverse UTM, eq. (2.12); returns (lat_deg, lon_deg)");
}

}  // namespace sensor_fusion::python
