#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "bindings.hpp"
#include "sensor_fusion/cloud.hpp"

namespace sensor_fusion::python {

namespace {

using RowsXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <int Dim>
using RowsD = Eigen::Matrix<double, Eigen::Dynamic, Dim, Eigen::RowMajor>;

template <int Dim>
cloud::Points<Dim> cols_from_rows(const RowsXd& rows) {
    if (rows.cols() != Dim) throw std::invalid_argument("expected points of shape (n, " + std::to_string(Dim) + ")");
    return rows.transpose();
}

/// Queries may be one point (shape (D,)) or many (shape (m, D)).
RowsXd as_queries(const py::array_t<double, py::array::c_style | py::array::forcecast>& q, int dim, bool& single) {
    single = q.ndim() == 1;
    if (q.size() % dim != 0) throw std::invalid_argument("query dimension does not match the tree");
    RowsXd m(q.size() / dim, dim);
    std::copy(q.data(), q.data() + q.size(), m.data());
    return m;
}

template <int Dim>
void bind_tree(py::module_& m, const char* name) {
    using Tree = cloud::KdTree<Dim>;
    py::class_<Tree>(m, name, "static k-d tree (Algorithm 4.2); use cloud.KdTree(points) to build one")
        .def(py::init([](const RowsXd& pts, int leaf) { return Tree(cols_from_rows<Dim>(pts), leaf); }),
             py::arg("points"), py::arg("leaf_size") = 10)
        .def_property_readonly("size", &Tree::size)
        .def_property_readonly("depth", &Tree::depth)
        .def_property_readonly("points", [](const Tree& t) { return RowsD<Dim>(t.points().transpose()); })
        .def(
            "knn",
            [](const Tree& t, const py::array_t<double, py::array::c_style | py::array::forcecast>& q,
               int k) -> py::tuple {
                bool single = false;
                const RowsXd Q = as_queries(q, Dim, single);
                const int kk = std::min(k, t.size());
                Eigen::Matrix<long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> I(Q.rows(), kk);
                RowsXd D(Q.rows(), kk);
                std::vector<int> idx;
                std::vector<double> d2;
                for (Eigen::Index i = 0; i < Q.rows(); ++i) {
                    t.knn(Q.row(i).transpose(), kk, idx, d2);
                    for (int j = 0; j < kk; ++j) {
                        I(i, j) = idx[static_cast<std::size_t>(j)];
                        D(i, j) = std::sqrt(d2[static_cast<std::size_t>(j)]);
                    }
                }
                if (single) return py::make_tuple(I.row(0).transpose().eval(), D.row(0).transpose().eval());
                return py::make_tuple(I, D);
            },
            py::arg("query"), py::arg("k"), "k nearest neighbours: (indices, distances), sorted by distance")
        .def(
            "radius",
            [](const Tree& t, const py::array_t<double, py::array::c_style | py::array::forcecast>& q, double r) {
                bool single = false;
                const RowsXd Q = as_queries(q, Dim, single);
                py::list out;
                std::vector<int> idx;
                std::vector<double> d2;
                for (Eigen::Index i = 0; i < Q.rows(); ++i) {
                    t.radius(Q.row(i).transpose(), r, idx, d2);
                    out.append(py::array_t<int>(static_cast<py::ssize_t>(idx.size()), idx.data()));
                }
                if (single) return py::object(out[0]);
                return py::object(out);
            },
            py::arg("query"), py::arg("r"), "indices within distance r (an array, or a list of arrays)");
}

template <int Dim>
py::tuple voxel_impl(const RowsXd& pts, double leaf, std::optional<Eigen::VectorXd> origin) {
    cloud::Point<Dim> o = cloud::Point<Dim>::Zero();
    if (origin) {
        if (origin->size() != Dim) throw std::invalid_argument("origin dimension mismatch");
        o = *origin;
    }
    const auto r = cloud::voxel_downsample<Dim>(cols_from_rows<Dim>(pts), leaf, o);
    Eigen::Matrix<long long, Eigen::Dynamic, Dim, Eigen::RowMajor> keys(static_cast<Eigen::Index>(r.keys.size()), Dim);
    for (std::size_t i = 0; i < r.keys.size(); ++i) keys.row(static_cast<Eigen::Index>(i)) = r.keys[i].transpose();
    return py::make_tuple(RowsD<Dim>(r.centroids.transpose()), py::array_t<int>(py::cast(r.voxel_of_point)),
                          py::array_t<int>(py::cast(r.count)), keys);
}

template <int Dim>
py::tuple normals_impl(const RowsXd& pts, int k, std::optional<Eigen::VectorXd> viewpoint, int leaf) {
    cloud::Point<Dim> vp = cloud::Point<Dim>::Zero();
    if (viewpoint) vp = *viewpoint;
    const cloud::Points<Dim> p = cols_from_rows<Dim>(pts);
    const cloud::KdTree<Dim> tree(p, leaf);
    const auto r = cloud::estimate_normals<Dim>(p, tree, k, vp);
    return py::make_tuple(RowsD<Dim>(r.normals.transpose()), r.variation);
}

template <int Dim>
py::tuple fit_plane_impl(const RowsXd& pts) {
    const auto f = cloud::fit_plane<Dim>(cols_from_rows<Dim>(pts));
    return py::make_tuple(Eigen::VectorXd(f.normal), Eigen::VectorXd(f.centroid), f.variation);
}

}  // namespace

void bind_cloud(py::module_ m) {
    bind_tree<2>(m, "KdTree2");
    bind_tree<3>(m, "KdTree3");
    m.def(
        "KdTree",
        [m](const RowsXd& pts, int leaf) -> py::object {
            if (pts.cols() == 2) return m.attr("KdTree2")(pts, leaf);
            if (pts.cols() == 3) return m.attr("KdTree3")(pts, leaf);
            throw std::invalid_argument("points must have shape (n, 2) or (n, 3)");
        },
        py::arg("points"), py::arg("leaf_size") = 10, "k-d tree over (n, 2) or (n, 3) points");

    m.def(
        "brute_force_knn",
        [](const RowsXd& pts, const Eigen::VectorXd& q, int k) {
            std::vector<int> idx;
            std::vector<double> d2;
            if (pts.cols() == 2)
                cloud::brute_force_knn<2>(cols_from_rows<2>(pts), q, k, idx, d2);
            else
                cloud::brute_force_knn<3>(cols_from_rows<3>(pts), q, k, idx, d2);
            for (auto& d : d2) d = std::sqrt(d);
            return py::make_tuple(idx, d2);
        },
        py::arg("points"), py::arg("query"), py::arg("k"));

    m.def(
        "voxel_downsample",
        [](const RowsXd& pts, double leaf, std::optional<Eigen::VectorXd> origin) {
            return pts.cols() == 2 ? voxel_impl<2>(pts, leaf, origin) : voxel_impl<3>(pts, leaf, origin);
        },
        py::arg("points"), py::arg("leaf"), py::arg("origin") = py::none(),
        "eqs. (4.2)-(4.3); returns (centroids, voxel_of_point, count, keys)");
    m.def(
        "estimate_normals",
        [](const RowsXd& pts, int k, std::optional<Eigen::VectorXd> viewpoint, int leaf) {
            return pts.cols() == 2 ? normals_impl<2>(pts, k, viewpoint, leaf) : normals_impl<3>(pts, k, viewpoint, leaf);
        },
        py::arg("points"), py::arg("k") = 10, py::arg("viewpoint") = py::none(), py::arg("leaf_size") = 10,
        "Algorithm 4.3; returns (normals, surface_variation)");
    m.def(
        "fit_plane", [](const RowsXd& pts) { return pts.cols() == 2 ? fit_plane_impl<2>(pts) : fit_plane_impl<3>(pts); },
        py::arg("points"), "PCA fit, eqs. (4.4)-(4.6): (normal, centroid, surface_variation)");

    m.def(
        "read_pcd",
        [](const std::string& path) {
            const auto d = cloud::read_pcd(path);
            py::dict out;
            for (std::size_t i = 0; i < d.fields.size(); ++i)
                out[py::str(d.fields[i])] = Eigen::VectorXd(d.values.row(static_cast<Eigen::Index>(i)).transpose());
            return out;
        },
        py::arg("path"), "read an ASCII or binary PCD file into {field: array}");
    m.def(
        "write_pcd",
        [](const std::string& path, const py::dict& fields, bool binary, const std::vector<std::string>& uint_fields) {
            cloud::PcdData d;
            Eigen::Index n = -1;
            std::vector<Eigen::VectorXd> cols;
            for (const auto& item : fields) {
                d.fields.push_back(py::cast<std::string>(item.first));
                cols.push_back(py::cast<Eigen::VectorXd>(item.second));
                if (n >= 0 && cols.back().size() != n) throw std::invalid_argument("all fields need the same length");
                n = cols.back().size();
            }
            d.values.resize(static_cast<Eigen::Index>(cols.size()), std::max<Eigen::Index>(n, 0));
            for (std::size_t i = 0; i < cols.size(); ++i) d.values.row(static_cast<Eigen::Index>(i)) = cols[i].transpose();
            cloud::write_pcd(path, d, binary, uint_fields);
        },
        py::arg("path"), py::arg("fields"), py::arg("binary") = true, py::arg("uint_fields") = std::vector<std::string>{});
}

}  // namespace sensor_fusion::python
