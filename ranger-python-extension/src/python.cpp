// ranger python bindings
// Copyright (c) 2016 Sven Peter
// sven.peter@iwr.uni-heidelberg.de or mail@svenpeter.me
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
// documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
// Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <exception>
#include <string>
#include <sstream>

namespace py = pybind11;

#include "globals.h"
#include "ForestClassification.h"
#include "ForestRegression.h"
#include "ForestSurvival.h"
#include "ForestProbability.h"

namespace
{

template <typename features_t, typename labels_t> class DataPython : public Data
{
  public:
    DataPython(py::array_t<features_t> &data)
    {
        py::buffer_info info = data.request();

        if (info.ndim != 2)
            throw std::runtime_error("X vector needs to have shape (n_instances, n_features)");
        if (info.strides[1] != sizeof(features_t))
            throw std::runtime_error("X vector needs to be c continuous.");
        if (info.strides[0] != sizeof(labels_t) * info.shape[1])
            throw std::runtime_error("y vector needs to be c continuous.");

        ptr = (features_t *)info.ptr;

        num_rows = info.shape[0];
        num_cols = num_cols_real = info.shape[1];

        stride_rows = info.strides[0] / sizeof(features_t);

        have_labels = false;
        fill_varnames(variable_names);

        num_cols_no_sparse = num_cols + 2;
    }

    void addLabels(py::array_t<labels_t> &labels)
    {
        py::buffer_info info = labels.request();

        if (info.ndim != 1)
            throw std::runtime_error("y vector has more than one dimension");
        if (info.shape[0] != num_rows)
            throw std::runtime_error("y vector has invalid number of labels");
        if (info.strides[0] != sizeof(labels_t))
            throw std::runtime_error("y vector needs to be c continuous.");

        num_cols++;
        have_labels = true;
        ptr_labels = (labels_t *)info.ptr;

        variable_names.push_back("labels");
    }

    virtual ~DataPython()
    {
    }

    double get(size_t row, size_t col) const
    {
        if (have_labels && col == num_cols_real)
            return ptr_labels[row];
        return ptr[row * stride_rows + col];
    }

    void set(size_t col, size_t row, double value, bool &error)
    {
        ptr[row * stride_rows + col] = (features_t)value;
        error = false;
    }

    void reserveMemory()
    {
    }

    void fill_varnames(std::vector<std::string> &v)
    {
        for (unsigned int i = 0; i < num_cols; ++i) {
            std::stringstream ss;
            ss << "feature" << i;
            v.push_back(ss.str());
        }
    }

  private:
    features_t *ptr;
    size_t stride_rows;
    size_t num_cols_real;

    bool have_labels;
    labels_t *ptr_labels;
};

template <typename T> struct MemoryModeMap {
};
template <> struct MemoryModeMap<float> {
    static const MemoryMode memory_mode = MEM_FLOAT;
};
template <> struct MemoryModeMap<double> {
    static const MemoryMode memory_mode = MEM_DOUBLE;
};

template <typename T> class RandomForestClassifier : public ForestClassification
{
  public:
    size_t n_trees;
    size_t mtry;
    size_t n_threads;
    size_t min_node_size;
    double sample_fraction;
    double alpha;
    double minprop;

    RandomForestClassifier(size_t n_trees = 10, size_t mtry = 0, size_t n_threads = 1, size_t min_node_size = 1,
                           double sample_fraction = 1, double alpha = DEFAULT_ALPHA, double minprop = DEFAULT_MINPROP)
        : n_trees(n_trees), mtry(mtry), n_threads(n_threads), min_node_size(min_node_size),
          sample_fraction(sample_fraction), alpha(alpha), minprop(minprop)
    {
    }

    virtual ~RandomForestClassifier()
    {
    }

    void fit(py::array_t<T> X, py::array_t<size_t> y)
    {
        DataPython<T, size_t> *dataptr;
        std::unique_ptr<DataPython<T, size_t>> uptr_data(new DataPython<T, size_t>(X));

        dataptr = uptr_data.get();
        dataptr->addLabels(y);

        std::vector<std::string> unordered_variable_names;
        init("labels", MemoryModeMap<T>::memory_mode, dataptr, mtry, "dummy", n_trees, 0, n_threads, IMP_GINI,
             min_node_size, "statusVariableName", false, true, unordered_variable_names, false, DEFAULT_SPLITRULE,
             false, sample_fraction, alpha, minprop);
        run(false);
    }

    py::array_t<size_t> predict(py::array_t<T> X)
    {
        DataPython<T, size_t> *dataptr;
        std::unique_ptr<DataPython<T, size_t>> uptr_data(new DataPython<T, size_t>(X));

        dataptr = uptr_data.get();

        std::vector<std::string> unordered_variable_names;
        init("labels", MemoryModeMap<T>::memory_mode, dataptr, mtry, "dummy", n_trees, 0, n_threads, IMP_GINI,
             min_node_size, "statusVariableName", true, true, unordered_variable_names, false, DEFAULT_SPLITRULE, false,
             sample_fraction, alpha, minprop);
        run(false);

        py::buffer_info info = X.request();
        auto result = py::array(py::buffer_info(nullptr, sizeof(size_t), py::format_descriptor<size_t>::value, 1,
                                                {info.shape[0]}, {sizeof(size_t)}));
        py::buffer_info info_out = result.request();

        size_t *ptr_predictions = (size_t *)info_out.ptr;
        for (unsigned int i = 0; i < info.shape[0]; ++i) {
            ptr_predictions[i] = (size_t)predictions[i][0];
        }

        return result;
    }
};
}

namespace py = pybind11;

PYBIND11_PLUGIN(pyranger)
{
    py::module m("pyranger", "ranger python bindings");

    py::class_<RandomForestClassifier<double>>(m, "RandomForestClassifier")
        .def(py::init<size_t, size_t, size_t, size_t, double, double, double>(), py::arg("n_trees") = 10,
             py::arg("mtry") = 0, py::arg("n_threads") = 1, py::arg("min_node_size") = 1,
             py::arg("sample_fraction") = 1, py::arg("alpha") = DEFAULT_ALPHA, py::arg("minprop") = DEFAULT_MINPROP)
        .def("fit", &RandomForestClassifier<double>::fit)
        .def("predict", &RandomForestClassifier<double>::predict)
        .def_readwrite("n_trees", &RandomForestClassifier<double>::n_trees)
        .def_readwrite("mtry", &RandomForestClassifier<double>::mtry)
        .def_readwrite("n_threads", &RandomForestClassifier<double>::n_threads)
        .def_readwrite("min_node_size", &RandomForestClassifier<double>::min_node_size)
        .def_readwrite("sample_fraction", &RandomForestClassifier<double>::sample_fraction)
        .def_readwrite("alpha", &RandomForestClassifier<double>::alpha)
        .def_readwrite("minprop", &RandomForestClassifier<double>::minprop);

    return m.ptr();
}