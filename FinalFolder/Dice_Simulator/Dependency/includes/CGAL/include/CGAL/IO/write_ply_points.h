// Copyright (c) 2015  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot

#ifndef CGAL_POINT_SET_PROCESSING_WRITE_PLY_POINTS_H
#define CGAL_POINT_SET_PROCESSING_WRITE_PLY_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/IO/helpers.h>
#include <CGAL/IO/PLY.h>

#include <CGAL/property_map.h>
#include <CGAL/assertions.h>
#include <CGAL/Iterator_range.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/version.hpp>

#include <iostream>
#include <fstream>
#include <iterator>
#include <tuple>
#include <type_traits>

namespace CGAL {

namespace IO {

#ifdef DOXYGEN_RUNNING // Document some parts from Stream_support here for convenience
  /**
     \ingroup PkgPointSetProcessing3IOPly

     Generates a %PLY property handler to write 3D points. Points are
     written as 3 %PLY properties of type `FT` and named `x`, `y` and
     `z`. `FT` is `float` if the points use
     `CGAL::Simple_cartesian<float>` and `double` otherwise.

     \tparam PointMap the property map used to store points.

     \sa `write_PLY_with_properties()`
     \sa \ref IOStreamPLY
  */
  template <typename PointMap>
  std::tuple<PointMap, PLY_property<FT>, PLY_property<FT>, PLY_property<FT> >
  make_ply_point_writer(PointMap point_map);

  /**
     \ingroup PkgPointSetProcessing3IOPly

     Generates a %PLY property handler to write 3D normal
     vectors. Vectors are written as 3 %PLY properties of type `FT`
     and named `nx`, `ny` and `nz`. `FT` is `float` if the vectors use
     `CGAL::Simple_cartesian<float>` and `double` otherwise.

     \tparam VectorMap the property map used to store vectors.

     \sa `write_PLY_with_properties()`
     \sa \ref IOStreamPLY
  */
  template <typename VectorMap>
  std::tuple<VectorMap, PLY_property<FT>, PLY_property<FT>, PLY_property<FT> >
  make_ply_normal_writer(VectorMap normal_map);
#endif

/**
   \ingroup PkgPointSetProcessing3IOPly

   \brief writes the range of `points` with properties using \ref IOStreamPLY.

   Properties are handled through a variadic list of property
   handlers. A `PropertyHandler` can either be:

   - A `std::pair<PropertyMap, PLY_property<T> >` if the user wants
   to write a scalar value T as a %PLY property (for example, writing
   an `int` variable as an `int` %PLY property).

   - A `std::tuple<PropertyMap, PLY_property<T>...>` if the
   user wants to write a complex object as several %PLY
   properties. In that case, a specialization of `Output_rep` must
   be provided for `PropertyMap::value_type` that handles both ASCII
   and binary output (see `CGAL::IO::get_mode()`).

   \attention To write to a binary file, the flag `std::ios::binary` must be set during the creation
              of the `ofstream`, and the \link PkgStreamSupportEnumRef `IO::Mode` \endlink
              of the stream must be set to `BINARY`.

   \tparam PointRange is a model of `ConstRange`. The value type of
                      its iterator is the key type of the `PropertyMap` objects provided
                      within the `PropertyHandler` parameter.
   \tparam PropertyHandler handlers to recover properties.

   \returns `true` if writing was successful, `false` otherwise.

   \sa \ref IOStreamPLY
   \sa `make_ply_point_writer()`
   \sa `make_ply_normal_writer()`
*/
template <typename PointRange,
          typename ... PropertyHandler>
  bool write_PLY_with_properties(std::ostream& os, ///< output stream.
                                 const PointRange& points, ///< input point range.
                                 PropertyHandler&& ... properties) ///< parameter pack of property handlers
{
  CGAL_precondition(points.begin() != points.end());

  if(!os)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  // Write header
  os << "ply" << std::endl
     << ((get_mode(os) == BINARY) ? "format binary_little_endian 1.0" : "format ascii 1.0") << std::endl
     << "comment Generated by the CGAL library" << std::endl
     << "element vertex " << points.size() << std::endl;

  internal::output_property_header (os, std::forward<PropertyHandler>(properties)...);

  os << "end_header" << std::endl;

  // Write positions + normals
  for(typename PointRange::const_iterator it = points.begin(); it != points.end(); it++)
    internal::output_properties (os, it, std::forward<PropertyHandler>(properties)...);

  return !os.fail();
}

/**
   \ingroup PkgPointSetProcessing3IOPly

   \brief writes the range of `points` (positions + normals, if available) using \ref IOStreamPLY.

   \attention To write to a binary file, the flag `std::ios::binary` must be set during the creation
              of the `ofstream`, and the \link PkgStreamSupportEnumRef `IO::Mode` \endlink
              of the stream must be set to `BINARY`.

   \tparam PointRange is a model of `ConstRange`. The value type of
                      its iterator is the key type of the named parameter `point_map`.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param os output stream
   \param points input point range
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point range}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Vector_3`}
       \cgalParamDefault{If this parameter is omitted, normals are not written in the output stream.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd

     \cgalParamNBegin{stream_precision}
       \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
       \cgalParamType{int}
       \cgalParamDefault{the precision of the stream `os`}
       \cgalParamExtra{This parameter is only meaningful while using \ascii encoding.}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if writing was successful, `false` otherwise.

   \sa `write_PLY_with_properties()`
*/
template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_PLY(std::ostream& os,
               const PointRange& points,
               const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
               , std::enable_if_t<internal::is_Range<PointRange>::value>* = nullptr
#endif
               )
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // basic geometric types
  typedef Point_set_processing_3_np_helper<PointRange, CGAL_NP_CLASS> NP_helper;
  typedef typename NP_helper::Const_point_map PointMap;
  typedef typename NP_helper::Normal_map NormalMap;

  const bool has_normals = NP_helper::has_normal_map(points, np);

  PointMap point_map = NP_helper::get_const_point_map(points, np);
  NormalMap normal_map = NP_helper::get_normal_map(points, np);

  if(!os)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  set_stream_precision_from_NP(os, np);

  if(has_normals)
    return write_PLY_with_properties(os, points,
                                     make_ply_point_writer(point_map),
                                     make_ply_normal_writer(normal_map));

  return write_PLY_with_properties(os, points, make_ply_point_writer(point_map));
}

/**
   \ingroup PkgPointSetProcessing3IOPly

   \brief writes the range of `points` (positions + normals, if available) using \ref IOStreamPLY.

   \tparam PointRange is a model of `ConstRange`. The value type of
                      its iterator is the key type of the named parameter `point_map`.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param filename the path to the output file
   \param points input point range
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{use_binary_mode}
       \cgalParamDescription{indicates whether data should be written in binary (`true`) or in \ascii (`false`)}
       \cgalParamType{Boolean}
       \cgalParamDefault{`true`}
     \cgalParamNEnd

     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point range}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Vector_3`}
       \cgalParamDefault{If this parameter is omitted, normals are not written in the output file.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd

     \cgalParamNBegin{stream_precision}
       \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
       \cgalParamType{int}
       \cgalParamDefault{`6`}
       \cgalParamExtra{This parameter is only meaningful while using \ascii encoding.}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if writing was successful, `false` otherwise.

   \sa `write_PLY_with_properties()`
*/
template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_PLY(const std::string& filename,
               const PointRange& points,
               const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
               , std::enable_if_t<internal::is_Range<PointRange>::value>* = nullptr
#endif
               )
{
  const bool binary = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, internal_np::use_binary_mode), true);
  if(binary)
  {
    std::ofstream os(filename, std::ios::binary);
    CGAL::IO::set_mode(os, CGAL::IO::BINARY);
    return write_PLY(os, points, np);
  }
  else
  {
    std::ofstream os(filename);
    CGAL::IO::set_mode(os, CGAL::IO::ASCII);
    return write_PLY(os, points, np);
  }
}

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

/// \cond SKIP_IN_MANUAL

template <typename ForwardIterator,
          typename PointMap,
          typename VectorMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_ply_points_and_normals(), please update your code")
bool write_ply_points_and_normals(std::ostream& os, ///< output stream.
                                  ForwardIterator first, ///< first input point.
                                  ForwardIterator beyond, ///< past-the-end input point.
                                  PointMap point_map, ///< property map: value_type of OutputIterator -> Point_3.
                                  VectorMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return IO::write_PLY(os, points, parameters::point_map(point_map)
                                              .normal_map(normal_map));
}

template <typename ForwardIterator,
          typename VectorMap>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_ply_points_and_normals(), please update your code")
bool write_ply_points_and_normals(std::ostream& os, ///< output stream.
                                  ForwardIterator first, ///< first input point.
                                  ForwardIterator beyond, ///< past-the-end input point.
                                  VectorMap normal_map) ///< property map: value_type of OutputIterator -> Vector_3.
{
  CGAL::Iterator_range<ForwardIterator> points(first, beyond);
  return IO::write_PLY(os, points, parameters::normal_map (normal_map));
}

template <typename ForwardIterator,
          typename PointMap >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_ply_points(), please update your code")
bool write_ply_points(std::ostream& os, ///< output stream.
                      ForwardIterator first, ///< first input point.
                      ForwardIterator beyond, ///< past-the-end input point.
                      PointMap point_map) ///< property map: value_type of OutputIterator -> Point_3.
{
  CGAL::Iterator_range<ForwardIterator> points(first, beyond);
  return IO::write_PLY(os, points, parameters::point_map(point_map));
}

template <typename ForwardIterator >
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::write_ply_points(), please update your code")
bool write_ply_points(std::ostream& os, ///< output stream.
                      ForwardIterator first, ///< first input point.
                      ForwardIterator beyond) ///< past-the-end input point.
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return IO::write_PLY(os, points);
}

/// \endcond


template <typename PointRange,
          typename ... PropertyHandler>
CGAL_DEPRECATED bool write_ply_points_with_properties(std::ostream& os, ///< output stream.
                                                      const PointRange& points, ///< input point range.
                                                      PropertyHandler&& ... properties) ///< parameter pack of property handlers
{
  return IO::write_PLY_with_properties(os, points, std::forward<PropertyHandler>(properties)...);
}


template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED bool write_ply_points(std::ostream& os, const PointRange& points, const CGAL_NP_CLASS& np = parameters::default_values())
{
  return IO::write_PLY(os, points, np);
}

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // CGAL_POINT_SET_PROCESSING_WRITE_PLY_POINTS_H