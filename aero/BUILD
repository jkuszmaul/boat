package(default_visibility = ['//visibility:public'])
cc_library(
  name="eigen_int",
  hdrs=["eigen_int.hpp"],
  deps=["@boost//:boost", "@eigen//:eigen"],
)
cc_library(
  name="dynamics",
  hdrs=["dynamics.h", "params.h"],
  srcs=["dynamics.cc"],
  deps=["@boost//:boost", "@eigen//:eigen"],
)
cc_library(
  name="control",
  hdrs=["control.h"],
  srcs=["control.cc"],
  deps=["@boost//:boost",
        "@eigen//:eigen",
        ":dynamics"],
)
cc_library(
  name="sim",
  hdrs=["sim.h"],
  deps=["@boost//:boost",
        "@eigen//:eigen",
        ":dynamics",
        ":control",
        ":eigen_int",
        "@glog//:glog",
        "@gflags//:gflags"],
  linkopts = ["-pthread"],
)
cc_binary(
  name="sim_main",
  srcs=["sim_main.cc"],
  deps=[":sim",
        "@glog//:glog",
        "@gflags//:gflags"],
)
cc_binary(
  name="integrate",
  srcs=["integrate.cc"],
  deps=["@boost//:boost",
        "@eigen//:eigen",
        ":eigen_int"],
)
