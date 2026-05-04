project(
  'Uno',
  'cpp', # 'c'
  version: '2.7.2',
  license: 'MIT',
  meson_version: '>= 0.63.0',
  default_options: [
    'buildtype=release',
    'libdir=lib',
    'default_library=shared',
    'warning_level=0',
    'c_std=c99',
    'cpp_std=c++11',
  ],
)

# Compilers
# cc = meson.get_compiler('c')
cxx = meson.get_compiler('cpp')

# Remove messages about deprecated Intel compilers
# if cc.get_id() == 'intel'
#   add_global_arguments('-diag-disable=10441', language : 'c')
#   add_global_link_arguments('-diag-disable=10441', language : 'c')
# endif
# if cc.get_id() == 'intel-cl'
#   add_global_arguments('/Qdiag-disable=10441', language : 'c')
#   add_global_link_arguments('/Qdiag-disable=10441', language : 'c')
# endif
if cxx.get_id() == 'intel'
  add_global_arguments('-diag-disable=10441', language : 'cpp')
  add_global_link_arguments('-diag-disable=10441', language : 'cpp')
endif
if cxx.get_id() == 'intel-cl'
  add_global_arguments('/Qdiag-disable=10441', language : 'cpp')
  add_global_link_arguments('/Qdiag-disable=10441', language : 'cpp')
endif

# Options
build_examples = get_option('examples')
build_tests = get_option('tests')

libblas_name = get_option('libblas')
libblas_path = get_option('libblas_path')

liblapack_name = get_option('liblapack')
liblapack_path = get_option('liblapack_path')

libasl_name = get_option('libasl')
libasl_path = get_option('libasl_path')

libbqpd_name = get_option('libbqpd')
libbqpd_path = get_option('libbqpd_path')

libhighs_name = get_option('libhighs')
libhighs_path = get_option('libhighs_path')

libmumps_name = get_option('libmumps')
libmumps_path = get_option('libmumps_path')

libhsl_name = get_option('libhsl')
libhsl_path = get_option('libhsl_path')

libspral_name = get_option('libspral')
libspral_path = get_option('libspral_path')

# Dependencies
if libblas_path == []
    libblas = dependency(libblas_name, required : false)
endif
if libblas_path != [] or not libblas.found()
    libblas = cxx.find_library(libblas_name, dirs : libblas_path, required : true)
endif

if liblapack_path == []
  liblapack = dependency(liblapack_name, required : false)
endif
if liblapack_path != [] or not liblapack.found()
  liblapack = cxx.find_library(liblapack_name, dirs : liblapack_path, required : true)
endif

if libasl_path == []
  libasl = dependency(libasl_name, required : false)
endif

if libasl_path != [] or not libasl.found()
  libasl = cxx.find_library(libasl_name, dirs : libasl_path, required : false)
endif

if libbqpd_path == []
  libbqpd = dependency(libbqpd_name, required : false)
endif
if libbqpd_path != [] or not libbqpd.found()
  libbqpd = cxx.find_library(libbqpd_name, dirs : libbqpd_path, required : false)
endif

if libhighs_path == []
  libhighs = dependency(libhighs_name, required : false)
endif
if libhighs_path != [] or not libhighs.found()
  libhighs = cxx.find_library(libhighs_name, dirs : libhighs_path, required : false)
endif

if libmumps_path == []
  libmumps = dependency(libmumps_name, required : false)
endif
if libmumps_path != [] or not libmumps.found()
  libmumps = cxx.find_library(libmumps_name, dirs : libmumps_path, required : false)
endif

if libhsl_path == []
  libhsl = dependency(libhsl_name, required : false)
endif
if libhsl_path != [] or not libhsl.found()
  libhsl = cxx.find_library(libhsl_name, dirs : libhsl_path, required : false)
endif

if libspral_path == []
  libspral = dependency(libspral_name, required : false)
endif
if libspral_path != [] or not libspral.found()
  libspral = cxx.find_library(libspral_name, dirs : libspral_path, required : false)
endif

libuno_deps = [libblas, liblapack, libasl, libbqpd, libhighs, libmumps, libhsl, libspral]

# Headers
uno_headers = []
libuno_include = []
libuno_include += include_directories('interfaces', 'uno')
# libspral_include += libhwloc_include
# libspral_include += libblas_include

# Sources
libuno_src = []
uno_examples = []
uno_tests = []

subdir('uno')
subdir('interfaces')
# subdir('unotest')

# Library
libuno = library('uno',
                 sources : libuno_src,
                 dependencies : libuno_deps,
                 include_directories: libuno_include,
                 install : true)

# Headers
install_headers(uno_headers)

# # OMP-related environment variables
omp_env = ['OMP_CANCELLATION=true', 'OMP_PROC_BIND=true']

# Examples
if build_examples
  foreach uno_example: uno_examples
    name = uno_example[0]
    file = uno_example[1]
    test(name,
         executable(name, file,
                    link_with : libuno,
                    dependencies : libuno_deps,
                    include_directories : libuno_include,
                    install : false),
         suite : ['examples'],
         timeout : 300, is_parallel : false, env: omp_env)
  endforeach
endif

# Tests
if build_tests
  foreach uno_test: uno_tests
    name = uno_test[0]
    file = uno_test[1]
    test(name,
         executable(name, file,
                    link_with : libuno,
                    dependencies : libuno_deps,
                    include_directories: libuno_include,
                    install: false),
         suite : ['tests'],
         timeout : 300, is_parallel : false, env: omp_env)
  endforeach
endif
