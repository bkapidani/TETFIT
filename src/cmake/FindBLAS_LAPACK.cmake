#
# This source file is part of EMT, the ElectroMagneticTool.
#
# Copyright (C) 2013-2015, Matteo Cicuttin - matteo.cicuttin@uniud.it
# Department of Electrical Engineering, University of Udine
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the University of Udine nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

include(FindPackageHandleStandardArgs)

## Try to find the system blas library
find_library(SYSTEM_BLAS_LIBRARIES
    NAMES   blas cblas libblas
    PATHS   /opt/local/lib/
            /usr/lib/
            /usr/local/lib
            ${SYSTEM_BLAS_LIB_SEARCH_DIRS}
)

find_library(SYSTEM_LAPACK_LIBRARIES
    NAMES   lapack liblapack
    PATHS   /opt/local/lib/
            /usr/lib/
            /usr/local/lib
            ${SYSTEM_LAPACK_LIB_SEARCH_DIRS}
)

if(SYSTEM_BLAS_LIBRARIES AND SYSTEM_LAPACK_LIBRARIES)
    set (HAVE_SYSTEM_BLAS_LAPACK TRUE)
    set (SYSTEM_BLAS_LAPACK_LIBRARIES "${SYSTEM_BLAS_LIBRARIES}" "${SYSTEM_LAPACK_LIBRARIES}")
endif ()

if (SYSTEM_BLAS_LAPACK_LIBRARIES)
    set(BLAS_LAPACK_LIBRARIES "${SYSTEM_BLAS_LAPACK_LIBRARIES}")
endif ()

if (NVIDIA_BLAS_LIBRARIES)
    set(BLAS_LAPACK_LIBRARIES "${NVIDIA_BLAS_LIBRARIES}" "${BLAS_LAPACK_LIBRARIES}")
endif ()

FIND_PACKAGE_HANDLE_STANDARD_ARGS(BLAS_LAPACK DEFAULT_MSG BLAS_LAPACK_LIBRARIES)

