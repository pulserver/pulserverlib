# csrc/CMakeLists.txt -- pulseqlib library build
#
# Public API surface (always built):
#   pulseqlib_error/math/parse/dedup/structure/core/waveforms/safety/
#   getters/freqmod/protocol/trajectory
#
# Opt-IN components (toggled by cache option below):
#   pulseqlib_cache.c          -- vendor-neutral binary cache format writer
#                                 (PULSEQLIB_BUILD_CACHE, default ON)
#
# Opt-IN components (toggled by seqdesc option below):
#   pulseqlib_seqdesc.c        -- sequence description for Bloch simulation
#   pulseqlib_cache_seqdesc.c  -- seqdesc binary cache writer
#                                 (PULSEQLIB_BUILD_SEQDESC, default OFF for the
#                                 static library; the pybind11 wrapper turns it
#                                 ON because the pge Python API exposes
#                                 pulseqlib_get_sequence_{description,parameters})
#
# These options exist so external vendor consumers (e.g., a PSD that wants the
# safety subsystem only) can compile a stripped variant. Default-ON consumers
# (pge, mrdserver) are byte-identical to the previous behaviour.

option(PULSEQLIB_BUILD_CACHE   "Build the binary cache writer (pulseqlib_cache.c)" ON)
option(PULSEQLIB_BUILD_SEQDESC "Build the sequence-description / Bloch-sim helpers" OFF)

# Always-on library sources
set(PULSEQLIB_SOURCES
    pulseqlib_error.c
    pulseqlib_math.c
    pulseqlib_parse.c
    pulseqlib_dedup.c
    pulseqlib_structure.c
    pulseqlib_core.c
    pulseqlib_waveforms.c
    pulseqlib_safety.c
    pulseqlib_getters.c
    pulseqlib_freqmod.c
    pulseqlib_protocol.c
    pulseqlib_trajectory.c
)

if(PULSEQLIB_BUILD_CACHE)
    list(APPEND PULSEQLIB_SOURCES pulseqlib_cache.c)
endif()

if(PULSEQLIB_BUILD_SEQDESC)
    list(APPEND PULSEQLIB_SOURCES pulseqlib_seqdesc.c)
    if(PULSEQLIB_BUILD_CACHE)
        list(APPEND PULSEQLIB_SOURCES pulseqlib_cache_seqdesc.c)
    endif()
endif()

# MD5 source (needed by pulseqlib_parse for signature verification)
set(MD5_SOURCES
    external_md5.c
)

# kissfft sources (co-located for self-containment)
set(KISSFFT_SOURCES
    external_kiss_fft.c
    external_kiss_fftr.c
)

# Public headers (for IDE / install)
set(PULSEQLIB_HEADERS
    pulseqlib.h
    pulseqlib_config.h
    pulseqlib_types.h
)

# Build the library
add_library(pulseqlib STATIC
    ${PULSEQLIB_SOURCES}
    ${KISSFFT_SOURCES}
    ${MD5_SOURCES}
)

# Include the current directory so #include "pulseqlib.h" works
target_include_directories(pulseqlib
    PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Bridge module (POSIX-only, uses fork/exec/pipe)
set(BRIDGE_SOURCES
    pulseqlib_bridge.c
)
add_library(pulseqlib_bridge STATIC ${BRIDGE_SOURCES})
target_link_libraries(pulseqlib_bridge PUBLIC pulseqlib)
target_include_directories(pulseqlib_bridge
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)