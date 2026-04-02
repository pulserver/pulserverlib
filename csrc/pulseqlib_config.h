/* pulseqlib_config.h -- platform configuration
 *
 * This header MUST be included (directly or transitively) before any
 * other pulseqlib header or source file.  It provides:
 *
 *   PULSEQLIB_VENDOR_*  -- vendor ID constants (used at runtime)
 *   PULSEQLIB_VENDOR    -- compile-time default vendor
 *   PULSEQLIB_ALLOC     -- heap allocator  (default: malloc)
 *   PULSEQLIB_FREE      -- heap deallocator (default: free)
 */

#ifndef PULSEQLIB_CONFIG_H
#define PULSEQLIB_CONFIG_H

#include <stdlib.h>

/* ------------------------------------------------------------------ */
/*  Vendor identifiers (runtime constants)                            */
/* ------------------------------------------------------------------ */
#define PULSEQLIB_VENDOR_SIEMENS        1
#define PULSEQLIB_VENDOR_GEHC           2
#define PULSEQLIB_VENDOR_PHILIPS        3
#define PULSEQLIB_VENDOR_UNITED_IMAGING 4
#define PULSEQLIB_VENDOR_BRUKER         5

/* Compile-time default (overrideable via -DPULSEQLIB_VENDOR=N) */
#ifndef PULSEQLIB_VENDOR
#define PULSEQLIB_VENDOR PULSEQLIB_VENDOR_GEHC
#endif

/* ------------------------------------------------------------------ */
/*  Acoustic peak-detection defaults                                  */
/* ------------------------------------------------------------------ */

#ifndef PULSEQLIB_PEAK_LOG10_THRESHOLD_DEFAULT
#define PULSEQLIB_PEAK_LOG10_THRESHOLD_DEFAULT 2.25f
#endif

#ifndef PULSEQLIB_PEAK_NORM_SCALE_DEFAULT
#define PULSEQLIB_PEAK_NORM_SCALE_DEFAULT 10.0f
#endif

#ifndef PULSEQLIB_PEAK_EPS_DEFAULT
#define PULSEQLIB_PEAK_EPS_DEFAULT 1e-30f
#endif

/* ------------------------------------------------------------------ */
/*  Allocator overrides                                               */
/* ------------------------------------------------------------------ */

/*
 * Override PULSEQLIB_ALLOC / PULSEQLIB_FREE *before* including this
 * header to use vendor-specific allocators, e.g.:
 *
 *   #define PULSEQLIB_ALLOC(sz)  MyVendorAlloc(sz)
 *   #define PULSEQLIB_FREE(ptr)  MyVendorFree(ptr)
 *   #include "pulseqlib_config.h"
 */
#ifndef PULSEQLIB_ALLOC
#define PULSEQLIB_ALLOC(sz) malloc(sz)
#endif

#ifndef PULSEQLIB_FREE
#define PULSEQLIB_FREE(ptr) free(ptr)
#endif

#endif /* PULSEQLIB_CONFIG_H */
