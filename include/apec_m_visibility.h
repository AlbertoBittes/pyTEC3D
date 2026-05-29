#ifndef APEC_M_VISIBILITY_H
#define APEC_M_VISIBILITY_H

#if defined(__GNUC__) || defined(__clang__)
#define APEC_M_HIDDEN __attribute__((visibility("hidden")))
#else
#define APEC_M_HIDDEN
#endif

#endif // APEC_M_VISIBILITY_H