#ifndef NUX_M1_MACRO_HXX
#define NUX_M1_MACRO_HXX

#define NUX_M1_SRC_EXPL 1  // explicit RHS
#define NUX_M1_SRC_IMPL 2  // implicit RHS (default)
#define NUX_M1_SRC_BOOST 3 // boost to fluid frame (approximate!)

#ifndef NUX_M1_SRC_METHOD
#define NUX_M1_SRC_METHOD NUX_M1_SRC_IMPL
#endif

#define SQ(X) ((X) * (X))

#endif
