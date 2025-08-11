#ifndef nuX_M1_MACRO_HXX
#define nuX_M1_MACRO_HXX

#define nuX_M1_SRC_EXPL 1  // explicit RHS
#define nuX_M1_SRC_IMPL 2  // implicit RHS (default)
#define nuX_M1_SRC_BOOST 3 // boost to fluid frame (approximate!)

#ifndef nuX_M1_SRC_METHOD
#define nuX_M1_SRC_METHOD nuX_M1_SRC_IMPL
#endif

#define SQ(X) ((X) * (X))

#endif
