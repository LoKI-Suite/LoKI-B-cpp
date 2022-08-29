#include "LoKI-B/LookupTable.h"
#include <iostream>
#include <exception>
#include <stdexcept>

using loki::LookupTable;
using Vector = LookupTable::Vector;

int nerrors=0;

#define test_equal(a,b) { \
std::cout << "Testing: (" << #a ") == (" #b << ")\t"; \
try { \
if ((a)==(b)) \
    std::cout << " OK"; \
else { \
    std::cout << "FAILED"; \
    ++nerrors; \
} \
std::cout << std::endl; \
} \
catch(std::exception& exc) \
{ \
  std::cout << "AN EXCEPTION WAS THROWN: " << exc.what() << std::endl; \
} \
}

#define test_expr(expr) test_equal(expr,true)

void test1()
{
    Vector x(4);
    x << 1, 2, 3, 4;
    const Vector y(2.0*x);
    const LookupTable lut(x,y);

    // check table properties
    test_expr( lut.size()==4 );
    test_expr( lut.x().size()==4 );
    test_expr( lut.y().size()==4 );
    test_expr( lut.xMin()==1 );
    test_expr( lut.xMax()==4 );

    test_expr( lut.interpolate(1)==2 );
    test_expr( lut.interpolate(2)==4 );
    test_expr( lut.interpolate(3)==6 );
    test_expr( lut.interpolate(4)==8 );

    // interpolation:
    test_expr( lut.interpolate(1.1)==2.2 );
    test_expr( lut.interpolate(2.2)==4.4 );
    test_expr( lut.interpolate(3.8)==7.6 );

    // linear extrapolation:
    test_expr( lut.interpolate(0)==0 );
    test_expr( lut.interpolate(5)==10 );

    // clip to boundaries instead extrapolating:
    test_expr( lut.interpolate_or_clip(0,true,true)==2 );
    test_expr( lut.interpolate_or_clip(1,true,true)==2 );
    test_expr( lut.interpolate_or_clip(2.5,true,true)==5 ); // interpolation
    test_expr( lut.interpolate_or_clip(4,true,true)==8 );
    test_expr( lut.interpolate_or_clip(5,true,true)==8 );
    // clip to boundaries instead of extrapolating:
    test_expr( lut.interpolate_or_set(0,-42,42)==-42 );
    test_expr( lut.interpolate_or_set(1,-42,42)==2 );
    test_expr( lut.interpolate_or_set(2.5,-42,42)==5 ); // interpolation
    test_expr( lut.interpolate_or_set(4,-42,42)==8 );
    test_expr( lut.interpolate_or_set(5,-42,42)==42 );
}

int main()
{
    test1();

    std::cout << "Number of errors: " << nerrors << std::endl;
    return nerrors;
}
