/* Simplified SWIG interface for B-spline debugging */

%{
#include <ql/math/interpolations/bsplineinterpolation/debug_inspection.hpp>
%}

namespace QuantLib {

    class BSplineDebugInspector {
    public:
        static Size getNumberOfSegments(const BSplineStructure& structure);
        static std::vector<Real> getSolutionCoefficients(const BSplineInterpolation& interp);
        static void printBasicInfo(const BSplineStructure& structure);
    };

}