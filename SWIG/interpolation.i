/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2002, 2003 Ferdinando Ametrano
 Copyright (C) 2003, 2004, 2008 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2018, 2020 Matthias Lungwitz
 
 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#ifndef quantlib_interpolation_i
#define quantlib_interpolation_i

%include termstructures.i
%include linearalgebra.i
%include optimizers.i

%{
using QuantLib::Discount;
using QuantLib::ZeroYield;
using QuantLib::ForwardRate;
using QuantLib::TermForwardRate;
using QuantLib::RateTime;
using QuantLib::BSplineModel;
using QuantLib::InterestRateIndex;
using QuantLib::YieldTermStructure;
%}

%{
// safe versions which copy their arguments
template <class I>
class SafeInterpolation {
  public:
    SafeInterpolation(const Array& x, const Array& y)
    : x_(x), y_(y), f_(x_.begin(),x_.end(),y_.begin()) {}
    Real operator()(Real x, bool allowExtrapolation=false) {
        return f_(x, allowExtrapolation);
    }
    Array x_, y_;
    I f_;
};
%}

%define make_safe_interpolation(T,Alias)
%{
typedef SafeInterpolation<QuantLib::T> Safe##T;
%}
%rename(Alias) Safe##T;
class Safe##T {
    #if defined(SWIGCSHARP)
    %rename(call) operator();
    #endif
  public:
    Safe##T(const Array& x, const Array& y);
    Real operator()(Real x, bool allowExtrapolation=false);
};
%enddef

make_safe_interpolation(LinearInterpolation,LinearInterpolation);
make_safe_interpolation(LogLinearInterpolation,LogLinearInterpolation);

make_safe_interpolation(BackwardFlatInterpolation,BackwardFlatInterpolation);
make_safe_interpolation(ForwardFlatInterpolation,ForwardFlatInterpolation);

make_safe_interpolation(CubicNaturalSpline,CubicNaturalSpline);
make_safe_interpolation(LogCubicNaturalSpline,LogCubicNaturalSpline);
make_safe_interpolation(MonotonicCubicNaturalSpline,MonotonicCubicNaturalSpline);
make_safe_interpolation(MonotonicLogCubicNaturalSpline,MonotonicLogCubicNaturalSpline);

make_safe_interpolation(KrugerCubic,KrugerCubic);
make_safe_interpolation(KrugerLogCubic,KrugerLogCubic);

make_safe_interpolation(FritschButlandCubic,FritschButlandCubic);
make_safe_interpolation(FritschButlandLogCubic,FritschButlandLogCubic);

make_safe_interpolation(Parabolic,Parabolic);
make_safe_interpolation(LogParabolic,LogParabolic);
make_safe_interpolation(MonotonicParabolic,MonotonicParabolic);
make_safe_interpolation(MonotonicLogParabolic,MonotonicLogParabolic);

make_safe_interpolation(LagrangeInterpolation,LagrangeInterpolation); 

%{
using QuantLib::Interpolation;
using QuantLib::BSplineModel;
using QuantLib::BSplineSegment;
using QuantLib::BSplineStructure;
using QuantLib::SplineConstraints;
using QuantLib::BSplineInterpolation;
using QuantLib::SpreadedInterpolationModel;
using QuantLib::SpreadedInterpolation;
using QuantLib::Linear;
using QuantLib::InterpolationMode;
using QuantLib::ModeSpan;
using QuantLib::generateVoronoiModeSpans;
%}

namespace std {
    %template(ConstraintTypeVector) std::vector<SplineConstraints::ConstraintType>;
}

// Create flattened container class for SplineConstraints enums
%rename (SplineConstraints_ConstraintType) SplineConstraints::ConstraintType;

// Expose the SplineConstraints class
%shared_ptr(SplineConstraints);

// Commented out to avoid SWIG shadowing warnings - conversion handled in Python
// %extend SplineConstraints {
//     SplineConstraints(Size nVariables,
//                       const std::vector<std::vector<Real>>& P,
//                       const std::vector<std::vector<Real>>& A,
//                       const std::vector<Real>& b,
//                       const std::vector<Real>& c,
//                       const std::vector<int>& constraintTypesInt,
//                       bool fitData,
//                       double epsAbsolute = 1e-12,
//                       double epsRelative = 1e-12,
//                       double epsInfeasible = 1e-13) {
//         // Convert integers to ConstraintType enum
//         std::vector<SplineConstraints::ConstraintType> types;
//         for (int typeInt : constraintTypesInt) {
//             types.push_back(static_cast<SplineConstraints::ConstraintType>(typeInt));
//         }
//         return new SplineConstraints(nVariables, P, A, b, c, types, fitData, 
//                                     epsAbsolute, epsRelative, epsInfeasible);
//     }
// }

class SplineConstraints {
public:
    // Expose the ConstraintType enum
    enum class ConstraintType {
        Equal,
        LessEqual
    };

    // New SCS-ordered constructor (constraints must be in SCS order: equalities first, then inequalities)
    SplineConstraints(Size nVariables,
                        const std::vector<std::vector<Real>>& P,
                        const std::vector<std::vector<Real>>& A,
                        const std::vector<Real>& b,
                        const std::vector<Real>& c,
                        Size numEqualities,
                        Size numInequalities,
                        bool fitData = false,
                        double epsAbsolute = 1e-12,
                        double epsRelative = 1e-12,
                        double epsInfeasible = 1e-13);

    SplineConstraints(Size nVariables = 0,
                        const std::vector<std::vector<Real>>& P = {},
                        const std::vector<std::vector<Real>>& A = {},
                        const std::vector<Real>& b = {},
                        const std::vector<Real>& c = {},
                        const std::vector<ConstraintType>& constraintTypes = {},
                        bool fitData = false,
                        double epsAbsolute = 1e-12,
                        double epsRelative = 1e-12,
                        double epsInfeasible = 1e-13);

    Integer get_num_variables() const;
    Integer get_num_constraints() const;
    Integer get_num_parameters() const;

    std::vector<std::vector<Real>> get_p_matrix() const;
    std::vector<std::vector<Real>> get_a_matrix() const;

    std::vector<Real> get_b_vector() const;
    std::vector<Real> get_c_vector() const;
    std::vector<Real> get_parameters() const;
    std::vector<ConstraintType> get_constraint_types() const;
    
    // Inspection methods for debugging constraint ordering
    std::vector<int> get_permutation() const;
    std::vector<std::vector<Real>> get_reordered_a_matrix() const;
    std::vector<Real> get_reordered_b_vector() const;
    bool is_ordered() const;
    Size get_num_equalities() const;
    Size get_num_inequalities() const;
};


// Expose the BSplineSegment class
%shared_ptr(BSplineSegment);

// Create flattened container classes following AnalyticHestonEngine pattern
%rename (BSpline_Side) BSplineSegment::Side;
%rename (BSpline_Type) BSplineSegment::InterpolationType;
%rename (BSpline_Transform) BSplineSegment::InterpolationTransform;
%rename (BSpline_Smoothness) BSplineSegment::InterpolationSmoothness;

class BSplineSegment {
public:
    // Simple enums like AnalyticHestonEngine for SWIG compatibility
    enum Side { 
        SideLeft, SideRight, SideAverage, SideActual, SideInside, SideDefault 
    };
    
    enum InterpolationType {
        TypeLinear,
        TypeLinearRT,
        TypeQuadraticSplineRT,
        TypeCubicSpline,
        TypeCubicSplineRT,
        TypeHermite,
        TypeHermiteRT,
        TypePiecewiseQuadraticContFwd,
        TypeConvexMonotone,
        TypeShortest,
        TypeClosest,
        TypeFollowing,
        TypePreceding
    };
    
    enum InterpolationTransform {
        TransformDefault,
        TransformLog,
        TransformExp,
        TransformRateTime,
        TransformRateTimeAnnualToContinuous,
        TransformContinuousToAnnual,
        TransformContinuousToSimple
    };
    
    enum InterpolationSmoothness {
        SmoothnessDiscontinuous, // Internal knots repeated k times for k-th degree spline, aka C^{-1}
        SmoothnessContinuous,    // aka C^0
        SmoothnessContinuouslyDifferentiable,      // aka C^1
        SmoothnessTwiceContinuouslyDifferentiable, // aka C^2
        SmoothnessHermite, // Internal knots are double, means C^{k-2} for a k-th degree spline
        SmoothnessDefault  // Internal knots are simple, means C^{k-1} for a k-th degree spline
    };
    
    // Compatibility typedefs for existing code
    using SideEnum = Side;
    using InterpolationTypeEnum = InterpolationType;
    using InterpolationTransformEnum = InterpolationTransform;
    using InterpolationSmoothnessEnum = InterpolationSmoothness;

    BSplineSegment(
        const std::vector<Real>& simpleKnots,
        Integer degree,
        const std::vector<Integer>& knotIndices,
        InterpolationSmoothnessEnum smoothness = SmoothnessDefault,
        InterpolationTransformEnum interpolationTransform = TransformDefault,
        SideEnum side = SideRight,
        Size requiredPoints = 1,
        bool isGlobal = true);


    #if defined(SWIGPYTHON)
        %rename(get_spline_segments) getSplineSegmentsSwig;
        %rename(evaluate_all) evaluateAllSwig;
        %rename(get_num_variables) getNumVariablesSwig;
    #endif

    std::vector<Real> evaluateAllSwig(Real x, Integer degree = -1, QuantLib::BSplineSegment::SideEnum side = QuantLib::BSplineSegment::SideRight) const;

    std::pair<Real, Real> range() const;
    std::vector<Real> knots() const;
    Size degree() const;
    Real value(const std::vector<Real>& coefficients, Real t, Integer nu, SideEnum side);
    std::vector<Real> value_functional(Real t, SideEnum side) const;
    std::vector<Real> derivative_functional(Real t, Integer nu = 1, Integer degree = -1, Real x0 = 0.0, SideEnum side = SideDefault) const;
    std::vector<std::vector<Real>> derivative_matrix(Integer nu = 1, Integer degree = -1, bool differenceOperator = false) const;
    std::vector<std::vector<Real>> anti_derivative_matrix(Integer nu = -1, Integer degree = -1,
                                                         const std::vector<Real> t0 = {},
                                                         bool differenceOperator = false) const;
    std::vector<std::vector<Real>> single_derivative_matrix(Integer degree = -1, bool differenceOperator = false) const;
    std::vector<std::vector<Real>> single_anti_derivative_matrix(Integer degree = -1, Real t0 = 0.0, bool differenceOperator = false) const;

    std::vector<Real> get_simple_knots() const;
    std::vector<Integer> get_knot_indices() const;
    Natural getNumVariablesSwig() const;

    InterpolationTransformEnum interpolationTransform() const;
    InterpolationSmoothnessEnum interpolationSmoothness() const;
    SideEnum side() const;
    std::string interpolation_transform() const;
    std::string interpolation_smoothness() const;
    std::string side_str() const;
};

namespace std {
    %template(BSplineSegmentVector) std::vector<ext::shared_ptr<BSplineSegment>>;
}

// Expose the BSplineStructure class
%shared_ptr(BSplineStructure);
class BSplineStructure {
  public:
    BSplineStructure(
            const std::vector<ext::shared_ptr<BSplineSegment>>& splineSegments,
            const ext::shared_ptr<SplineConstraints>& splineConstraints,
            bool useSegmentNodes = false, bool rejectZeroNode = true);

#if defined(SWIGPYTHON)
    %rename(get_num_variables) getNumVariablesSwig;
    %rename(get_spline_segments) getSplineSegmentsSwig;
    %rename(evaluate_all) evaluateAllSwig;
#endif

    Integer getNumVariablesSwig() const;
    std::vector<ext::shared_ptr<BSplineSegment>> getSplineSegmentsSwig() const;
    std::vector<Real> evaluateAllSwig(Real x, BSplineSegment::SideEnum side) const;
};

// Expose InterpolationMode enum
enum class InterpolationMode {
    HARD,   // Exact interpolation constraints
    LS,     // Least squares fitting
    AUTO    // Automatically determine based on DOF
};

// Expose ModeSpan struct
struct ModeSpan {
    Real start;
    Real end;
    InterpolationMode mode;
    
    ModeSpan(Real s, Real e, InterpolationMode m);
};

namespace std {
    %template(ModeSpanVector) std::vector<ModeSpan>;
}

// Expose Voronoi span generation function
std::vector<ModeSpan> generateVoronoiModeSpans(
    const std::vector<Real>& dataPoints,
    Real domainStart,
    Real domainEnd,
    Real densityThreshold = 2.0,
    Real boundaryWidth = 0.1
);

// Expose the BSplineModel class
%shared_ptr(BSplineModel);
class BSplineModel {
public:
    BSplineModel(ext::shared_ptr<BSplineStructure>& splineStructure);

    ext::shared_ptr<BSplineInterpolation> interpolate(const std::vector<Real>& x, const std::vector<Real> y) const;
    Real getStartPoint() const;
    Real getEndPoint() const;
    ext::shared_ptr<BSplineStructure> get_structure() const;
};

// Expose the Interpolation class, but without constructor
%shared_ptr(Interpolation)
class Interpolation {
  public:
    Real operator()(Real x, bool allowExtrapolation = false) const;
    Real primitive(Real x, bool allowExtrapolation = false) const;
    Real derivative(Real x, bool allowExtrapolation = false) const;
    Real secondDerivative(Real x, bool allowExtrapolation = false) const;
    Real xMin() const;
    Real xMax() const;
    bool isInRange(Real x) const;
    void update();
};

// Expose the BSplineInterpolation class
%shared_ptr(BSplineInterpolation);
class BSplineInterpolation : public Interpolation {
public:
    BSplineInterpolation(
      const std::vector<Time>& x,
      const std::vector<Real>& y,
      const ext::shared_ptr<BSplineStructure>& splineStructure
    );
    Real operator()(Real x, bool allowExtrapolation = false) const;
    Real primitive(Real x, bool allowExtrapolation = false) const;
    Real derivative(Real x, bool allowExtrapolation = false) const;
    Real secondDerivative(Real x, bool allowExtrapolation = false) const;
    Real xMin() const;
    Real xMax() const;
    bool isInRange(Real x) const;
    void update();

    ext::shared_ptr<BSplineStructure> get_structure() const;
    std::vector<Real> get_coefficients() const;
};

%inline %{
    ext::shared_ptr<BSplineInterpolation> as_bspline_interpolation(const ext::shared_ptr<Interpolation>& ip) {
        return ext::dynamic_pointer_cast<BSplineInterpolation>(ip);
    }
%}

// Create Python-accessible enum containers
%pythoncode %{
# Create enum containers that mimic the AnalyticHestonEngine_OptimalAlpha pattern
class BSpline_Side:
    Left = BSplineSegment.SideLeft
    Right = BSplineSegment.SideRight
    Average = BSplineSegment.SideAverage
    Actual = BSplineSegment.SideActual
    Inside = BSplineSegment.SideInside
    Default = BSplineSegment.SideDefault

class BSpline_Type:
    Linear = BSplineSegment.TypeLinear
    LinearRT = BSplineSegment.TypeLinearRT
    QuadraticSplineRT = BSplineSegment.TypeQuadraticSplineRT
    CubicSpline = BSplineSegment.TypeCubicSpline
    CubicSplineRT = BSplineSegment.TypeCubicSplineRT
    Hermite = BSplineSegment.TypeHermite
    HermiteRT = BSplineSegment.TypeHermiteRT
    PiecewiseQuadraticContFwd = BSplineSegment.TypePiecewiseQuadraticContFwd
    ConvexMonotone = BSplineSegment.TypeConvexMonotone
    Shortest = BSplineSegment.TypeShortest
    Closest = BSplineSegment.TypeClosest
    Following = BSplineSegment.TypeFollowing
    Preceding = BSplineSegment.TypePreceding

class BSpline_Transform:
    Default = BSplineSegment.TransformDefault
    Log = BSplineSegment.TransformLog
    Exp = BSplineSegment.TransformExp
    RateTime = BSplineSegment.TransformRateTime
    RateTimeAnnualToContinuous = BSplineSegment.TransformRateTimeAnnualToContinuous
    ContinuousToAnnual = BSplineSegment.TransformContinuousToAnnual
    ContinuousToSimple = BSplineSegment.TransformContinuousToSimple

class BSpline_Smoothness:
    Default = BSplineSegment.SmoothnessDefault
    Discontinuous = BSplineSegment.SmoothnessDiscontinuous
    Continuous = BSplineSegment.SmoothnessContinuous
    ContinuouslyDifferentiable = BSplineSegment.SmoothnessContinuouslyDifferentiable
    TwiceContinuouslyDifferentiable = BSplineSegment.SmoothnessTwiceContinuouslyDifferentiable
    Hermite = BSplineSegment.SmoothnessHermite

# Create flat middle container class (like AnalyticHestonEngine_OptimalAlpha pattern)
# SWIG generates SplineConstraints.SplineConstraints_ConstraintType_Equal at class level
class SplineConstraints_ConstraintType:
    Equal = SplineConstraints.SplineConstraints_ConstraintType_Equal
    LessEqual = SplineConstraints.SplineConstraints_ConstraintType_LessEqual

# Also create nested access for SplineConstraints.ConstraintType.Equal pattern
class _ConstraintType:
    Equal = SplineConstraints.SplineConstraints_ConstraintType_Equal
    LessEqual = SplineConstraints.SplineConstraints_ConstraintType_LessEqual

# Attach as nested class to SplineConstraints
SplineConstraints.ConstraintType = _ConstraintType
%}


// %{
// // safe versions which copy their arguments
// class SafeBSplineInterpolation {
//   public:
//     SafeBSplineInterpolation(const std::vector<Time>& x,
//       const std::vector<Real>& y,  const ext::shared_ptr<BSplineStructure>& splineStructure)
//     : x_(x), y_(y), f_(x_.begin(), x_.end(), y_.begin(), splineStructure) {}
//     Real operator()(Real x, bool allowExtrapolation=false) {
//         return f_(x, allowExtrapolation);
//     }
//     std::vector<Time> x_;
//     std::vector<Real> y_;
//     BSplineInterpolation f_;
// };
// %}

// %rename(BSplineInterpolation) SafeBSplineInterpolation;
// class SafeBSplineInterpolation {
//     #if defined(SWIGCSHARP)
//     %rename(call) operator();
//     #endif
//   public:
//     SafeBSplineInterpolation(const std::vector<Time>& x,
//       const std::vector<Real>& y, const ext::shared_ptr<BSplineStructure>& splineStructure);
//     Real operator()(Real x, bool allowExtrapolation=false);
// };

// %extend SafeBSplineInterpolation {
//     Real derivative(Real x, bool extrapolate = false) {
//         return self->f_.derivative(x,extrapolate);
//     }
//     Real secondDerivative(Real x, bool extrapolate = false) {
//         return self->f_.secondDerivative(x,extrapolate);
//     }
//     Real primitive(Real x, bool extrapolate = false) {
//         return self->f_.primitive(x,extrapolate);
//     }
// };

// Expose the SpreadedInterpolationModel class
// %shared_ptr(SpreadedInterpolationModel);
%shared_ptr(SpreadedInterpolationModel<ZeroYield, Linear>);
%shared_ptr(SpreadedInterpolationModel<ZeroYield, BSplineModel>);
%shared_ptr(SpreadedInterpolationModel<TermForwardRate, BSplineModel>);

template <class Traits, class Interpolator>
class SpreadedInterpolationModel {
  public:
      // SpreadedInterpolationModel(Interpolator& factory, const ext::shared_ptr<YieldTermStructure>& baseCurve,
      //   const std::optional<ext::shared_ptr<InterestRateIndex>>& index = std::nullopt);

      %extend{
        SpreadedInterpolationModel(Interpolator& factory, const ext::shared_ptr<YieldTermStructure>& baseCurve,
          const ext::shared_ptr<InterestRateIndex>& index) {
          return new SpreadedInterpolationModel<Traits, Interpolator>(factory, baseCurve, index);
        }

        SpreadedInterpolationModel(Interpolator& factory, const ext::shared_ptr<YieldTermStructure>& baseCurve) {
          return new SpreadedInterpolationModel<Traits, Interpolator>(factory, baseCurve, std::nullopt);
        }
      };
};

%template(SpreadedZeroLinearModel) SpreadedInterpolationModel<ZeroYield, Linear>;
%template(SpreadedZeroBSplineModel) SpreadedInterpolationModel<ZeroYield, BSplineModel>;
%template(SpreadedTermForwardBSplineModel) SpreadedInterpolationModel<TermForwardRate, BSplineModel>;
%define extend_spline(T)
%extend Safe##T {
    Real derivative(Real x, bool extrapolate = false) {
        return self->f_.derivative(x,extrapolate);
    }
    Real secondDerivative(Real x, bool extrapolate = false) {
        return self->f_.secondDerivative(x,extrapolate);
    }
    Real primitive(Real x, bool extrapolate = false) {
        return self->f_.primitive(x,extrapolate);
    }
}
%enddef

extend_spline(CubicNaturalSpline);
extend_spline(LogCubicNaturalSpline);
extend_spline(MonotonicCubicNaturalSpline);
extend_spline(MonotonicLogCubicNaturalSpline);

extend_spline(KrugerCubic);
extend_spline(KrugerLogCubic);

extend_spline(FritschButlandCubic);
extend_spline(FritschButlandLogCubic);

extend_spline(Parabolic);
extend_spline(LogParabolic);
extend_spline(MonotonicParabolic);
extend_spline(MonotonicLogParabolic);

%{
// safe versions which copy their arguments
template <class I>
class SafeInterpolation2D {
  public:
    SafeInterpolation2D(const Array& x, const Array& y, const Matrix& m)
    : x_(x), y_(y), m_(m), f_(x_.begin(),x_.end(),y_.begin(),y_.end(),m_) {}
    Real operator()(Real x, Real y, bool allowExtrapolation=false) {
        return f_(x,y, allowExtrapolation);
    }
  protected:
    Array x_, y_;
    Matrix m_;
    I f_;
};
%}

%define make_safe_interpolation2d(T,Alias)
%{
typedef SafeInterpolation2D<QuantLib::T> Safe##T;
%}
%rename(Alias) Safe##T;
class Safe##T {
    #if defined(SWIGCSHARP)
    %rename(call) operator();
    #endif
  public:
    Safe##T(const Array& x, const Array& y, const Matrix& m);
    Real operator()(Real x, Real y, bool allowExtrapolation=false);
};
%enddef

make_safe_interpolation2d(BilinearInterpolation,BilinearInterpolation);
make_safe_interpolation2d(BicubicSpline,BicubicSpline);


// interpolation traits

%{
using QuantLib::CubicInterpolation;
using QuantLib::MixedInterpolation;
using QuantLib::BackwardFlat;
using QuantLib::ForwardFlat;
using QuantLib::Linear;
using QuantLib::LogLinear;
using QuantLib::Cubic;
using QuantLib::Bicubic;
using QuantLib::ConvexMonotone;
using QuantLib::DefaultLogCubic;
using QuantLib::MonotonicLogCubic;
using QuantLib::KrugerLog;
using QuantLib::BSplineInterpolation;

class MonotonicCubic : public Cubic {
  public:
    MonotonicCubic()
    : Cubic(CubicInterpolation::Spline, true,
            CubicInterpolation::SecondDerivative, 0.0,
            CubicInterpolation::SecondDerivative, 0.0) {}
};

class SplineCubic : public Cubic {
  public:
    SplineCubic()
    : Cubic(CubicInterpolation::Spline, false,
            CubicInterpolation::SecondDerivative, 0.0,
            CubicInterpolation::SecondDerivative, 0.0) {}
};

class Kruger : public Cubic {
  public:
    Kruger()
    : Cubic(CubicInterpolation::Kruger) {}
};

class SplineLogCubic : public QuantLib::LogCubic {
  public:
    SplineLogCubic()
    : QuantLib::LogCubic(CubicInterpolation::Spline, false,
                         CubicInterpolation::SecondDerivative, 0.0,
                         CubicInterpolation::SecondDerivative, 0.0) {}
};

// class LogMixedLinearCubic : public QuantLib::LogMixedLinearCubic {
//   public:
//     // We add defaults for all constructor arguments because wrappers for
//     // InterpolatedDiscountCurve and PiecewiseYieldCurve assume that all
//     // interpolators have default constructors.
//     LogMixedLinearCubic(
//         Size n = 0,
//         MixedInterpolation::Behavior behavior = MixedInterpolation::ShareRanges,
//         CubicInterpolation::DerivativeApprox da = CubicInterpolation::Spline,
//         bool monotonic = true)
//     : QuantLib::LogMixedLinearCubic(n, behavior, da, monotonic) {}
// };

// class MixedRateTimeLinearParabolic : public QuantLib::MixedRateTimeLinearParabolic {
//   public:
//     MixedRateTimeLinearParabolic(
//         Size n = 0)
//     : QuantLib::MixedRateTimeLinearParabolic(n) {}
// };

class ParabolicCubic : public QuantLib::Cubic {
  public:
    ParabolicCubic()
    : QuantLib::Cubic(CubicInterpolation::Parabolic, false,
                      CubicInterpolation::SecondDerivative, 0.0,
                      CubicInterpolation::SecondDerivative, 0.0) {}
};

class MonotonicParabolicCubic : public QuantLib::Cubic {
  public:
    MonotonicParabolicCubic()
    : QuantLib::Cubic(CubicInterpolation::Parabolic, true,
                      CubicInterpolation::SecondDerivative, 0.0,
                      CubicInterpolation::SecondDerivative, 0.0) {}
};

class LogParabolicCubic : public QuantLib::LogCubic {
  public:
    LogParabolicCubic()
    : QuantLib::LogCubic(CubicInterpolation::Parabolic, false,
                         CubicInterpolation::SecondDerivative, 0.0,
                         CubicInterpolation::SecondDerivative, 0.0) {}
};

class MonotonicLogParabolicCubic : public QuantLib::LogCubic {
  public:
    MonotonicLogParabolicCubic()
    : QuantLib::LogCubic(CubicInterpolation::Parabolic, true,
                         CubicInterpolation::SecondDerivative, 0.0,
                         CubicInterpolation::SecondDerivative, 0.0) {}
};
%}

%nodefaultctor CubicInterpolation;
struct CubicInterpolation {
    enum DerivativeApprox {
        Spline,
        SplineOM1,
        SplineOM2,
        FourthOrder,
        Parabolic,
        FritschButland,
        Akima,
        Kruger,
        Harmonic,
    };
};

%nodefaultctor MixedInterpolation;
struct MixedInterpolation {
    enum Behavior { ShareRanges, SplitRanges };
};

struct BackwardFlat {};
struct ForwardFlat {};
struct Linear {};
struct LogLinear {};
struct Cubic {};
struct Bicubic {};
struct MonotonicCubic {};
struct DefaultLogCubic {};
struct MonotonicLogCubic {};
struct SplineCubic {};
struct SplineLogCubic {};
struct Kruger {};
struct KrugerLog {};
struct ConvexMonotone {
    ConvexMonotone(Real quadraticity = 0.3,
                   Real monotonicity = 0.7,
                   bool forcePositive = true);
};
struct ParabolicCubic {};
struct MonotonicParabolicCubic {};
struct LogParabolicCubic {};
struct MonotonicLogParabolicCubic {};

// struct LogMixedLinearCubic {
//     #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
//     %feature("kwargs") LogMixedLinearCubic;
//     #endif
//     LogMixedLinearCubic(
//         Size n = 0,
//         MixedInterpolation::Behavior behavior = MixedInterpolation::ShareRanges,
//         CubicInterpolation::DerivativeApprox da = CubicInterpolation::Spline,
//         bool monotonic = true);
// };


%{
using QuantLib::RichardsonExtrapolation;
%}

class RichardsonExtrapolation {
  public:
    Real operator()(Real t=2.0) const;
    Real operator()(Real t, Real s) const;
    
#if defined(SWIGPYTHON)
    %extend {
        RichardsonExtrapolation(
            PyObject* fct, Real delta_h, Real n = Null<Real>()) {
        
            UnaryFunction f(fct);
            return new RichardsonExtrapolation(f, delta_h, n); 
        }
    }
#elif defined(SWIGJAVA) || defined(SWIGCSHARP)
    %extend {
        RichardsonExtrapolation(
            UnaryFunctionDelegate* fct, Real delta_h, Real n = Null<Real>()) {
        
            UnaryFunction f(fct);
            return new RichardsonExtrapolation(f, delta_h, n); 
        }
    }
#else
  private:
    RichardsonExtrapolation();
#endif
};


%{
class SafeConvexMonotoneInterpolation {
  public:
    SafeConvexMonotoneInterpolation(const Array& x, const Array& y,
                                    Real quadraticity = 0.3,
                                    Real monotonicity = 0.7,
                                    bool forcePositive = true)
    : x_(x), y_(y), f_(x_.begin(), x_.end(), y_.begin(),
                       quadraticity, monotonicity, forcePositive) {}
    Real operator()(Real x, bool allowExtrapolation=false) {
        return f_(x, allowExtrapolation);
    }
    Array x_, y_;
    QuantLib::ConvexMonotoneInterpolation<Array::const_iterator, Array::const_iterator> f_;
};
%}


%{
using QuantLib::ChebyshevInterpolation;
%}

class ChebyshevInterpolation {
    #if defined(SWIGCSHARP)
    %rename(call) operator();
    #endif

  public:
    enum PointsType {FirstKind, SecondKind};
    ChebyshevInterpolation(const Array& f, PointsType pointsType = SecondKind);
#if defined(SWIGPYTHON)
    %extend {
        ChebyshevInterpolation(
            Size n, PyObject* fct, PointsType pointsType = SecondKind) {
        
            UnaryFunction f(fct);
            return new ChebyshevInterpolation(n, f, pointsType); 
        }
    }
#elif defined(SWIGJAVA) || defined(SWIGCSHARP)
    %extend {
        ChebyshevInterpolation(
            Size n, UnaryFunctionDelegate* fct, PointsType pointsType = SecondKind) {
        
            UnaryFunction f(fct);
            return new ChebyshevInterpolation(n, f, pointsType); 
        }
    }
#endif
    
    Real operator()(Real z, bool allowExtrapolation=false) const;
    static Array nodes(Size n, PointsType pointsType);
};


%rename(ConvexMonotoneInterpolation) SafeConvexMonotoneInterpolation;
class SafeConvexMonotoneInterpolation {
    #if defined(SWIGCSHARP)
    %rename(call) operator();
    #endif
  public:
    SafeConvexMonotoneInterpolation(const Array& x, const Array& y,
                                    Real quadraticity = 0.3,
                                    Real monotonicity = 0.7,
                                    bool forcePositive = true);
    Real operator()(Real x, bool allowExtrapolation=false);
};

// Include debug inspection interface
%include bspline_debug.i

// Add debug extensions to B-spline classes
%extend BSplineStructure {
    void printDebugInfo() {
        QuantLib::BSplineDebugInspector::printBasicInfo(*$self);
    }
}

#endif
