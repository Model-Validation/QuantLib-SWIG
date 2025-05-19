
/*
 Copyright (C) 2005, 2006 StatPro Italia srl

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

#ifndef quantlib_forward_curve_i
#define quantlib_forward_curve_i

%include termstructures.i
%include interpolation.i

%{
using QuantLib::InterpolatedForwardCurve;
using QuantLib::InterpolatedTermForwardCurve;
using QuantLib::InterpolatedInstantaneousForwardCurve;
%}

%shared_ptr(InterpolatedForwardCurve<BackwardFlat>);
%shared_ptr(InterpolatedForwardCurve<ForwardFlat>);
%shared_ptr(InterpolatedForwardCurve<Linear>);
%shared_ptr(InterpolatedForwardCurve<Cubic>);
%shared_ptr(InterpolatedForwardCurve<SplineCubic>);
%shared_ptr(InterpolatedForwardCurve<ParabolicCubic>);
%shared_ptr(InterpolatedForwardCurve<MonotonicParabolicCubic>);
%shared_ptr(InterpolatedTermForwardCurve<Linear>);
%shared_ptr(InterpolatedInstantaneousForwardCurve<Linear>);
            
template <class Interpolator>
class InterpolatedForwardCurve : public YieldTermStructure {
  public:
    InterpolatedForwardCurve(const std::vector<Date>& dates,
                             const std::vector<Rate>& forwards,
                             const DayCounter& dayCounter,
                             const Calendar& calendar = Calendar(),
                             const Interpolator& i = Interpolator());
    const std::vector<Date>& dates() const;
    const std::vector<Rate>& forwards() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Rate> > nodes() const;
    #endif
};

%template(ForwardCurve) InterpolatedForwardCurve<BackwardFlat>;
%template(ForwardFlatForwardCurve) InterpolatedForwardCurve<ForwardFlat>;
%template(LinearForwardCurve) InterpolatedForwardCurve<Linear>;
%template(CubicForwardCurve) InterpolatedForwardCurve<Cubic>;
%template(NaturalCubicForwardCurve) InterpolatedForwardCurve<SplineCubic>;
%template(ParabolicCubicForwardCurve) InterpolatedForwardCurve<ParabolicCubic>;
%template(MonotonicParabolicCubicForwardCurve) InterpolatedForwardCurve<MonotonicParabolicCubic>;

template <class Interpolator>
class InterpolatedTermForwardCurve : public YieldTermStructure {
  public:
    InterpolatedTermForwardCurve(const std::vector<Date>& dates,
                             const std::vector<Rate>& forwards,
                             const ext::shared_ptr<IborIndex>& index,
                             const Compounding comp = Simple,
                             const Frequency freq = NoFrequency,
                             const Interpolator& i = Interpolator());
    const std::vector<Date>& dates() const;
    const std::vector<Rate>& forwards() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Rate> > nodes() const;
    #endif
};

%template(TermForwardCurveLinear) InterpolatedTermForwardCurve<Linear>;

template <class Interpolator>
class InterpolatedInstantaneousForwardCurve : public YieldTermStructure {
  public:
    InterpolatedInstantaneousForwardCurve(const std::vector<Date>& dates,
                             const std::vector<Rate>& forwards,
                             const DayCounter& dayCounter,
                             const Calendar& calendar = Calendar(),
                             const Interpolator& i = Interpolator());
    const std::vector<Date>& dates() const;
    const std::vector<Rate>& forwards() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Rate> > nodes() const;
    #endif
};

%template(InstantaneousForwardCurveLinear) InterpolatedInstantaneousForwardCurve<Linear>;

#endif
