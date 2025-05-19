/*
 Copyright (C) 2008, 2009 StatPro Italia srl
 Copyright (C) 2018 Matthias Lungwitz

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

#ifndef quantlib_default_probability_structures_i
#define quantlib_default_probability_structures_i

%include common.i
%include types.i
%include date.i
%include calendars.i
%include daycounters.i
%include scheduler.i
%include observer.i
%include marketelements.i
%include interpolation.i
%include termstructures.i
%include piecewiseyieldcurve.i
%include bonds.i

%{
using QuantLib::DefaultProbabilityTermStructure;
using QuantLib::BSplineModel;
%}

%shared_ptr(DefaultProbabilityTermStructure);
class DefaultProbabilityTermStructure : public TermStructure {
  private:
    DefaultProbabilityTermStructure();
  public:
    Probability defaultProbability(const Date&, bool extrapolate = false);
    Probability defaultProbability(Time, bool extrapolate = false);
    Probability defaultProbability(const Date&, const Date&,
                                   bool extrapolate = false);
    Probability defaultProbability(Time, Time, bool extrapolate = false);

    Probability survivalProbability(const Date&, bool extrapolate = false);
    Probability survivalProbability(Time, bool extrapolate = false);

    Real defaultDensity(const Date&, bool extrapolate = false);
    Real defaultDensity(Time, bool extrapolate = false);

    Real hazardRate(const Date&, bool extrapolate = false);
    Real hazardRate(Time, bool extrapolate = false);
};


%template(DefaultProbabilityTermStructureHandle)
Handle<DefaultProbabilityTermStructure>;
%template(RelinkableDefaultProbabilityTermStructureHandle)
RelinkableHandle<DefaultProbabilityTermStructure>;


// concrete curves


// flat forward curve

%{
using QuantLib::FlatHazardRate;
%}

%shared_ptr(FlatHazardRate);
class FlatHazardRate : public DefaultProbabilityTermStructure {
  public:
    FlatHazardRate(Integer settlementDays,
                   const Calendar& calendar,
                   const Handle<Quote>& hazardRate,
                   const DayCounter& dayCounter);
    FlatHazardRate(const Date& todaysDate,
                   const Handle<Quote>& hazardRate,
                   const DayCounter& dayCounter);
};


%{
using QuantLib::InterpolatedHazardRateCurve;
%}

// add other instantiations both here and below the class
%shared_ptr(InterpolatedHazardRateCurve<BackwardFlat>);
%shared_ptr(InterpolatedHazardRateCurve<BSplineModel>);

template <class Interpolator>
class InterpolatedHazardRateCurve : public DefaultProbabilityTermStructure {
  public:
    InterpolatedHazardRateCurve(const std::vector<Date>& dates,
                                const std::vector<Real>& hazardRates,
                                const DayCounter& dayCounter,
                                const Calendar& calendar = Calendar(),
                                const Interpolator& i = Interpolator());
    const std::vector<Date>& dates() const;
    const std::vector<Real>& hazardRates() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Real> > nodes() const;
    #endif
};

%template(HazardRateCurve) InterpolatedHazardRateCurve<BackwardFlat>;
%template(BSplineHazardRateCurve) InterpolatedHazardRateCurve<BSplineModel>;

%{
using QuantLib::InterpolatedDefaultDensityCurve;
%}

// add other instantiations both here and below the class
%shared_ptr(InterpolatedDefaultDensityCurve<Linear>);
%shared_ptr(InterpolatedDefaultDensityCurve<BSplineModel>);

template <class Interpolator>
class InterpolatedDefaultDensityCurve : public DefaultProbabilityTermStructure {
  public:
    InterpolatedDefaultDensityCurve(const std::vector<Date>& dates,
                                    const std::vector<Real>& densities,
                                    const DayCounter& dayCounter,
                                    const Calendar& calendar = Calendar(),
                                    const Interpolator& i = Interpolator());
    const std::vector<Date>& dates() const;
    const std::vector<Real>& defaultDensities() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Real> > nodes() const;
    #endif
};

%template(DefaultDensityCurve) InterpolatedDefaultDensityCurve<Linear>;
%template(BSplineDefaultDensityCurve) InterpolatedDefaultDensityCurve<BSplineModel>;


%{
using QuantLib::InterpolatedSurvivalProbabilityCurve;
%}

// add other instantiations both here and below the class
%shared_ptr(InterpolatedSurvivalProbabilityCurve<Linear>);
%shared_ptr(InterpolatedSurvivalProbabilityCurve<LogLinear>);
%shared_ptr(InterpolatedSurvivalProbabilityCurve<BSplineModel>);

template <class Interpolator>
class InterpolatedSurvivalProbabilityCurve : public DefaultProbabilityTermStructure {
  public:
    InterpolatedSurvivalProbabilityCurve(const std::vector<Date>& dates,
                                         const std::vector<Probability>& probabilities,
                                         const DayCounter& dayCounter,
                                         const Calendar& calendar = Calendar(),
                                         const Interpolator& i = Interpolator());
    const std::vector<Date>& dates() const;
    const std::vector<Probability>& survivalProbabilities() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Real> > nodes() const;
    #endif
};

%template(SurvivalProbabilityCurve) InterpolatedSurvivalProbabilityCurve<Linear>;
%template(LogSurvivalProbabilityCurve) InterpolatedSurvivalProbabilityCurve<LogLinear>;
%template(BSplineProbabilityCurve) InterpolatedSurvivalProbabilityCurve<BSplineModel>;

%{
using QuantLib::DefaultProbabilityHelper;
using QuantLib::SpreadCdsHelper;
using QuantLib::UpfrontCdsHelper;
%}

// rate helpers for curve bootstrapping

%shared_ptr(DefaultProbabilityHelper)
class DefaultProbabilityHelper : public Observable {
  public:
    Handle<Quote> quote() const;
    Date latestDate() const;
    Date earliestDate() const;
    Date maturityDate() const;
    Date latestRelevantDate() const;
    Date pillarDate() const;
    Real impliedQuote() const;
    Real quoteError() const;
  private:
    DefaultProbabilityHelper();
};

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<DefaultProbabilityHelper> )
#endif
namespace std {
    %template(DefaultProbabilityHelperVector)
    vector<ext::shared_ptr<DefaultProbabilityHelper> >;
}


%shared_ptr(SpreadCdsHelper)
class SpreadCdsHelper : public DefaultProbabilityHelper {
  public:
    SpreadCdsHelper(const Handle<Quote>& runningSpread,
                    const Period& tenor,
                    Integer settlementDays,
                    const Calendar& calendar,
                    Frequency frequency,
                    BusinessDayConvention paymentConvention,
                    DateGeneration::Rule rule,
                    const DayCounter& dayCounter,
                    Real recoveryRate,
                    const Handle<YieldTermStructure>& discountCurve,
                    bool settlesAccrual = true,
                    CreditDefaultSwap::ProtectionPaymentTime protectionPaymentTime =
                        CreditDefaultSwap::ProtectionPaymentTime::atDefault,
                    const Date& startDate = Date(),
                    const DayCounter& lastPeriodDayCounter = DayCounter(),
                    bool rebatesAccrual = true,
                    CreditDefaultSwap::PricingModel model = CreditDefaultSwap::Midpoint);

    SpreadCdsHelper(Rate runningSpread,
                    const Period& tenor,
                    Integer settlementDays,
                    const Calendar& calendar,
                    Frequency frequency,
                    BusinessDayConvention paymentConvention,
                    DateGeneration::Rule rule,
                    const DayCounter& dayCounter,
                    Real recoveryRate,
                    const Handle<YieldTermStructure>& discountCurve,
                    bool settlesAccrual = true,
                    CreditDefaultSwap::ProtectionPaymentTime protectionPaymentTime =
                        CreditDefaultSwap::ProtectionPaymentTime::atDefault,
                    const Date& startDate = Date(),
                    const DayCounter& lastPeriodDayCounter = DayCounter(),
                    bool rebatesAccrual = true,
                    CreditDefaultSwap::PricingModel model = CreditDefaultSwap::Midpoint);

    ext::shared_ptr<CreditDefaultSwap> swap() const;
    void setTermStructure(DefaultProbabilityTermStructure*);
};


%shared_ptr(UpfrontCdsHelper)
class UpfrontCdsHelper : public DefaultProbabilityHelper {
  public:
    /*! \note the upfront must be quoted in fractional units. */
    UpfrontCdsHelper(const Handle<Quote>& upfront,
                      Rate runningSpread,
                      const Period& tenor,
                      Integer settlementDays,
                      const Calendar& calendar,
                      Frequency frequency,
                      BusinessDayConvention paymentConvention,
                      DateGeneration::Rule rule,
                      const DayCounter& dayCounter,
                      Real recoveryRate,
                      const Handle<YieldTermStructure>& discountCurve,
                      Natural upfrontSettlementDays = 3,
                      bool settlesAccrual = true,
                      CreditDefaultSwap::ProtectionPaymentTime protectionPaymentTime =
                          CreditDefaultSwap::ProtectionPaymentTime::atDefault,
                      const Date& startDate = Date(),
                      const DayCounter& lastPeriodDayCounter = DayCounter(),
                      bool rebatesAccrual = true,
                      CreditDefaultSwap::PricingModel model = CreditDefaultSwap::Midpoint);

    /*! \note the upfront must be quoted in fractional units. */
    UpfrontCdsHelper(Rate upfront,
                      Rate runningSpread,
                      const Period& tenor,
                      Integer settlementDays,
                      const Calendar& calendar,
                      Frequency frequency,
                      BusinessDayConvention paymentConvention,
                      DateGeneration::Rule rule,
                      const DayCounter& dayCounter,
                      Real recoveryRate,
                      const Handle<YieldTermStructure>& discountCurve,
                      Natural upfrontSettlementDays = 3,
                      bool settlesAccrual = true,
                      CreditDefaultSwap::ProtectionPaymentTime protectionPaymentTime =
                          CreditDefaultSwap::ProtectionPaymentTime::atDefault,
                      const Date& startDate = Date(),
                      const DayCounter& lastPeriodDayCounter = DayCounter(),
                      bool rebatesAccrual = true,
                      CreditDefaultSwap::PricingModel model = CreditDefaultSwap::Midpoint);

    ext::shared_ptr<CreditDefaultSwap> swap() const;
    void setTermStructure(DefaultProbabilityTermStructure*);
};



// bootstrap traits

%{
using QuantLib::HazardRate;
using QuantLib::DefaultDensity;
using QuantLib::SurvivalProbability;
%}

struct HazardRate {};
struct DefaultDensity {};

// curve

%{
using QuantLib::PiecewiseDefaultCurve;
%}

/* We have to resort to a macro, because the R implementation of shared_ptr
   can't take class templates with two or more template arguments. */

%define export_piecewise_default_curve(Name,Traits,Interpolator)

%{
typedef PiecewiseDefaultCurve<Traits, Interpolator> Name;
%}

%shared_ptr(Name);
class Name : public DefaultProbabilityTermStructure {
  public:
    %extend {
        Name(const Date& referenceDate,
             const std::vector<ext::shared_ptr<DefaultProbabilityHelper> >& instruments,
             const DayCounter& dayCounter,
             const Interpolator& i = Interpolator(),
             const _IterativeBootstrap& b = _IterativeBootstrap()) {
            return new Name(referenceDate, instruments, dayCounter,
                            i, Name::bootstrap_type(b.accuracy, b.minValue, b.maxValue));
        }
        Name(Integer settlementDays, const Calendar& calendar,
             const std::vector<ext::shared_ptr<DefaultProbabilityHelper> >& instruments,
             const DayCounter& dayCounter,
             const Interpolator& i = Interpolator(),
             const _IterativeBootstrap& b = _IterativeBootstrap()) {
            return new Name(settlementDays, calendar, instruments, dayCounter,
                            i, Name::bootstrap_type(b.accuracy, b.minValue, b.maxValue));
        }
        Name(const Date& referenceDate,
             const std::vector<ext::shared_ptr<DefaultProbabilityHelper> >& instruments,
             const DayCounter& dayCounter,
             const _IterativeBootstrap& b) {
            return new Name(referenceDate, instruments, dayCounter,
                            Interpolator(), Name::bootstrap_type(b.accuracy, b.minValue, b.maxValue));
        }
        Name(Integer settlementDays, const Calendar& calendar,
             const std::vector<ext::shared_ptr<DefaultProbabilityHelper> >& instruments,
             const DayCounter& dayCounter,
             const _IterativeBootstrap& b) {
            return new Name(settlementDays, calendar, instruments, dayCounter,
                            Interpolator(), Name::bootstrap_type(b.accuracy, b.minValue, b.maxValue));
        }
    }
    const std::vector<Date>& dates() const;
    const std::vector<Time>& times() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Real> > nodes() const;
    #endif
};

%enddef

// add other instantiations if you need them
export_piecewise_default_curve(PiecewiseFlatHazardRate, HazardRate, BackwardFlat);
export_piecewise_default_curve(PiecewiseForwardFlatHazardRate, HazardRate, ForwardFlat);
export_piecewise_default_curve(PiecewiseLogLinearSurvivalProbability, SurvivalProbability, LogLinear);


%define export_iterative_bspline_default_curve(Name, Traits)

%{
typedef PiecewiseDefaultCurve<Traits, BSplineModel> Name;
%}

%shared_ptr(Name);
class Name: public DefaultProbabilityTermStructure {
  public:
    %extend {
        Name(const Date& referenceDate,
                        const std::vector<ext::shared_ptr<DefaultProbabilityHelper> >& instruments,
                        const DayCounter& dayCounter,
                        const BSplineModel& bsplineModel,
                        const _IterativeBootstrap& b = _IterativeBootstrap()) {
                return new Name(
                    referenceDate, instruments, dayCounter, bsplineModel,
                    Name::bootstrap_type(b.accuracy, b.minValue, b.maxValue));
        };

        Name(const Date& referenceDate,
            const std::vector<ext::shared_ptr<DefaultProbabilityHelper> >& instruments,
            const DayCounter& dayCounter,
            const std::vector<Handle<Quote> >& jumps,
            const std::vector<Date>& jumpDates,
            const BSplineModel& bsplineModel,
            const _IterativeBootstrap& b = _IterativeBootstrap()) {
            return new Name(
                    referenceDate, instruments, dayCounter, jumps, jumpDates, bsplineModel,
                    Name::bootstrap_type(b.accuracy));
            }
    }

    const std::vector<Time>& times() const;
    const std::vector<Date>& dates() const;
    const std::vector<Real>& data() const;
    std::vector<std::pair<Date, Real> > nodes() const;

    const Interpolation getInterpolation() const {
        return ext::make_shared<Interpolation>(interpolation_);
    };
};

%enddef

%{
// global bootstrapper
class DefaultAdditionalErrors {
  std::vector<ext::shared_ptr<DefaultProbabilityHelper> > additionalHelpers_;
  public:
    DefaultAdditionalErrors(const std::vector<ext::shared_ptr<DefaultProbabilityHelper> >& additionalHelpers)
    : additionalHelpers_(additionalHelpers) {}
    Array operator()() const {
        Array errors(additionalHelpers_.size() - 2);
        Real a = additionalHelpers_.front()->impliedQuote();
        Real b = additionalHelpers_.back()->impliedQuote();
        for (Size k = 0; k < errors.size(); ++k) {
            errors[k] = (static_cast<Real>(errors.size()-k) * a + static_cast<Real>(1+k) * b) / static_cast<Real>(errors.size()+1)
                - additionalHelpers_.at(1+k)->impliedQuote();
        }
        return errors;
    }
};

struct _DefaultGlobalBootstrap {
    std::vector<ext::shared_ptr<DefaultProbabilityHelper> > additionalHelpers;
    std::vector<Date> additionalDates;
    double accuracy;
    _DefaultGlobalBootstrap(double accuracy = Null<double>())
    : accuracy(accuracy) {}
   _DefaultGlobalBootstrap(const std::vector<ext::shared_ptr<DefaultProbabilityHelper> >& additionalHelpers,
                     const std::vector<Date>& additionalDates,
                     double accuracy = Null<double>())
    : additionalHelpers(additionalHelpers), additionalDates(additionalDates), accuracy(accuracy) {}
};
%}

%rename(DefaultGlobalBootstrap) _DefaultGlobalBootstrap;
struct _DefaultGlobalBootstrap {
    _DefaultGlobalBootstrap(doubleOrNull accuracy = Null<double>());
    _DefaultGlobalBootstrap(const std::vector<ext::shared_ptr<DefaultProbabilityHelper> >& additionalHelpers,
                     const std::vector<Date>& additionalDates,
                     doubleOrNull accuracy = Null<double>());
};


%define export_global_bspline_default_curve(Name, Traits)

%{
typedef PiecewiseDefaultCurve<Traits, BSplineModel, QuantLib::GlobalBootstrap> Name;
%}

%shared_ptr(Name);
class Name: public DefaultProbabilityTermStructure {
  public:
    %extend {
      Name(const Date& referenceDate,
        const std::vector<ext::shared_ptr<DefaultProbabilityHelper> >& instruments,
        const DayCounter& dayCounter,
        const BSplineModel& bsplineModel,
        const _DefaultGlobalBootstrap& b
      ) {
        if (b.additionalHelpers.empty()) {
          return new Name(
            referenceDate, instruments, dayCounter, bsplineModel,
            Name::bootstrap_type(b.accuracy)
          );
        } else {
          return new Name(
            referenceDate, instruments, dayCounter, bsplineModel,
            Name::bootstrap_type(b.additionalHelpers,
              AdditionalDates(b.additionalDates),
              DefaultAdditionalErrors(b.additionalHelpers),
              b.accuracy
            )
          );
        }
      };

      Name(const Date& referenceDate,
        const std::vector<ext::shared_ptr<DefaultProbabilityHelper> >& instruments,
        const DayCounter& dayCounter,
        const std::vector<Handle<Quote> >& jumps,
        const std::vector<Date>& jumpDates,
        const BSplineModel& bsplineModel,
        const _DefaultGlobalBootstrap& b
      ) {
        if (b.additionalHelpers.empty()) {
          return new Name(
            referenceDate, instruments, dayCounter, bsplineModel,
            Name::bootstrap_type(b.accuracy)
          );
        } else {
          return new Name(
            referenceDate, instruments, dayCounter, bsplineModel,
            Name::bootstrap_type(b.additionalHelpers,
              AdditionalDates(b.additionalDates),
              DefaultAdditionalErrors(b.additionalHelpers),
              b.accuracy
            )
          );
        }
      };

    }

    const std::vector<Time>& times() const;
    const std::vector<Date>& dates() const;
    const std::vector<Real>& data() const;
    std::vector<std::pair<Date, Real> > nodes() const;

    const Interpolation getInterpolation() const {
        return ext::make_shared<Interpolation>(interpolation_);
    };
};

%enddef

export_iterative_bspline_default_curve(PiecewiseBSplineHazardCurve, HazardRate);
export_iterative_bspline_default_curve(PiecewiseBSplineSurvivalProbabilityCurve, SurvivalProbability);
export_iterative_bspline_default_curve(PiecewiseBSplineDefaultDensityCurve, DefaultDensity);

export_global_bspline_default_curve(GlobalPiecewiseBSplineHazardCurve, HazardRate);
export_global_bspline_default_curve(GlobalPiecewiseBSplineSurvivalProbabilityCurve, SurvivalProbability);
export_global_bspline_default_curve(GlobalPiecewiseBSplineDefaultDensityCurve, DefaultDensity);

// bond engine based on default probability

%{
using QuantLib::RiskyBondEngine;
%}

%shared_ptr(RiskyBondEngine)
class RiskyBondEngine : public PricingEngine {
  public:
    RiskyBondEngine(const Handle<DefaultProbabilityTermStructure>& defaultCurve,
                    Real recoveryRate,
                    const Handle<YieldTermStructure>& riskFreeCurve);
};


#endif
