/*
 Copyright (C) 2024 QuantLib Contributors

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

#ifndef quantlib_experimental_isda_i
#define quantlib_experimental_isda_i

%include instruments.i
%include creditdefaultswap.i
%include termstructures.i
%include date.i

#ifdef QL_ENABLE_ISDA_CDS

%{
using QuantLib::IsdaCdsInterface;
%}

#if defined(SWIGPYTHON)
%feature("docstring") IsdaCdsInterface "
Interface to ISDA CDS Standard Model C library.

This class provides a bridge between QuantLib objects and the
ISDA CDS Standard Model C library, allowing for ISDA-compliant
CDS pricing and analytics.

Example usage:
    # Convert QuantLib curves to ISDA format
    isda_yield_curve = IsdaCdsInterface.convertYieldCurve(yield_curve, base_date)
    
    # Price CDS using ISDA methodology
    price = IsdaCdsInterface.priceCds(cds, yield_curve, default_curve, recovery_rate, today)
    
    # Calculate par spread
    interface = IsdaCdsInterface()
    par_spread = interface.calculateParSpread(schedule, yield_curve, default_curve, 
                                            recovery_rate, today, protection_start)
    compila
    # Bootstrap default curve from CDS quotes
    default_curve = interface.bootstrapDefaultCurve(yield_curve, cds_quotes, tenors, 
                                                  recovery_rate, today)
    
    # Remember to free ISDA curves when done
    IsdaCdsInterface.freeTCurve(isda_yield_curve)
";
#endif

#if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
%feature("kwargs") IsdaCdsInterface::priceCds;
%feature("kwargs") IsdaCdsInterface::calculateParSpread;
%feature("kwargs") IsdaCdsInterface::bootstrapDefaultCurve;
#endif

// Include ISDA headers for SWIG
%{
#include <isda/tcurve.h>
#include <isda/bastypes.h>
#include <isda/cdate.h>
%}

// Python-friendly TCurve representation
%pythoncode %{
class TCurveData:
    """
    Python-friendly representation of ISDA TCurve structure.
    
    This class provides read-only access to all TCurve fields for inspection.
    All data is copied from the original C structure, so no memory management is needed.
    """
    
    def __init__(self, num_items, dates, rates, base_date, basis, day_count_conv):
        self._num_items = num_items
        self._dates = dates
        self._rates = rates
        self._base_date = base_date
        self._basis = basis
        self._day_count_conv = day_count_conv
    
    @property
    def num_items(self):
        """Number of rate points in the curve."""
        return self._num_items
    
    @property
    def dates(self):
        """List of dates (as TDate long integers) for rate points."""
        return self._dates.copy()  # Return copy to prevent modification
    
    @property
    def rates(self):
        """List of rates/discount factors for rate points."""
        return self._rates.copy()  # Return copy to prevent modification
    
    @property
    def base_date(self):
        """Base/discount date for the curve (as TDate long integer)."""
        return self._base_date
    
    @property
    def basis(self):
        """Compounding basis (number of periods per year)."""
        return self._basis
    
    @property
    def day_count_conv(self):
        """Day count convention identifier."""
        return self._day_count_conv
    
    def get_rate_point(self, index):
        """Get rate point at given index as (date, rate) tuple."""
        if 0 <= index < self._num_items:
            return (self._dates[index], self._rates[index])
        else:
            raise IndexError(f"Index {index} out of range [0, {self._num_items})")
    
    def get_date_range(self):
        """Get date range as (min_date, max_date) tuple."""
        if self._num_items > 0:
            return (min(self._dates), max(self._dates))
        else:
            return (None, None)
    
    def get_rate_range(self):
        """Get rate range as (min_rate, max_rate) tuple."""
        if self._num_items > 0:
            return (min(self._rates), max(self._rates))
        else:
            return (None, None)
    
    def day_count_conv_name(self):
        """Get human-readable name for day count convention."""
        day_count_names = {
            1: "ACT/365",
            2: "ACT/365F", 
            3: "ACT/360",
            4: "30/360",
            5: "30E/360",
            8: "EFFECTIVE_RATE"
        }
        return day_count_names.get(self._day_count_conv, f"Unknown({self._day_count_conv})")
    
    def basis_name(self):
        """Get human-readable name for compounding basis."""
        if self._basis == 0:
            return "Simple"
        elif self._basis == 1:
            return "Annual"
        elif self._basis == 5000:
            return "Continuous"
        elif self._basis == -2:
            return "Discount Factor"
        else:
            return f"Custom({self._basis})"
    
    def __str__(self):
        return (f"TCurveData(num_items={self._num_items}, "
                f"base_date={self._base_date}, "
                f"basis={self.basis_name()}, "
                f"day_count={self.day_count_conv_name()})")
    
    def __repr__(self):
        return self.__str__()
    
    def summary(self):
        """Get a summary of the curve data."""
        if self._num_items == 0:
            return "Empty curve"
        
        date_range = self.get_date_range()
        rate_range = self.get_rate_range()
        
        return f"""TCurve Summary:
  Number of points: {self._num_items}
  Base date: {self._base_date}
  Date range: {date_range[0]} to {date_range[1]}
  Rate range: {rate_range[0]:.6f} to {rate_range[1]:.6f}
  Basis: {self.basis_name()} ({self._basis})
  Day count: {self.day_count_conv_name()} ({self._day_count_conv})
"""

# Helper function to convert TDate to Python date
def tdate_to_python_date(tdate):
    """Convert ISDA TDate (long) to Python date object."""
    # ISDA TDate is days since Jan 1, 1900
    # Python date needs days since Jan 1, 1900 but with different epoch
    import datetime
    try:
        # TDate 0 corresponds to Jan 1, 1900 in ISDA
        base_date = datetime.date(1900, 1, 1)
        return base_date + datetime.timedelta(days=int(tdate))
    except (ValueError, OverflowError):
        return None

def python_date_to_tdate(python_date):
    """Convert Python date object to ISDA TDate (long)."""
    import datetime
    try:
        base_date = datetime.date(1900, 1, 1)
        delta = python_date - base_date
        return int(delta.days)
    except (ValueError, OverflowError):
        return None
%}

class IsdaCdsInterface {
  public:
    // Default constructor
    IsdaCdsInterface();
    
    // Destructor
    ~IsdaCdsInterface();

    // Static method for CDS pricing
    static Real priceCds(
        const CreditDefaultSwap& cds,
        const Handle<YieldTermStructure>& yieldCurve,
        const Handle<DefaultProbabilityTermStructure>& defaultCurve,
        Real recoveryRate,
        const Date& today);

    // Instance methods
    Real calculateParSpread(
        const Schedule& schedule,
        const Handle<YieldTermStructure>& yieldCurve,
        const Handle<DefaultProbabilityTermStructure>& defaultCurve,
        Real recoveryRate,
        const Date& today,
        const Date& protectionStart) const;

    // Date conversion utilities
    static long dateToTDate(const Date& date);
    static Date tDateToDate(long tDate);
};

#if defined(SWIGPYTHON)
%pythoncode %{
# Python convenience functions for ISDA CDS interface

def create_isda_yield_curve(yield_curve, base_date):
    """
    Create ISDA yield curve data for inspection.
    
    Note: This function is temporarily disabled due to compilation issues.
    Use priceCds directly for now.
    """
    raise NotImplementedError("Curve conversion temporarily disabled - use priceCds directly")

def create_isda_default_curve(default_curve, base_date, recovery_rate=0.4):
    """
    Create ISDA default curve data for inspection.
    
    Note: This function is temporarily disabled due to compilation issues.
    Use priceCds directly for now.
    """
    raise NotImplementedError("Curve conversion temporarily disabled - use priceCds directly")

def price_cds_isda(cds, yield_curve, default_curve, recovery_rate, today):
    """
    Price CDS using ISDA methodology.
    
    Args:
        cds: QuantLib CreditDefaultSwap instrument
        yield_curve: QuantLib yield term structure handle
        default_curve: QuantLib default probability term structure handle
        recovery_rate: Recovery rate assumption
        today: Valuation date
        
    Returns:
        CDS price using ISDA methodology
    """
    return IsdaCdsInterface.priceCds(cds, yield_curve, default_curve, recovery_rate, today)

class IsdaCdsHelper:
    """
    Helper class for ISDA CDS operations.
    
    Note: Curve conversion is temporarily disabled. Use pricing functions directly.
    """
    
    def __init__(self):
        self.interface = IsdaCdsInterface()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        pass
    
    def convert_yield_curve(self, yield_curve, base_date):
        """Convert yield curve - temporarily disabled."""
        raise NotImplementedError("Curve conversion temporarily disabled - use priceCds directly")
    
    def convert_default_curve(self, default_curve, base_date, recovery_rate=0.4):
        """Convert default curve - temporarily disabled."""
        raise NotImplementedError("Curve conversion temporarily disabled - use priceCds directly")
    
    def bootstrap_default_curve(self, yield_curve, cds_quotes, tenors, recovery_rate, today):
        """Bootstrap default curve - temporarily disabled."""
        raise NotImplementedError("Bootstrap temporarily disabled - use priceCds directly")
    
    def calculate_par_spread(self, schedule, yield_curve, default_curve, recovery_rate, today, protection_start):
        """Calculate par spread using ISDA methodology."""
        return self.interface.calculateParSpread(schedule, yield_curve, default_curve, 
                                               recovery_rate, today, protection_start)
    
    def price_cds(self, cds, yield_curve, default_curve, recovery_rate, today):
        """Price CDS using ISDA methodology."""
        return IsdaCdsInterface.priceCds(cds, yield_curve, default_curve, recovery_rate, today)

%}
#endif

#endif // QL_ENABLE_ISDA_CDS

#endif // quantlib_experimental_isda_i