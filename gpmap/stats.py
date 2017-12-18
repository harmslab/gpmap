import numpy as np
import math

# -----------------------------------------------------------------------
# Unbiased calculations of sample statistics to error statistics
# -----------------------------------------------------------------------


def c4_correction(n_samples):
    """Return the correction scalar for calculating standard deviation from a normal distribution. """
    k1 = round(n_samples / 2.0, 4)
    k2 = round((n_samples - 1) / 2.0, 4)

    # If the number of samples is < 100, calculate a correction scalar.
    if n_samples < 100:

        if k1.is_integer():
            c4 = np.sqrt(2.0 / (math.pi * (2 * k1 - 1))) * ((2**(2 * k1 - 2) *
                                                             math.factorial(k1 - 1)**2) / math.factorial(2 * k1 - 2))
        elif k2.is_integer():
            c4 = np.sqrt(math.pi / k2) * (math.factorial(2 * k2 - 1) /
                                          (2**(2 * k2 - 1) * math.factorial(k2 - 1)**2))
        else:
            raise Exception("""Non-integer value for correction term, c4.""")

    # Else this scalar is 1
    else:
        c4 = 1.0

    return c4


def unbiased_var(x, axis=None):
    """This enforces that the unbias estimate for variance is calculated"""
    # First make sure that the samples are in a numpy array
    return np.std(x, axis=1, ddof=1)


def unbiased_std(x, axis=None):
    """ A correction to numpy's standard deviation calculation.
    Calculate the unbiased estimation of standard deviation, which includes a correction
    factor for sample sizes < 100.
    """
    # First make sure that the samples are in a numpy array
    x = np.array(x)

    # What are the number of samples for this axis of interest
    if axis is None:
        n_samples = x.size
    else:
        n_samples = x.shape[axis]

    # If only 1 sample is given, just return numpy's normal standard deviation
    if n_samples is 1:
        return np.std(x, axis=axis)

    # Calculate the correction scalar
    c4 = c4_correction(n_samples)

    corrected_s = c4 * np.std(x, axis=axis, ddof=1)
    return corrected_s


def unbiased_sterror(x, axis=None):
    """ Unbiased error. """
    # First make sure that the samples are in a numpy array
    x = np.array(x)

    # What are the number of samples for this axis of interest
    if axis is None:
        n_samples = x.size
    else:
        n_samples = x.shape[axis]

    # Calculate biased standard deviation
    std = np.std(x, axis=axis)

    # Correct with c4 scalar
    c4 = c4_correction(n_samples)
    return std * np.sqrt(1 - c4)


# -----------------------------------------------------------------------
# Correction to sample statistics
# -----------------------------------------------------------------------

def corrected_std(var, n_samples=2):
    """Calculate the unbiased standard deviation from a biased standard deviation. """
    _var = np.array(var)
    _std = np.sqrt(_var)

    # If the sample size is 1, no correction applies
    if n_samples > 100:
        c4 = 1
    else:
        c4 = c4_correction(n_samples)
    return _std / c4


def corrected_sterror(var, n_samples=2):
    """Calculate an unbiased standard error from a BIASED standard deviation. """
    _var = np.array(var)

    # If sample size is 1, no correction applies
    if n_samples > 100:
        _std = np.sqrt(_var)
        sterr = _std / np.sqrt(n_samples)
    else:
        _std = np.sqrt(_var)
        sterr = _std / np.sqrt(n_samples)
        #_std = corrected_std(_var, n_samples=n_samples)
        #sterr = _std/np.sqrt(n_samples)

        #c4 = c4_correction(n_samples)
        #sterr = _std * np.sqrt(c4**(-2)-1)
    return sterr
