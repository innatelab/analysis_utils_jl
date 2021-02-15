module MSInstrument

export MSErrorModel,
    signal_std, signal_precision, signal_snr, signal_std,
    detection_likelihood_log

using Distributions
using Printf: @printf

"""
Calibrated MS instrument noise/detection model parameters.
"""
mutable struct MSErrorModel
    # log(signal) to zscore transformation
    zScale::Float64
    zShift::Float64

    # intercept and scale for zscore conversion into the argument of Bernoulli-Logit
    zDetectionFactor::Float64 # real<lower=0>
    zDetectionIntercept::Float64
    #eightDetectionMax::Float64 # real<lower=0,upper=1>

    detectionMax::Float64 # probability to detect signal from an intensive peptide

    # signal heteroscedascity model parameters
    sigmaScaleHi::Float64 # <lower=0>
    sigmaScaleLo::Float64 # <lower=0>
    sigmaOffset::Float64
    sigmaBend::Float64
    sigmaSmooth::Float64 # <lower=0>

    # intercept and scale of bernlogit() conversion for log(signal)
    signalLogDetectionFactor::Float64
    signalLogDetectionIntercept::Float64

    logDetectionMax::Float64 # log(detectionMax)
    log1mDetectionMax::Float64 # log(1-detectionMax)

    MSErrorModel(zScale::Number, zShift::Number) = new(zScale, zShift)

    function MSErrorModel(
        zScale::Float64, zShift::Float64,
        zDetectionFactor::Float64, zDetectionIntercept::Float64, detectionMax::Float64,
        sigmaScaleHi::Float64, sigmaScaleLo::Float64, sigmaOffset::Float64,
        sigmaBend::Float64, sigmaSmooth::Float64)
        zScale > 0 || throw(ArgumentError("zScale must be positive (found $zScale)"))
        zDetectionFactor > 0 || throw(ArgumentError("zDetectionFactor must be positive (found $zDetectionFactor)"))
        sigmaScaleHi > 0 || throw(ArgumentError("sigmaScaleHi must be positive (found $sigmaScaleHi)"))
        sigmaScaleLo > 0 || throw(ArgumentError("sigmaScaleLo must be positive (found $sigmaScaleLo)"))
        sigmaSmooth >= 0 || throw(ArgumentError("sigmaSmooth must be non-negative (found $sigmaSmooth)"))
        0.0 < detectionMax <= 1 || throw(ArgumentError("detectionMax should be in range (0,1]"))
        res = new(zScale, zShift,
                  zDetectionFactor, zDetectionIntercept, detectionMax,
                  sigmaScaleHi, sigmaScaleLo, sigmaOffset, sigmaBend, sigmaSmooth)
        sync!(res)
    end
end

MSErrorModel(dict::AbstractDict) =
    MSErrorModel(dict["zScale"], dict["zShift"],
            dict["zDetectionFactor"], dict["zDetectionIntercept"], dict["detectionMax"],
            dict["sigmaScaleHi"], dict["sigmaScaleLo"], dict["sigmaOffset"],
            dict["sigmaBend"], dict["sigmaSmooth"])

function Base.show(io::IO, instr::MSErrorModel)
    @printf(io, "Z=(log(signal)%+.3f)*%.3f\n", -instr.zShift, instr.zScale)
    @printf(io, "P(detect|Z)=%.4f Logit⁻¹(%.4f Z%+.4f)\n",
            instr.detectionMax, instr.zDetectionFactor, instr.zDetectionIntercept)
    @printf(io, "std(signal)=Exp(½(%.4f%+.4f)(Z%+.4f)+½(%.4f%+.4f)√((Z%+.4f)²%+.4f)%+.4f)\n",
            instr.sigmaScaleHi, instr.sigmaScaleLo, -instr.sigmaBend,
            instr.sigmaScaleHi, -instr.sigmaScaleLo, -instr.sigmaBend,
            instr.sigmaSmooth, instr.sigmaOffset)
end

# update dependent parameters
function sync!(params::MSErrorModel)
    params.signalLogDetectionFactor = params.zScale * params.zDetectionFactor
    params.signalLogDetectionIntercept = params.zDetectionIntercept - params.zShift * params.zScale * params.zDetectionFactor
    params.logDetectionMax = log(params.detectionMax)
    params.log1mDetectionMax = log1p(-params.detectionMax)
    return params
end

# convert signal to zscore
signal2zscore(params::MSErrorModel, signal::Number) =
    ifelse(isfinite(signal), (log(signal)-params.zShift)*params.zScale, NaN)

# inverse of sigma: FIXME should be inverse of var
function signal_precision(params::MSErrorModel, signal::Float64)
    if isfinite(signal)
        res = zscore_precision(params, signal2zscore(params, signal))
        if !isfinite(res)
            throw(BoundsError("$res precision for signal $signal (zscore = $(signal2zscore(params, signal)))"))
        end
        return res
    else
        return NaN
    end
end

function signal_snr(params::MSErrorModel, signal::Float64)
    if isfinite(signal)
        zscore = signal2zscore(params, signal)
        res = exp(log(signal) - zscore_logstd(params, zscore))
        if !isfinite(res)
            throw(BoundsError("$res SNR for signal $signal (zscore = $zscore)"))
        end
        return res
    else
        return NaN
    end
end

function zscore_logstd(params::MSErrorModel, z::Float64)
    zd = z - params.sigmaBend;
    return 0.5 * (params.sigmaScaleHi + params.sigmaScaleLo) * zd +
           0.5 * (params.sigmaScaleHi - params.sigmaScaleLo) * sqrt(zd*zd + params.sigmaSmooth) +
           params.sigmaOffset
end

zscore_std(params::MSErrorModel, z::Float64) = exp(zscore_logstd(params, z))
zscore_precision(params::MSErrorModel, z::Float64) = exp(-zscore_logstd(params, z)) # FIXME should be -2zscore_logstd

signal_std(params::MSErrorModel, signal::Number) =
    ifelse(isfinite(signal), zscore_std(params, signal2zscore(params, signal)), NaN)
signal_logstd(params::MSErrorModel, signal::Number) =
    ifelse(isfinite(signal), zscore_logstd(params, signal2zscore(params, signal)), NaN)

function signal_likelihood_log(signal::Float64, expected::Float64, signalPrecision::Float64)
    err = abs(expected-signal) * signalPrecision
    return -err # -Distributions.logtwo + log(signalPrecision) #this part is relatively expensive, but is not dependent on expected, so ignore
end

signal_likelihood_log(params::MSErrorModel, signal::Float64, expected::Float64) =
    signal_likelihood_log(signal, expected, signal_precision(params, signal))

function detection_likelihood_log(params::MSErrorModel, is_detected::Bool, expected_log::Float64)
    z = muladd(expected_log, params.signalLogDetectionFactor, params.signalLogDetectionIntercept)
    return ( is_detected
        ? -Distributions.log1pexp(-z) #+params.logDetectionMax # invlogit(z)*detMax
        : -Distributions.log1pexp(z) ) #logsumexp( -Distributions.log1pexp(z)+params.logDetectionMax, params.log1mDetectionMax ) ) # invlogit(-z)*detMax+(1-detMax)
end

function rand_missing_log_intensities(params::MSErrorModel,
                                      signals::Vector{Float64};
                                      nbins::Int = 10,
                                      mean_shift::Float64 = 0.0,
                                      std_scale::Float64 = 3.0,
                                      det_factor_scale::Float64 = 10.0)
    log_signals = log(signals)
    log_qsignals = quantile(log_signals, (1:nbins)/nbins)
    signal_bins = map(x -> max(1, findfirst(q -> isless(x, q), log_qsignals)), log_signals)
    nbinsignals = countmap(signal_bins)

    z_signal_prior = Normal()
    function logprob_posterior_undetected(log_signal::Float64, log_intens)
        # llh of the signal that is not detected
        z_signal = (log_intens-log_signal)*params.zScale*std_scale
        z_missed = map(x -> muladd(x, params.signalLogDetectionFactor * det_factor_scale, params.signalLogDetectionIntercept), log_intens)
        return logpdf(z_signal_prior, z_signal) - sum(Distributions.log1pexp, z_missed)
    end
    rand_signals = map(1:length(nbinsignals)) do bin_ix
        nsignals = nbinsignals[bin_ix]
        run(model(x -> logprob_posterior_undetected(log_qsignals[bin_ix], x), init = log_qsignals[bin_ix]),
            SliceSampler(), SerialMC(nsteps=nsignals*2+10, burnin=10, thinning=2)).samples
    end
    res = similar(signals)
    for i in 1:length(rand_signals)
        res[signal_bins .== i] = rand_signals[i]
    end
    return res
end

function Base.rand(params::MSErrorModel,
                   intensities::AbstractArray,
                   signals::AbstractArray;
                   missing_std_scale::Float64 = 3.0,
                   missing_det_factor_scale::Float64 = 10.0
)
    quant_mask = isfinite.(intensities)
    res = similar(intensities)
    if !all(quant_mask)
        res[!quant_mask] = rand_missing_log_intensities(params, signals[!quant_mask],
                                                        std_scale=missing_std_scale,
                                                        det_factor_scale=missing_det_factor_scale) |> exp
    end
    if any(quant_mask)
        res[quant_mask] = max.(intensities[quant_mask] +
            randn(sum(quant_mask)).*signal_std.(params, intensities[quant_mask]),
            res[quant_mask])
    end
    return res
end

end
