module MSInstrument

export MSErrorModel,
    signal_std, signal_precision, signal_snr, signal_std,
    detection_loglikelihood

using Distributions, StatsFuns, SpecialFunctions, LogExpFunctions
using Printf: @printf, @sprintf

# default probability that MS measurement is outlier
const DefaultOutlierProb = 1E-4

"""
Mixture of *Normal* and *Cauchy* distributions.

*Normal* distribution handles the noise of the proper *signal*,
while *Cauchy* models the *outliers*.
The `outlierProb` parameter defines the probability of the outlier, i.e.
the *weight* of Cauchy distribution.

Both distributions have the same `center` and `scale` parameters.
"""
struct GaussCauchyMix <: AbstractMixtureModel{Univariate, Continuous, Union{Normal, Cauchy}}
    signal::Normal
    outlier::Cauchy
    mix::Categorical

    signalLogprob::Float64
    outlierLogprob::Float64
    outlierThreshold::Float64 # the point where normal*(1-outlierProb) = cauchy*outlierProb

    function GaussCauchyMix(center::Float64 = 0.0, scale::Float64 = 1.0;
                            outlierProb::Float64 = DefaultOutlierProb)
        signal = Normal(center, scale)
        outlier = Cauchy(center, scale)
        0 < outlierProb < 1 || throw(ArgumentError("outlierProb should be in (0,1) range ($(outlierProb) given)"))
        new(signal, outlier,
            Categorical([1 - outlierProb, outlierProb]),
            log1p(-outlierProb), log(outlierProb),
            sqrt(-1 - 2*lambertw(outlierProb/(outlierProb-1)/sqrt(2*π*ℯ), -1)))
    end
end

outlierProb(d::GaussCauchyMix) = last(probs(d.mix))
signalProb(d::GaussCauchyMix) = first(probs(d.mix))
outlierLogprob(d::GaussCauchyMix) = d.outlierLogprob
signalLogprob(d::GaussCauchyMix) = d.signalLogprob
# where the weighted signal pdf equals outlier pdf
outlierThreshold(d::GaussCauchyMix) = d.outlierThreshold
Distributions.ncomponents(::GaussCauchyMix) = 2
Distributions.probs(d::GaussCauchyMix) = probs(d.mix)
Distributions.component(d::GaussCauchyMix, i::Integer) =
    i == 1 ? d.signal : d.outlier

zvalue(d::GaussCauchyMix, x::Float64) =
    StatsFuns.zval(location(d.signal), scale(d.signal), x)

_stdpdf(d::GaussCauchyMix, z::Float64) =
    signalProb(d) * normpdf(z) +
    outlierProb(d) / (π * (1 + abs2(z)))

Distributions.pdf(d::GaussCauchyMix, x::Float64) =
    _stdpdf(d, zvalue(d, x)) / std(d.signal)

_stdlogpdf(d::GaussCauchyMix, z::Float64) =
    logaddexp(signalLogprob(d) + normlogpdf(z),
              outlierLogprob(d) - log1psq(z) - Distributions.logπ)

Distributions.logpdf(d::GaussCauchyMix, x::Float64) =
    _stdlogpdf(d, zvalue(d, x)) - log(std(d.signal))

"""
Calibrated MS instrument noise/detection model parameters.
"""
mutable struct MSErrorModel{B}
    # log(signal) to zscore transformation
    zScale::Float64
    zShift::Float64
    # distribution of the signal normalized to its std
    distr::GaussCauchyMix

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

    MSErrorModel{B}(zScale::Number, zShift::Number;
                    outlierProb::Float64 = DefaultOutlierProb) where B =
        new{B}(zScale, zShift,
               GaussCauchyMix(outlierProb = outlierProb))

    MSErrorModel(zScale::Number, zShift::Number;
                 logintensityBase::Number = 2,
                 outlierProb::Float64 = DefaultOutlierProb) =
        MSErrorModel{logintensityBase}(zScale, zShift; outlierProb = outlierProb)

    function MSErrorModel(
        zScale::Float64, zShift::Float64,
        zDetectionFactor::Float64, zDetectionIntercept::Float64, detectionMax::Float64,
        sigmaScaleHi::Float64, sigmaScaleLo::Float64, sigmaOffset::Float64,
        sigmaBend::Float64, sigmaSmooth::Float64;
        logintensityBase::Number = 2,
        outlierProb::Float64 = DefaultOutlierProb)
        logintensityBase > 1 || throw(ArgumentError("logintensityBase must be > 1 (found $logintensityBase)"))
        zScale > 0 || throw(ArgumentError("zScale must be positive (found $zScale)"))
        0.0 < outlierProb < 1 || throw(ArgumentError("outlierProb must be in range (0,1) ($outlierProb given)"))
        zDetectionFactor > 0 || throw(ArgumentError("zDetectionFactor must be positive (found $zDetectionFactor)"))
        sigmaScaleHi > 0 || throw(ArgumentError("sigmaScaleHi must be positive (found $sigmaScaleHi)"))
        sigmaScaleLo > 0 || throw(ArgumentError("sigmaScaleLo must be positive (found $sigmaScaleLo)"))
        sigmaSmooth >= 0 || throw(ArgumentError("sigmaSmooth must be non-negative (found $sigmaSmooth)"))
        0.0 < detectionMax <= 1 || throw(ArgumentError("detectionMax should be in range (0,1]"))
        res = new{logintensityBase}(zScale, zShift,
                  GaussCauchyMix(outlierProb = outlierProb),
                  zDetectionFactor, zDetectionIntercept, detectionMax,
                  sigmaScaleHi, sigmaScaleLo, sigmaOffset, sigmaBend, sigmaSmooth)
        sync!(res)
    end
end

MSErrorModel(dict::AbstractDict) =
    MSErrorModel(dict["zScale"], dict["zShift"],
            dict["zDetectionFactor"], dict["zDetectionIntercept"], dict["detectionMax"],
            dict["sigmaScaleHi"], dict["sigmaScaleLo"], dict["sigmaOffset"],
            dict["sigmaBend"], dict["sigmaSmooth"],
            logintensityBase = get(dict, "logintensityBase", MathConstants.ℯ),
            outlierProb = get(dict, "outlierProb", DefaultOutlierProb))

Base.convert(::Type{<:Dict}, params::MSErrorModel) = Dict{String, Any}(
    "logintensityBase" => logintensityBase(params),
    "zScale" => params.zScale, "zShift" => params.zShift,
    "zDetectionFactor" => params.zDetectionFactor,
    "zDetectionIntercept" => params.zDetectionIntercept,
    "sigmaScaleHi" => params.sigmaScaleHi,
    "sigmaScaleLo" => params.sigmaScaleLo,
    "sigmaOffset" => params.sigmaOffset,
    "sigmaBend" => params.sigmaBend,
    "sigmaSmooth" => params.sigmaSmooth,
    "outlierProb" => outlierProb(params))

logintensity_transform(::MSErrorModel{MathConstants.ℯ}) = log
logintensity_transform(::MSErrorModel{2}) = log2
logintensity_transform(::MSErrorModel{10}) = log10
logintensity_transform(::MSErrorModel{B}) where B = x -> log(B, x)

logintensity_transform_label(::MSErrorModel{MathConstants.ℯ}) = "Log"
logintensity_transform_label(::MSErrorModel{2}) = "Log₂"
logintensity_transform_label(::MSErrorModel{10}) = "Log₁₀"
logintensity_transform_label(::MSErrorModel{B}) where B <: Integer = @sprintf "Log[%d]" B
logintensity_transform_label(::MSErrorModel{B}) where B = @sprintf "Log[%.3f]" B

logintensity_invtransform(::MSErrorModel{MathConstants.ℯ}) = exp
logintensity_invtransform(::MSErrorModel{2}) = exp2
logintensity_invtransform(::MSErrorModel{10}) = exp10
logintensity_invtransform(::MSErrorModel{B}) where B = x -> pow(B, x)

logintensity_invtransform_label(::MSErrorModel{MathConstants.ℯ}) = "Exp"
logintensity_invtransform_label(::MSErrorModel{B}) where B =
    B isa Integer ? @sprintf("%d^", B) : @sprintf("%.3f^", B)

function Base.show(io::IO, params::MSErrorModel)
    @printf(io, "Z=(%s(signal)%+.3f)*%.3f\n",
            logintensity_transform_label(params), -params.zShift, params.zScale)
    @printf(io, "P(detect|Z)=%.4f Logit⁻¹(%.4f Z%+.4f)\n",
            params.detectionMax, params.zDetectionFactor, params.zDetectionIntercept)
    @printf(io, "std(signal)=%s(½(%.4f%+.4f)(Z%+.4f)+½(%.4f%+.4f)√((Z%+.4f)²%+.4f)%+.4f) P(outlier)=%.4g\n",
            logintensity_invtransform_label(params),
            params.sigmaScaleHi, params.sigmaScaleLo, -params.sigmaBend,
            params.sigmaScaleHi, -params.sigmaScaleLo, -params.sigmaBend,
            params.sigmaSmooth, params.sigmaOffset, outlierProb(params))
end

logintensityBase(::MSErrorModel{B}) where B = B
outlierProb(params::MSErrorModel) = outlierProb(params.distr)

"""
Create a copy of the MS error model for the given log-intensity `base`.
"""
function change_logintensityBase(params::MSErrorModel, base::Number)
    old_base = logintensityBase(params)
    (old_base == base) && return copy(params) # same model

    k = log(base, old_base)
    return MSErrorModel(
        params.zScale / k, params.zShift * k,
        params.zDetectionFactor, params.zDetectionIntercept, params.detectionMax,
        params.sigmaScaleHi * k, params.sigmaScaleLo * k, params.sigmaOffset * k,
        params.sigmaBend, params.sigmaSmooth * k^2,
        logintensityBase = base, outlierProb = outlierProb(params)
    )
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
    ifelse(isfinite(signal), (logintensity_transform(params)(signal)-params.zShift)*params.zScale, NaN)

signal2zscore(params::MSErrorModel, signal::Number) =
    ifelse(isfinite(signal), (logintensity_transform(params)(signal)-params.zShift)*params.zScale, NaN)

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
        res = logintensity_invtransform(params)(
                logintensity_transform(params)(signal) - zscore_logstd(params, zscore))
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

zscore_std(params::MSErrorModel, z::Float64) =
    logintensity_invtransform(params)(zscore_logstd(params, z))
zscore_precision(params::MSErrorModel, z::Float64) =
    logintensity_invtransform(params)(-zscore_logstd(params, z)) # FIXME should be -2zscore_logstd

signal_std(params::MSErrorModel, signal::Number) =
    ifelse(isfinite(signal), zscore_std(params, signal2zscore(params, signal)), NaN)
signal_logstd(params::MSErrorModel, signal::Number) =
    ifelse(isfinite(signal), zscore_logstd(params, signal2zscore(params, signal)), NaN)

@inline function zscore_loglikelihood(params::MSErrorModel,
            delta::Number, z::Number)
    logstd = zscore_logstd(params, z)
    _normlogpdf(params.distr, delta*exp(-logstd)) - (logstd - params.zShift) * log(logintensityBase(params))
end

@inline function zscore_loglikelihood_and_class(params::MSErrorModel,
                                      delta::Number, z::Number)
    logstd = zscore_logstd(params, z)
    normdelta = delta*logintensity_invtransform(params)(-logstd)
    logpdf(params.distr, normdelta) - (logstd - params.zShift) * log(logintensityBase(params)),
    abs(normdelta) <= outlierThreshold(params.distr)
end

signal_loglikelihood(params::MSErrorModel, signal::Number, expected::Number) =
    zscore_loglikelihood(params, signal-expected, signal2zscore(params, signal))

function detection_loglikelihood(params::MSErrorModel, is_detected::Bool, expected_log::Float64)
    z = muladd(expected_log, params.signalLogDetectionFactor, params.signalLogDetectionIntercept)
    return ( is_detected
        ? -Distributions.log1pexp(-z) #+params.logDetectionMax # invlogit(z)*detMax
        : -Distributions.log1pexp(z) ) #logsumexp( -Distributions.log1pexp(z)+params.logDetectionMax, params.log1mDetectionMax ) ) # invlogit(-z)*detMax+(1-detMax)
end

function rand_missing_log_intensities(params::MSErrorModel,
                                      signals::AbstractVector{Float64};
                                      nbins::Int = 10,
                                      mean_shift::Float64 = 0.0,
                                      std_scale::Float64 = 3.0,
                                      det_factor_scale::Float64 = 10.0)
    log_signals = logintensity_transform(params).(signals)
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
