module MSInstrumentCalibration

using StatsBase, Distributions, DataFrames, BlackBoxOptim, Printf, JSON, JLD2

const MSInstrument = Main.MSInstrument
const BBOUtils = Main.BBOUtils

import .MSInstrument.MSErrorModel

function prior_probability_log(model::MSErrorModel, intensity_range::AbstractArray)
  res = 0.0
  zPrior = Cauchy(0, 1)
  res += logpdf(zPrior, model.zDetectionIntercept)
  res += logpdf(zPrior, model.zDetectionFactor)
  stdPrior = Laplace(0, 0.01)
  # minimize std over the intensity range
  res += sum(x -> logpdf(stdPrior, MSInstrument.signal_std(model, x)/x),
             intensity_range)
  res += logpdf(Beta(10.0, 1.1), model.detectionMax)
  return res
end

function reference_intensities(intensities::AbstractArray{<:Real};
                               quantile_range::Union{Tuple{Number, Number}, Number} = 0.99,
                               nbins::Integer=20)
    if quantile_range isa Number
        qtile_min = (1.0 - quantile_range)*0.5
        qtile_max = 1.0 - (1.0 - quantile_range)*0.5
    else
        qtile_min, qtile_max = quantile_range
    end
    logintens_min = floor(log(quantile(intensities, qtile_min)))
    logintens_max = ceil(log(quantile(intensities, qtile_max)))
    return exp.(LinRange(logintens_min, logintens_max, nbins))
end

"""
Data for `MSErrorCalibrationProblem`.
"""
struct MSErrorCalibrationData
    zScale::Float64
    zShift::Float64
    intensity_bins::Vector{Float64}

    intensity::Vector{Float64}
    logintensity_predicted::Vector{Float64}
    intensity_zscore_predicted::Vector{Float64}
    intensity_error::Vector{Float64}

    intensities_df::DataFrame
    shifts_df::DataFrame
end

Base.length(data::MSErrorCalibrationData) = length(data.zscore)

function MSErrorCalibrationData(data, mschannels::AbstractDataFrame;
        object::Symbol = :protgroup,
        kwargs...
)
    if object === :protgroup
        return import_protgroup_data(data, mschannels; kwargs...)
    elseif object === :pepmodstate
        return import_pepmodstate_data(data, mschannels; kwargs...)
    else
        throw(ArgumentError("Unsupported object type: $object"))
    end
end

function import_protgroup_data(data, mschannels::AbstractDataFrame;
    intensity_prefix::AbstractString="Intensity ",
    kwargs...
)
    @info "Preparing protein groups instrument calibration data..."
    protgroupsWide = data[:protgroups]::DataFrame
    protgroups = protgroupsWide[[:id, Symbol("Protein IDs"), Symbol("Majority protein IDs")]]#@where(:is_contaminant.==0 & :is_reverse.==0)
    rename!(protgroups, :id => :protgroup_id)
    protgroups.use_for_calib = .!protgroups.is_contaminant .& .!protgroups.is_reverse
    intensity_cols = filter(col -> startswith(string(col), intensity_prefix),
                            propertynames(protgroupsWide))
    intensities = vcat([DataFrame(protgroup_id=protgroupsWide.id,
                                  intensity=protgroupsWide[col],
                                  mschannel=replace(string(col), intensity_prefix, ""))
                        for col in intensity_cols]...)
    MSErrorCalibrationData(intensities, mschannels;
                           object_col=:protgroup_id, kwargs...)
end

function import_pepmodstate_data(data, mschannels::AbstractDataFrame;
    intensity_col::Symbol=:intensity,
    allowed_ident_types::Union{Nothing, AbstractArray{<:AbstractString}} = nothing,#["ISO-MSMS", "MULTI-MSMS", "MULTI-SECPEP"]
    kwargs...
)
    @info "Preparing pepmodstate instrument calibration data..."

    # skip contaminants as they are not labeled and might differ between technical replicates
    pepmodstates_df = copy(data.pepmodstates)
    if !hasproperty(pepmodstates_df, :is_contaminant)
        pepmodstates_df = innerjoin(pepmodstates_df, data.peptides[!, [:peptide_id, :is_contaminant, :is_reverse]],
                                    on=:peptide_id)
    end
    pepmodstates_df = filter(r -> !r.is_contaminant && !r.is_reverse, pepmodstates_df)
    # allow peaks identified&quantified only with MULTI-MSMS (+ISO-MSMS for label-free)
    # to be used for noise model inference
    # this should exclude variability due to incorrect feature matching
    # FIXME add keyword to enable/disable this behaviour
    allowed_ident_types_set = isnothing(allowed_ident_types) ? nothing : Set(allowed_ident_types)
    pepmodstates_df = innerjoin(pepmodstates_df, combine(groupby(data.observations, :pepmodstate_id)) do pms_obs_df
        DataFrame(use_for_calib = all(r -> ((isnothing(allowed_ident_types_set) ||
                                             get(r.ident_type) ∈ allowed_ident_types_set)) ||
                                            ismissing(r[intensity_col]),
                                      eachrow(pms_obs_df)))
    end, on=:pepmodstate_id)
    MSErrorCalibrationData(copy(data.observations), pepmodstates_df, mschannels;
            object_col=:pepmodstate_id,
            intensity_col=intensity_col,
            kwargs...)
end

# low-level constructor of calibration data
function MSErrorCalibrationData(intensities::DataFrame,
    objects::DataFrame, mschannels::DataFrame;
    object_col::Symbol = :protgroup_id,
    obj_use_col::Symbol = :use_for_calib,
    msrun_col::Symbol = :msrun,
    exp_col::Symbol = :experiment,
    mstag_col::Union{Symbol,Nothing} = nothing,
    intensity_col::Symbol = :intensity,
    nbins::Integer = 20,
    bin_oversize::Number = 2,
    quantile_range::Union{Tuple{Number, Number}, Number} = 0.99,
    max_objXexp::Union{Integer, Nothing} = nothing,
    missed_intensity_factor::Number = 1.0, # how much to adjust the predicted intensity for each missing value
    error_scale::Number = 1.0 # how much the predicted error is scaled (would like to scale down if data are not tech. replicates)
)
    mschannel_idcols = [msrun_col]
    isnothing(mstag_col) || push!(mschannel_idcols, mstag_col)
    mschannel_cols = vcat(mschannel_idcols, [exp_col])
    mschannels = mschannels[!, mschannel_cols]
    if !isnothing(mstag_col)
        expchannel_col = Symbol(string(exp_col), "X", string(mstag_col))
        mschannels[!, expchannel_col] = categorical(coalesce.(mschannels[!, exp_col], "") .*
                                                    "#" .* coalesce.(mschannels[!, mstag_col], ""))
    else
        expchannel_col = exp_col
    end
    expchannels = combine(groupby(mschannels, expchannel_col)) do df
        DataFrame(n_mschannels = length(unique(df[!, msrun_col])))
    end
    expchannels_used = Set(expchannels[expchannels.n_mschannels .> 1, expchannel_col])
    mschannels.is_used = mschannels[!, expchannel_col] .∈ Ref(expchannels_used)
    @info "Using $(sum(mschannels.is_used)) of $(size(mschannels, 1)) mschannels(s) (tech. replicates)"

    intensities = innerjoin(semijoin(intensities[[!ismissing(i) && isfinite(i) && i > 0.0 for i in intensities[!, intensity_col]],
                                        vcat([intensity_col, object_col], mschannel_idcols)],
                            objects[objects[!, obj_use_col], [object_col]], on=object_col),
                       mschannels[mschannels.is_used, vcat(mschannel_idcols, [expchannel_col])],
                       on=mschannel_idcols)
    @assert !any(ismissing, intensities[!, intensity_col])
    # all objectXexperiment pairs, where at least one quantification available
    objXexp = unique(intensities[!, [object_col, expchannel_col]])
    objects = semijoin(objects, objXexp, on=object_col)
    @info "Using $(nrow(intensities)) quantification(s) of $(nrow(objects)) object(s) in $(
        nrow(unique(intensities[!, mschannel_idcols]))) mschannel(s)"

    intens_qtl_col = Symbol(string(intensity_col), "_qtl")
    intens_rank_col = Symbol(string(intensity_col), "_rank")
    (intensity_col != :intensity) && rename!(intensities, intensity_col => :intensity)
    intensities.logintensity = log.(intensities.intensity)

    @info "Calculating MS run shifts..."
    shifts_df = combine(groupby(semijoin(intensities, objects, on=object_col), expchannel_col)) do df
        df = combine(groupby(df, msrun_col)) do tech_df
            res = copy(tech_df)
            res.intensity_rank = invperm(sortperm(res.logintensity))
            res.intensity_qtl = values(res.intensity_rank) ./ size(res, 1)
            select(res, Not(msrun_col))
        end
        tech_replicates = sort(unique(values(df[!, msrun_col])))
        # data of the first tech. replicate
        tech_repl_1st = length(tech_replicates) > 0 ? tech_replicates[1] : nothing
        df_1st = df[df[!, msrun_col] .== tech_repl_1st,
                    [msrun_col, expchannel_col, object_col,
                     :logintensity, :intensity_qtl]]
        msrun1_col = Symbol(string(msrun_col), "_1st")
        rename!(df_1st, msrun_col => msrun1_col,
                        :logintensity => :logintensity_1st,
                        :intensity_qtl => :intensity_qtl_1st)
        # compare first replicate with all others
        paired_df = innerjoin(df, df_1st, on=[expchannel_col; object_col])
        combine(groupby(paired_df, [msrun_col])) do subdf
            w = weights(0.5.*sqrt.(subdf.intensity_qtl .+ subdf.intensity_qtl_1st))
            logintensity_shifts = subdf.logintensity .- subdf.logintensity_1st
            med_shift = median(logintensity_shifts, w)
            DataFrame(logshift = med_shift,
                      logshift_sd = std(logintensity_shifts .- med_shift, w),
                      n = size(subdf, 1))
        end
    end

    shifts_df.mult = exp.(shifts_df.logshift)
    @info "MS experiments normalization:\n$shifts_df"

    intensities_expanded = innerjoin(semijoin(objXexp, objects[objects[!, obj_use_col], [object_col, obj_use_col]],
                                              on=object_col),
                                     shifts_df[!, [expchannel_col, msrun_col, :mult]], on=[expchannel_col])
    join_cols = [expchannel_col, msrun_col, object_col]
    intensities_expanded = leftjoin(intensities_expanded,
                                    intensities[!, [join_cols; [:intensity, :logintensity]]],
                                    on=join_cols)
    intensities_expanded.intensity_norm = intensities_expanded.intensity ./ intensities_expanded.mult
    @info "Missing quantifications in the extended intensities table: $(
            count(ismissing, intensities_expanded.intensity)))"

    # predict average object intensities in each experiment
    expchanXobj_df = combine(groupby(intensities_expanded, [expchannel_col, object_col])) do expobjs_df
        all(ismissing, expobjs_df.intensity_norm) && return DataFrame() # no intensities for given feature (in a given label), skip
        missed_factor = missed_intensity_factor^count(ismissing, expobjs_df.intensity)
        # FIXME == 1 is a workaround for mapreduce() bug
        predicted = size(expobjs_df, 1) == 1 ? expobjs_df.intensity_norm[1] :
                    median(skipmissing(expobjs_df.intensity_norm))
        DataFrame(intensity_norm_predicted = predicted,
                  intensity_norm_predicted_corr = predicted .* missed_factor)
    end
    @info "Found $(nrow(expchanXobj_df)) unique experiment×channel×object usable for calibration"

    ref_intensities = reference_intensities(expchanXobj_df.intensity_norm_predicted_corr,
                                            quantile_range=quantile_range, nbins=nbins)
    @info "Using intensities range $(extrema(ref_intensities))"
    expchanXobj_df.intensity_bin = searchsortedfirst.(Ref(ref_intensities), expchanXobj_df.intensity_norm_predicted_corr)

    intensities_expanded = leftjoin(intensities_expanded, expchanXobj_df, on=[expchannel_col, object_col])
    intensities_expanded.intensity_predicted = intensities_expanded.intensity_norm_predicted_corr .* intensities_expanded.mult
    intensities_expanded.logintensity_predicted = log.(intensities_expanded.intensity_predicted)
    intensities_expanded.intensity_prediction_error = intensities_expanded.intensity_predicted .- intensities_expanded.intensity

    usable_perbin = countmap(expchanXobj_df.intensity_bin)
    min_nbin = quantile(values(usable_perbin), 0.1)
    use_perbin = ceil(Int, bin_oversize*min_nbin)
    if !isnothing(max_objXexp) && use_perbin*length(usable_perbin) > max_objXexp
        use_perbin = fld1(max_objXexp, length(usable_perbin))
    end
    expchanXobj_df[!, :use_for_calib] .= false
    expchanXobj_df = combine(groupby(expchanXobj_df, :intensity_bin)) do binobjs_df
        binobjs_df[sample(1:nrow(binobjs_df), min(nrow(binobjs_df), use_perbin)), :use_for_calib] .= true
        return binobjs_df
    end
    @info "Using $(use_perbin) experiment×channel×object tuple(s) for each of $(length(usable_perbin)) intensity bin(s): $(sum(expchanXobj_df.use_for_calib)) in total"

    intensities_used = semijoin(intensities_expanded,
            select(filter(r -> r.use_for_calib, expchanXobj_df), [expchannel_col, object_col]),
            on=[expchannel_col, object_col])

    # WARNING used data has different (log-uniform) distribution than original data (log-normal),
    # so the zscores are biased
    zShift = mean(skipmissing(intensities_used.logintensity))
    zScale = 1.0/std(skipmissing(intensities_used.logintensity))
    intensities_used.intensity_zscore = (intensities_used.logintensity .- zShift) .* zScale
    intensities_used.intensity_prediction_error = intensities_used.intensity_predicted .- intensities_used.intensity
    intensities_used.intensity_zscore_predicted = (intensities_used.logintensity_predicted .- zShift) .* zScale

    return MSErrorCalibrationData(zScale, zShift, ref_intensities,
                                  coalesce.(intensities_used.intensity, NaN),
                                  coalesce.(intensities_used.logintensity_predicted, NaN),
                                  coalesce.(intensities_used.intensity_zscore_predicted, NaN),
                                  coalesce.(error_scale .* intensities_used.intensity_prediction_error, NaN),
                                  intensities_used,
                                  shifts_df)
end

reference_intensities(data::MSErrorCalibrationData) =
    reference_intensities(data.intensity_bins)

function likelihood_log(model::MSErrorModel, data::MSErrorCalibrationData)
    res = 0.0
    noiseDistr = Laplace(0.0, 1.0)
    min_detect_logprob = log(1E-3) # TODO configurable
    @inbounds for i in eachindex(data.intensity_zscore_predicted)
        if isfinite(data.intensity_error[i])
            logprec = -MSInstrument.zscore_logstd(model, data.intensity_zscore_predicted[i])
            if data.intensity_error[i] != 0 # skip 0 errors (single quant for an object)
                res += logprec + logpdf(noiseDistr, exp(logprec) * data.intensity_error[i])
            end
            res += max(min_detect_logprob, MSInstrument.detection_likelihood_log(model, true, data.logintensity_predicted[i]))
        else
            res += max(min_detect_logprob, MSInstrument.detection_likelihood_log(model, false, data.logintensity_predicted[i]))
        end
        #print( "error=", data.error[i] ," zscore=", data.zscore[i], " expected=", data.expected_log[i], " τ=", zscore_precision(model, data.zscore[i]), " llh=", res, " " )
        #if !isfinite(res) error("Infinite likelihood at $i: $model") end
    end
    return res
end

struct MSErrorModelFactory
    # log(signal) to zscore transformation
    zScale::Float64
    zShift::Float64
    param_bounds::Vector{ParamBounds}

    function MSErrorModelFactory(
        zScale::Float64,
        zShift::Float64;
        detectionMaxMin::Number = 0.0
    )
        zScale > 0 || throw(ArgumentError("zScale should be positive ($zScale found)"))
        new(zScale, zShift,
            [(0.0, 100.0), (-50.0, 50.0), # zDetectionFactor & intercept
             (Float64(detectionMaxMin), 1.0), # detectionMax
             (0.0, 1.0/zScale), (0.0, 1.0/zScale), # scale hi/lo
             (0.0, 50.0), (-50.0, 50.0), # offset & bend
             (0, 50.0)] # smooth
            )
    end
end

nparams(factory::MSErrorModelFactory) = length(factory.param_bounds)

# uninitialized model
emptymodel(factory::MSErrorModelFactory) = MSErrorModel(factory.zScale, factory.zShift)

function params2model!(model::MSErrorModel, factory::MSErrorModelFactory, params::AbstractVector)
    length(params)==nparams(factory) || throw(DimensionMismatch("Length of parameters vector and factory model parameters differ"))
    model.zDetectionFactor = params[1]
    model.zDetectionIntercept = params[2]
    model.detectionMax = params[3]
    model.sigmaScaleHi = params[4]
    model.sigmaScaleLo = params[5]
    model.sigmaOffset = params[6]
    model.sigmaBend = params[7]
    model.sigmaSmooth = params[8]
    MSInstrument.sync!(model)
end

params2model(factory::MSErrorModelFactory, params::AbstractVector) = params2model!(emptymodel(factory), factory, params)

struct BayesianLogprobAggregator
    prior_weight::Float64

    BayesianLogprobAggregator(prior_weight::Real = 1.0) = new(prior_weight)
end

(agg::BayesianLogprobAggregator)(score::NTuple{2, Float64}) =
    agg.prior_weight * score[1] + score[2]

const BayesianLogprobMaximization =
    ParetoFitnessScheme{2, Float64, false, BayesianLogprobAggregator}

BayesianLogprobMaximization(; prior_weight::Real=1.0) =
    ParetoFitnessScheme{2}(fitness_type=Float64, is_minimizing=false,
                           aggregator=BayesianLogprobAggregator(prior_weight))

prior_weight(fitscheme::BayesianLogprobMaximization) = aggregator(fitscheme).prior_weight

"""
`BlackBoxOptim` wrapper for fitting `MSErrorModel`.
"""
mutable struct MSErrorCalibrationProblem{FS} <: OptimizationProblem{FS}
    fitness_scheme::FS
    factory::MSErrorModelFactory
    data::MSErrorCalibrationData
    intensity_range::Vector{Float64}
    search_space::ContinuousRectSearchSpace
    models_pool_lock::Threads.SpinLock
    models_pool::Vector{MSErrorModel}

    MSErrorCalibrationProblem(factory::MSErrorModelFactory, data::MSErrorCalibrationData;
                              fitness_scheme::FitnessScheme = BayesianLogprobMaximization(),
                              intensity_range = reference_intensities(data)) =
        new{typeof(fitness_scheme)}(fitness_scheme, factory, data,
            intensity_range,
            RectSearchSpace(factory.param_bounds),
            Threads.SpinLock(), Vector{MSErrorModel}())
end

BlackBoxOptim.name(::MSErrorCalibrationProblem{FS}) where FS = "MSErrorCalibrationProblem{$FS}"

Base.copy(problem::MSErrorCalibrationProblem) =
    MSErrorCalibrationProblem(problem.factory, problem.data, fitness_scheme=fitness_scheme(problem))

BlackBoxOptim.show_fitness(io::IO, score::NTuple{2,Float64}, problem::MSErrorCalibrationProblem) =
    @printf(io, "(prior=%.3f llh=%.3f: agg=%.3f)",
            score[1], score[2],
            BlackBoxOptim.aggregate(score, fitness_scheme(problem)))

BlackBoxOptim.show_fitness(io::IO, score::IndexedTupleFitness{2,Float64}, problem::MSErrorCalibrationProblem) =
    BlackBoxOptim.show_fitness(io, score.orig, problem)

function prior_and_likelihood_logs(x, p::MSErrorCalibrationProblem)
    lock(p.models_pool_lock)
    tmp_model = isempty(p.models_pool) ? emptymodel(p.factory) : pop!(p.models_pool)
    unlock(p.models_pool_lock)
    params2model!(tmp_model, p.factory, x)
    res = (prior_probability_log(tmp_model, p.intensity_range),
           likelihood_log(tmp_model, p.data))
    lock(p.models_pool_lock)
    push!(p.models_pool, tmp_model)
    unlock(p.models_pool_lock)
    return res
end

# Evaluate fitness of a candidate solution on the 1st objective function of a problem.
BlackBoxOptim.fitness(x, p::MSErrorCalibrationProblem{BayesianLogprobMaximization}) =
    prior_and_likelihood_logs(x, p)

BlackBoxOptim.fitness(x, p::MSErrorCalibrationProblem{<:ScalarFitnessScheme}) =
    sum(prior_and_likelihood_logs(x, p))

function BlackBoxOptim.bbsetup(data::MSErrorCalibrationData;
                      fitness_scheme=BayesianLogprobMaximization(),
                      detectionMaxMin::Number=0.0,
                      kwargs...)
    problem = MSErrorCalibrationProblem(
        MSErrorModelFactory(data.zScale, data.zShift, detectionMaxMin=detectionMaxMin), data,
        fitness_scheme=fitness_scheme)
    return bbsetup(BBOUtils.DefaultFactory(problem; kwargs...); kwargs...)
end

calibfileprefix(info) =
    "mscalib_$(info.msinstrument)_$(info.quanttype)_$(info.quantobj)_$(info.id)_$(info.calib_ver)"

function save_result(filename::AbstractString, problem::MSErrorCalibrationProblem,
                     result::BlackBoxOptim.OptimizationResults;
                     data::Union{MSErrorCalibrationData, Nothing}=nothing,
                     info=nothing, verbose::Bool=false)
    fileext = splitext(filename)[2]
    if fileext == ".json"
        format = :json
    elseif fileext == ".jld2"
        format = :jld2
    else
        throw(ArgumentError("Unknown/unsupported file format: $fileext"))
    end

    model = params2model(problem.factory, best_candidate(result))
    verbose && @info "Best instrument params: $model"

    verbose && @info "Saving optimization results into $filename..."
    if format == :jld2
        @save(filename, info, model, result, data)
    elseif format == :json
        open(filename, "w") do io JSON.print(io, Dict{Symbol,Any}(
            :info => info,
            :mscalib => model,
            :fitness => best_fitness(result),
            :raw_params => best_candidate(result),
        )) end
    end
    nothing
end

end
