module OptCoverUtils

using DataFrames, OptEnrichedSetCover

const FrameUtils = Main.FrameUtils

function collection2mosaic(sets::AbstractDict{S, Set{T}},
                           ref_sets::Union{AbstractDict{S, Set{TR}}, Nothing} = nothing,
                           ref_objs::Union{Set{TR}, Nothing} = nothing;
                           setXset_frac_extra_elms::Real=0) where {S, T, TR}
    all_objs = foldl(union!, values(sets), init=Set{T}())
    if ref_sets !== nothing
        # all annotated + observed ACs (for set relevance calculation)
        if ref_objs === nothing
            ref_objs = Set{TR}()
        end
        all_refs = ref_sets !== nothing ?
                   foldl(union!, values(ref_sets), init=copy(ref_objs)) :
                   nothing
        set_relevances = Dict(begin
            ref_set = ref_sets[set_id]
            set_id => set_relevance(length(intersect(ref_set, ref_objs)),
                                    length(ref_set), length(ref_objs), length(all_refs))
        end for set_id in keys(sets))
    else
        set_relevances = nothing
    end
    SetMosaic(sets, all_objs, set_relevances,
              setXset_nextra_elms=round(Int, setXset_frac_extra_elms * length(all_objs)))
end

function collections2mosaics(colls::AbstractDict{Symbol, <:AbstractDict},
                             ref_colls::Union{AbstractDict{Symbol, <:AbstractDict}, Nothing} = nothing,
                             ref_objs::Union{Set, Nothing} = nothing;
                             setXset_frac_extra_elms::Real=0,
                             verbose::Bool = false)
    Dict(begin
        verbose && @info("preparing $coll_id mosaic...")
        mosaic = collection2mosaic(sets, ref_colls !== nothing ? ref_colls[coll_id] : nothing, ref_objs,
                                   setXset_frac_extra_elms=setXset_frac_extra_elms)
        verbose && @info("done $coll_id ($(length(sets)) set(s), $(ntiles(mosaic)) tile(s))")
        coll_id => mosaic
    end for (coll_id, sets) in pairs(colls))
end

# automatically detect max_log_pvalue threshold
function automask(mosaic::SetMosaic, masks;
                  max_sets::Integer=2500, min_nmasked::Integer=2,
                  max_setsize::Union{Integer, Nothing} = nothing,
                  max_min_pvalue::Number=0.01, max_pvalue::Number = min(max_min_pvalue*10, 1.0),
                  verbose::Bool = false)
    max_log10_pvalue = log10(max_pvalue)
    max_min_log10_pvalue = log10(max_min_pvalue)
    while true
        masked_mosaic = mask(mosaic, masks,
                             min_nmasked=min_nmasked, max_setsize=max_setsize,
                             max_overlap_logpvalue=max_log10_pvalue*log(10),
                             max_min_overlap_logpvalue=max_min_log10_pvalue*log(10))
        verbose && @info("$(nsets(masked_mosaic)) set(s) (max set⋂mask log₁₀(P-value)=$max_log10_pvalue min(log₁₀(P-value))=$max_min_log10_pvalue)")
        if nsets(masked_mosaic) <= max_sets
            return masked_mosaic
        end
        max_log10_pvalue -= 0.5
        max_min_log10_pvalue -= 0.5
    end
end

# helper function for unnest() with the explicit key type of the unnested dict
function unnest(key::Function, ::Type{K}, ddict::Dict) where K
    V = valtype(valtype(ddict))
    nvals = sum(length, values(ddict))
    keys = sizehint!(Vector{K}(), nvals)
    vals = sizehint!(Vector{V}(), nvals)
    for (dict_key, dict) in ddict
        for (el_key, el) in dict
            k = key(dict_key, el_key)
            push!(keys, k)
            push!(vals, el)
        end
    end
    return Dict(zip(keys, vals))
end

"""
    unnest([key::Function, ]ddict::Dict{KO, Dict{KI, V}}) where {KO, KI, V}

Convert dict of dicts into a dict with a keys returned by `key(outerkey::KO, innerkey::KI)`.
`key` defaults to a tuple `(outerkey, innerkey)`.
"""
function unnest(key::Function, ddict::Dict{KO, Dict{KI, V}}) where {KO, KI, V}
    dd1k = first(keys(ddict))
    d1k = first(keys(first(values(ddict))))
    K = typeof(key(dd1k, d1k))
    unnest(key, K, ddict)
end

unnest(ddict::Dict) = unnest((dk, ek)->(dk, ek), ddict)

tuplecol(df::AbstractDataFrame, cols::NTuple{N, Symbol},
         rowixs = 1:size(df, 1)) where N =
    NTuple{N, Any}[ntuple(i -> df[j, cols[i]], N) for j in rowixs]

"""
    filter_multicover(multicover_df::AbstractDataFrame;
                      term_cols=[:term_id], set_cols=[:cluster],
                      max_term_pvalue=1E-4, max_set_pvalue=max_term_pvalue,
                      max_entry_pvalue=max_set_pvalue)

Filters the multicover report in 3 steps:
 * keep only the terms that have at least one entry with
   `set_overlap_log10pvalue` ≤ `log10(max_term_pvalue)`
 * keep only the sets that have at least one entry with
   `set_overlap_log10pvalue` ≤ `log10(max_set_pvalue)` for one of the selected terms
 * keep only report entries for the selected terms and sets that have
   `set_overlap_log10pvalue` ≤ `log10(max_entry_pvalue)`
"""
function filter_multicover(multicover_df::AbstractDataFrame;
                           term_cols=[:term_id], set_cols=[:cluster],
                           max_term_pvalue::Number=1E-4, min_term_overlap::Union{Integer, Nothing}=2,
                           max_set_pvalue::Union{Number, Nothing}=max_term_pvalue, min_set_overlap::Union{Integer, Nothing}=min_term_overlap,
                           max_entry_pvalue::Union{Number, Nothing}=nothing, min_entry_overlap::Union{Integer, Nothing}=nothing)
    minolap_term = isnothing(min_term_overlap) ? 0 : min_term_overlap
    minolap_set = isnothing(min_set_overlap) ? 0 : min_set_overlap
    minolap_entry = isnothing(min_entry_overlap) ? 0 : min_entry_overlap
    maxpval_term = isnothing(max_term_pvalue) ? 1.0 : max_term_pvalue
    maxpval_set = isnothing(max_set_pvalue) ? 1.0 : max_set_pvalue
    maxpval_entry = isnothing(max_entry_pvalue) ? 1.0 : max_entry_pvalue

    (maxpval_set < maxpval_term) &&
        @warn("max_set_pvalue ($max_set_pvalue) is more stringent than max_term_pvalue ($max_term_pvalue)")
    (maxpval_entry < maxpval_set) &&
        @warn("max_entry_pvalue ($max_entry_pvalue) is more stringent than max_set_pvalue ($max_set_pvalue)")
    (minolap_set > minolap_term) &&
        @warn("min_set_overlap ($min_set_overlap) is more stringent than min_term_overlap ($min_term_overlap)")
    (minolap_entry > minolap_set) &&
        @warn("min_entry_overlap ($min_entry_overlap) is more stringent than min_set_overlap ($min_set_overlap)")

    terms_v = tuplecol(multicover_df, tuple(term_cols...))
    sets_v = tuplecol(multicover_df, tuple(set_cols...))
    term_ids = Set(terms_v[(multicover_df.nmasked .>= max(minolap_term, minolap_set, minolap_entry)) .&
                           (multicover_df.set_overlap_log10pvalue .<= log10(min(maxpval_term, maxpval_set, maxpval_entry)))])
    set_ids = Set(sets_v[(terms_v .∈ Ref(term_ids)) .&
                         (multicover_df.nmasked .>= max(minolap_set, minolap_entry)) .&
                         (multicover_df.set_overlap_log10pvalue .<= log10(min(maxpval_set, maxpval_entry)))])
    @info "  $(length(term_ids))/$(length(Set(terms_v))) significant terms, $(length(set_ids))/$(length(Set(sets_v))) significant sets"
    return multicover_df[(sets_v .∈ Ref(set_ids)) .& (terms_v .∈ Ref(term_ids)) .&
                         (multicover_df.nmasked .>= minolap_entry) .&
                         (multicover_df.set_overlap_log10pvalue .<= log10(maxpval_entry)), :]
end

include(joinpath(@__DIR__, "covers_report.jl"))

end
