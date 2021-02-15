module ProtgroupXMatch

using DataFrames
using Printf: @sprintf

"""
References a protein group by its `row_ix` in a specific MQ dataset.

`pg_rank` defines the protgroup-related quality of the reference, the lower the better.
`ac_rank` defines the AC-related quality of the reference, the lower the better.
"""
struct ProtgroupRef
    row_ix::Int
    pg_rank::Int
    ac_rank::Int

    ProtgroupRef(row_ix::Int=0, pg_rank::Int=-1, ac_rank::Int=-1) =
        new(row_ix, pg_rank, ac_rank)
end

Base.isempty(ref::ProtgroupRef) = ref.row_ix <= 0

pg_rank(ref::ProtgroupRef) = ref.pg_rank
ac_rank(ref::ProtgroupRef) = ref.ac_rank

"""
Matches protein groups between `n` MQ datasets.
"""
struct ProtgroupCrossref
    refs::Vector{ProtgroupRef}

    ProtgroupCrossref(n::Int) = new(fill(ProtgroupRef(), n))
end

Base.length(xr::ProtgroupCrossref) = length(xr.refs)

Base.@propagate_inbounds Base.getindex(xr::ProtgroupCrossref, i::Int) = xr.refs[i]

Base.@propagate_inbounds Base.setindex!(xr::ProtgroupCrossref, ref::ProtgroupRef, i::Int) =
    xr.refs[i] = ref

# highest rankf(ref) of non-empty references
function worst_rank(rankf, xr::ProtgroupCrossref)
    res = -1
    for ref in xr.refs
        isempty(ref) && continue
        r = rankf(ref)
        (r > res) && (res = r)
    end
    return res
end

# lowest rankf(ref) of non-empty references
function best_rank(rankf, xr::ProtgroupCrossref)
    res = -1
    for ref in xr.refs
        isempty(ref) && continue
        r = rankf(ref)
        if (res == -1) || (r < res)
            res = r
        end
    end
    return res
end

const UniprotIsoformRegex = r"^(\w+)-(?:PRO_)?(\d+)$"

function strip_uniprot_isoform(ac::AbstractString)
    acmatch = match(UniprotIsoformRegex, ac)
    acmatch === nothing ? ac : convert(String, acmatch[1])
end

# ranks uniprot DB entry
rank_uniprot(uprot::DataFrameRow) =
    ifelse(ismissing(uprot.gene_name), 10, 0) +
    begin
        acmatch = match(UniprotIsoformRegex, uprot.protein_ac)
        isnothing(acmatch) ? 0 : parse(Int, acmatch[2]) - 1
    end +
    ifelse(uprot.src_db == "sp", 0, 3) +
    convert(Int, coalesce(uprot.protein_existence, 7)) - 1

"""
Checks if non-empty protein group references in `xr1` and `xr2` are the same.
"""
function iscompatible(xr1::ProtgroupCrossref, xr2::ProtgroupCrossref)
    @assert length(xr1.refs) == length(xr2.refs)
    has_overlap = false
    for i in eachindex(xr1.refs)
        @inbounds ref1 = xr1.refs[i]
        @inbounds ref2 = xr2.refs[i]
        if !isempty(ref1) && !isempty(ref2)
            has_overlap = ref1.row_ix == ref2.row_ix
            has_overlap || return false
        end
    end
    return has_overlap
end

"""
Copies non-empty `xr2` protein groups references into `xr1`.
`xr1` and `xr2` should be compatible.
"""
function update!(xr1::ProtgroupCrossref, xr2::ProtgroupCrossref)
    @assert length(xr1.refs) == length(xr2.refs)
    @inbounds for i in eachindex(xr1.refs)
        have_1 = !isempty(xr1.refs[i])
        have_2 = !isempty(xr2.refs[i])
        @assert !have_1 || !have_2 || (xr1.refs[i].row_ix == xr2.refs[i].row_ix)
        if !have_1 && have_2
            xr1.refs[i] = xr2.refs[i]
        elseif have_1 && have_2 &&
            (xr1.refs[i].row_ix == xr2.refs[i].row_ix) &&
            (xr1.refs[i].pg_rank >= xr2.refs[i].pg_rank) &&
            (xr1.refs[i].ac_rank >= xr2.refs[i].ac_rank)
            # improve match rank
            xr1.refs[i] = xr2.refs[i]
        end
    end
    return xr1
end

"""
Applies `f` over `cols` data of the best-pg-ranked proteing groups references of `xr` and
reduces the result using `op` binary function.

`op` should return a tuple of the reduction result and a boolean flag
indicating whether to continue the reduction.
"""
function mapreduce(f, op, xr::ProtgroupCrossref, cols::AbstractVector, init,
                   max_pg_rank=best_rank(pg_rank, xr))
    res = init
    @assert length(cols) == length(xr)
    @inbounds for (i, ref) in enumerate(xr.refs)
        isempty(ref) && continue
        (pg_rank(ref) > max_pg_rank) && continue
        col_val = cols[i][ref.row_ix]
        f_val = f(col_val)
        res, inprogress = op(res, f_val)
        inprogress || break
    end
    res
end

"""
Collect all unique values that are referred by all best-ranked matches of `xr`.
The values are stored as `sep`-separated strings in `cols`.
"""
function intersect_dlm_values(xr::ProtgroupCrossref, cols::AbstractVector,
                              val_scores::Union{Dict, Nothing} = nothing,
                              def_score=0, max_pg_rank=best_rank(pg_rank, xr),
                              sep=";")
    T = nonmissingtype(eltype(first(cols)))
    res = mapreduce(identity, function (acc, vals_str)
        ismissing(vals_str) && return (acc, true)
        vals = [T(val) for val in split(vals_str, sep)]
        if isempty(acc) # first acc assignment
            return (vals, true)
        else
            _vals = Set{T}(vals)
            filter!(val -> val ∈ _vals, acc)
            return (acc, !isempty(acc)) # don't continue if the intersection is zero
        end
    end, xr, cols, Vector{T}(), max_pg_rank)
    if !isempty(res) && (val_scores !== nothing)
        sort!(res, by = x -> (get(val_scores, x, def_score), x))
    end
    isempty(res) ? missing : join(res, sep)
end

"""
Check if any best-ranked matches of `xr` have `true` in `cols`.
"""
Base.any(xr::ProtgroupCrossref, cols::Vector) =
    mapreduce(identity, function (acc, val)
        ismissing(val) ? (acc, true) : (val, !val)
    end, xr, cols, missing)

# build a dict of ACxPGR -> highest protgroup_ref
function protgrouprefs(protgroups_df::AbstractDataFrame,
                       protein_acs_penalties::Dict{String})
    refs = Dict{Tuple{String, Int}, ProtgroupRef}()
    miss_penalty = isempty(protein_acs_penalties) ? 1 : maximum(values(protein_acs_penalties))+1
    for (ac_col, acs_col) in enumerate(intersect([:majority_protein_acs, :protein_acs], propertynames(protgroups_df)))
        acs_data = protgroups_df[!, acs_col]
        pgr = ac_col-1
        for (row_ix, acs_str) in enumerate(acs_data)
            if !ismissing(acs_str)
                acs = split(acs_str, ";")
                for ac in acs
                    new_ref = ProtgroupRef(row_ix, pgr, get(protein_acs_penalties, ac, miss_penalty))
                    # put the ACxPGR if not yet there (with the higher rank)
                    old_ref = get!(refs, (ac, pgr), new_ref)
                    if old_ref.ac_rank > new_ref.ac_rank
                        refs[(ac, pgr)] = new_ref
                    end
                end
            end
        end
    end
    return refs
end

function match_protgroups(protgroup_dfs::AbstractVector{Pair{K, T}},
                          protein_acs_penalties::Dict{String}) where {K, T}
    # build ac-to-pgs map
    @info("Building $(length(protgroup_dfs))-way AC↔protgroup_id map...")
    n_ds = length(protgroup_dfs)
    ds_names = first.(protgroup_dfs)
    pg_dfs = last.(protgroup_dfs)
    colnames = reduce(union!, propertynames.(pg_dfs), init=Vector{Symbol}())
    common_colnames = reduce(intersect!, propertynames.(pg_dfs), init=copy(colnames))
    ac2xrefs = Dict{Tuple{String, Int, Int}, ProtgroupCrossref}()
    for (ds_ix, pg_df) in enumerate(pg_dfs)
        ds_refs = protgrouprefs(pg_df, protein_acs_penalties)
        # append pg_ref to all-datasets dict with the same match ranks
        for ((ac, pgr), ref) in ds_refs
            xr = get!(ac2xrefs, (ac, ref.pg_rank, ref.ac_rank)) do
                ProtgroupCrossref(n_ds)
            end
            xr[ds_ix] = ref
        end
    end
    function check_column(colname)
        if colname in common_colnames
            return true
        else
            @warn "$colname is not present in all datasets, skipping"
            return false
        end
    end

    @info("$(length(ac2xrefs)) AC→PG reference(s)")
    @info("Stripping ACs from protgroup_id↔protgroup_id↔… map and removing duplicates...")
    # eliminate duplicated and sub- matches
    xrefs = Vector{ProtgroupCrossref}()
    for new_xr in values(ac2xrefs)
        absorbed = false
        for xr in xrefs
            if iscompatible(new_xr, xr)
                update!(xr, new_xr)
                absorbed = true
            end
        end
        # the new match is not compatible to any others
        absorbed || push!(xrefs, new_xr)
    end

    @info("Converting the protgroup matches into DataFrame...")
    res = DataFrame(best_pg_rank = [best_rank(pg_rank, xr) for xr in xrefs],
                    worst_pg_rank = [worst_rank(pg_rank, xr) for xr in xrefs],
                    best_ac_rank = [best_rank(ac_rank, xr) for xr in xrefs],
                    worst_ac_rank = [worst_rank(ac_rank, xr) for xr in xrefs])
    @inbounds for (ds_ix, ds_name) in enumerate(ds_names)
        col_ix = [isempty(xr[ds_ix]) ? missing : xr[ds_ix].row_ix for xr in xrefs]
        col_pg_rank = [isempty(xr[ds_ix]) ? missing : pg_rank(xr[ds_ix]) for xr in xrefs]
        col_ac_rank = [isempty(xr[ds_ix]) ? missing : ac_rank(xr[ds_ix]) for xr in xrefs]
        res[!, Symbol("rowix_", ds_name)] = col_ix
        res[!, Symbol("pgrank_", ds_name)] = col_pg_rank
        res[!, Symbol("acrank_", ds_name)] = col_ac_rank
    end
    # remove duplicates (FIXME how to avoid them)
    uniq_mask = .!nonunique(res)
    res = res[uniq_mask, :]
    xrefs = xrefs[uniq_mask]
    @info("$(length(xrefs)) distinct protgroup matches")

    @info("Calculating the protein group data intersection:")
    @info("    * ACs")
    max_score = isempty(protein_acs_penalties) ? 0 : maximum(values(protein_acs_penalties))
    for (coli, coln) in enumerate(intersect(colnames, [:majority_protein_acs, :protein_acs]))
        cols = getindex.(pg_dfs, !, coln)
        res[!, coln] = intersect_dlm_values.(xrefs, Ref(cols), Ref(protein_acs_penalties), max_score+1, coli-1)
    end
    @info("    * gene names")
    for coln in intersect(colnames, [:gene_names, :protein_names])
        if check_column(coln)
            cols = getindex.(pg_dfs, !, coln)
            res[!, coln] = intersect_dlm_values.(xrefs, Ref(cols))
        end
    end
    @info("    * contaminants and reverse sequences")
    for coln in intersect(colnames, [:is_contaminant, :is_reverse])
        if check_column(coln)
            cols = getindex.(pg_dfs, !, coln)
            res[!, coln] = any.(xrefs, Ref(cols))
        end
    end
    @info("    * organism")
    if :organism ∈ colnames && check_column(:organism)
        cols = getindex.(pg_dfs, !, :organism)
        res.organism = intersect_dlm_values.(xrefs, Ref(cols))
    end
    return res
end

end
