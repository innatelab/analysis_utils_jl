module FrameUtils

using DataFrames, CategoricalArrays
using Printf: @sprintf

function categorical!(df::AbstractDataFrame, col::Union{Symbol, String}; levels=nothing, ordered::Bool=false)
    df[!, col] .= categorical(convert(levels === nothing ? Vector : Vector{eltype(levels)}, df[!, col]),levels=levels, ordered=ordered)
    return df
end

categorical!(df::AbstractDataFrame, cols::AbstractVector) =
    reduce((df, col) -> categorical!(df, col), cols, init=df)

# note: modifies dest
function matchcategorical(dest::AbstractCategoricalArray, src::AbstractCategoricalArray;
                          allowmissing::Bool=false)
    @show levels(dest) levels(src)
    if levels(dest) != levels(src)
        levels!(dest, levels(src), allowmissing=allowmissing)
    end
    if isordered(dest) != isordered(src)
        ordered!(dest, isordered(src))
    end
    return dest
end

function matchcategorical(dest::AbstractArray, src::AbstractCategoricalArray;
                          allowmissing::Bool=false)
    return categorical(dest, levels=levels(src), ordered=isordered(src))
end

function matchcategoricals!(dest::AbstractDataFrame, src::AbstractDataFrame;
                            cols::Union{Nothing, AbstractVector{Symbol}}=nothing,
                            allowmissing::Bool=false,
                            verbose::Bool=true)
    @show cols
    if cols === nothing
        _cols = filter(propertynames(src)) do col
            return isa(src[!, col], AbstractCategoricalArray) && (col in propertynames(dest))
        end
    else
        _cols = convert(Vector{Symbol}, cols)
    end
    for col in _cols
        verbose && @info("Matching categorical :$col")
        @show col typeof(dest[!, col]) typeof(src[!, col])
        dest[!, col] = matchcategorical(dest[!, col], src[!, col], allowmissing=allowmissing)
    end
    @info "matchingcategorical done"
    return dest
end

"""
    dropcategoricals!(df::DataFrame)

Convert all categorical columns of `df` into non-categorical equivalents.
"""
function dropcategoricals!(df::DataFrame)
    for col in names(df)
        if df[!, col] isa AbstractCategoricalVector
            df[!, col] = unwrap.(df[!, col])
        end
    end
    return df
end

dropcategoricals(df::AbstractDataFrame; copycols=true) = dropcategoricals!(copy(df, copycols=copycols))

AT(::Type{A}) where {A<:AbstractArray} = A
AT(::Type{A}) where {A<:AbstractArray{T,N}} where {T>:Missing,N} = Array{nonmissingtype(T), N}
AT(::Type{A}) where {A<:AbstractArray{Any,N}} where {N} = A
AT(::Type{A}) where {A<:AbstractCategoricalArray{T,N}} where {T>:Missing,N} = CategoricalArray{nonmissingtype(T), N}
AT(::Type{A}) where {A<:AbstractCategoricalArray{T,N}} where {T, N} = A

# FIXME replace in favor of dropmissing()?
convert2nonmissing(a::AbstractArray) = convert(AT(typeof(a)), a)

function indexunique(::Type{T}, vals::AbstractVector;
                     start::T=one(T), missingindex::Union{T,Missing}=zero(T),
                     sort::Bool=true) where T
    uniquevals = unique(vals)
    sort && sort!(uniquevals)
    offset = start - 1
    val2ix = Dict(val => T(i + offset) for (i, val) in enumerate(uniquevals))
    return get.(Ref(val2ix), vals, missingindex)
end

indexunique(vals::AbstractVector; kwargs...) = indexunique(Int, vals; kwargs...)

# multivalue version of unstack()
function DataFrames.unstack(df::AbstractDataFrame,
                            rowkeys::Union{AbstractVector{<:Union{Integer, Symbol}}, Symbol, Integer},
                            colkey::Union{Integer, Symbol},
                            values::AbstractVector{<:Union{Integer, Symbol}};
                            sep::Union{String,Char}='_',
                            namewidecols::Function = (valcol, colkey, sep) -> string(valcol, sep, colkey),
                            widecol_groups::Symbol=:colkey)
    wide_df = DataFrame()
    widecol_groups ∈ [:value, :colkey] || throw(ArgumentError("Unsupported value widecol_groups=$(widecol_groups)"))
    colgroup_poses = Vector{Int}()
    for val in values
        val_df = unstack(df, rowkeys, colkey, val,
                         renamecols = col->namewidecols(val, col, sep))
        if !isempty(wide_df) # exclude rowkeys cols from all but 1st wide_df
            isequal(wide_df[!, rowkeys], val_df[!, rowkeys]) || error("Non-matching rowkeys")
            select!(val_df, Not(rowkeys))
            if widecol_groups == :value # just append at the back of wide_df
                for col in names(val_df) # should be faster than hcat()?
                    insertcols!(wide_df, col => val_df[!, col])
                end
            elseif widecol_groups == :colkey
                for (i, col) in enumerate(names(val_df)) # should be faster than hcat()?
                    insertcols!(wide_df, colgroup_poses[i], col => val_df[!, col], after=true)
                    colgroup_poses[i] += i # adjust by this and all the preceding columns inserted
                end
            end
        else
            wide_df = val_df
            colgroup_poses = collect(length(rowkeys)+1:size(wide_df, 2))
        end
    end
    return wide_df
end

# multivariable version of DataFrames.stack()
function pivot_longer(df::AbstractDataFrame,
                      id_cols::Union{AbstractVector{<:Union{Integer, Symbol}}, Symbol, Integer, Nothing} = nothing;
                      measure_vars_regex::Regex=r"^(?<value>[^.]+)\.(?<var>[^.]+)$",
                      var_col::Symbol=:variable, value_col::Union{Symbol,Nothing}=nothing)
    has_value_capture = Base.PCRE.substring_number_from_name(measure_vars_regex.regex, :value) >= 0
    mes_matches = match.(Ref(measure_vars_regex), names(df))
    if all(isnothing, mes_matches)
        @warn "No column matches measure_vars_regex"
        return select(df, id_cols)
    end
    # collect different values and corresponding columns in long format
    value2cols = Dict{Symbol, Vector{Pair{Symbol, Symbol}}}()
    for mes_match in mes_matches
        mes_match === nothing && continue
        val_col = has_value_capture ? Symbol(mes_match[:value]) : value_col
        val_col_map = get!(() -> Vector{Pair{Symbol, Symbol}}(), value2cols, val_col)
        push!(val_col_map, Symbol(mes_match.match) => Symbol(mes_match[:var]))
    end
    # stack individual measures
    idcols = id_cols !== nothing ? id_cols : propertynames(df)[mes_matches .== nothing]
    val_long_dfs = [stack(rename!(df[!, i == 1 ? vcat(idcols, first.(wide_cols)) : first.(wide_cols)], wide_cols),
                          last.(wide_cols), variable_name=var_col, value_name=val_col)
                    for (i, (val_col, wide_cols)) in enumerate(value2cols)]
    # remove variable_name column from all but first stacked frames
    vars = val_long_dfs[1][!, var_col]
    for i in 2:length(val_long_dfs)
        @assert isequal(vars, val_long_dfs[i][!, var_col])
        select!(val_long_dfs[i], Not(var_col))
    end
    # concatenate horizontally
    return hcat(val_long_dfs...)
end

# convert N-dim array into dataframe
function array2frame(arr::AbstractArray{<:Number, N},
                     axes::Vararg{AbstractDataFrame, N};
                     data_col=:signal) where N
    df = DataFrame(data_col => vec(arr))
    ninner = 1
    for (i, axis_df) in enumerate(axes)
        axis_rowixs = repeat(1:size(axis_df, 1), inner=ninner, outer=size(df, 1)÷(ninner*size(axis_df, 1)))
        for axis_col in names(axis_df)
            df[!, axis_col] = axis_df[axis_rowixs, axis_col]
        end
        ninner *= size(axis_df, 1)
    end
    return df
end

# helper type-stable kernel function for frame2array
# that fills the resulting N-dimensional array with data_v values
function _frame2array(data_v::AbstractVector{T1},
                      axes_dfs::NTuple{N, DataFrame},
                      axes_ixs::Matrix{Int},
                      default::T2) where {N, T1, T2}
    T = Union{T1, T2}
    @assert size(axes_ixs, 2) == length(data_v)
    res = fill!(Array{T, N}(undef, ntuple(i -> max(size(axes_dfs[i], 1), 1), N)), default)
    @inbounds for i in eachindex(data_v)
        res[CartesianIndex{N}(ntuple(j -> axes_ixs[j, i], N))] = data_v[i]
    end
    return res
end

"""
Convert data frame into N-D array with `data_col` values and
return it along with the axes dataframes (as a 2-tuple).
"""
function frame2array(data_df::AbstractDataFrame,
                     axes::AbstractVector{Vector{Symbol}};
                     data_col::Union{Symbol, Nothing}=:signal,
                     missed=nothing, default=missing,
                     verbose::Bool=false)
    verbose && @info("Extracting the axes and indexing the frame...")
    axes_dfs = Vector{DataFrame}()
    axes_ixs = fill(0, (length(axes), size(data_df, 1)))
    for (i, ax_cols) in enumerate(axes)
        if isempty(ax_cols) # empty degenerated axis
            push!(axes_dfs, DataFrame())
            axes_ixs[i, :] .= 1
            continue
        end
        # group data_df by axis
        data_df_axgrp = groupby(data_df, ax_cols, sort=true)
        # unique axis values
        axis_df = data_df[data_df_axgrp.idx[data_df_axgrp.starts], ax_cols]
        # index data_df along the axis
        axis_i = view(axes_ixs, i, :)
        @inbounds for j in eachindex(data_df_axgrp.starts)
            axis_i[data_df_axgrp.idx[data_df_axgrp.starts[j]:data_df_axgrp.ends[j]]] .= j
        end
        push!(axes_dfs, axis_df)
    end
    data_v = isnothing(data_col) ?
        fill(true, nrow(data_df)) :
        isnothing(missed) ? data_df[!, data_col] : coalesce.(data_df[!, data_col], missed)
    verbose && @info("Reshaping into $(length(axes))-tensor...")
    res = _frame2array(data_v, ntuple(i -> axes_dfs[i], length(axes)), axes_ixs, default)
    return res, axes_dfs
end

function frame2collection!(coll::Dict{S,Set{T}}, df::AbstractDataFrame;
    obj_col::Symbol, set_col::Symbol,
    min_size::Integer=1
) where {T,S}
    gdf = groupby(df, set_col, sort=false)
    sizehint!(coll, length(gdf))
    for set_df in gdf
        size(set_df , 1) >= min_size || continue
        objs = Set{T}(skipmissing(set_df[!, obj_col]))
        if length(objs) >= min_size
            coll[set_df[1, set_col]] = objs
        end
    end
    return coll
end

function frame2collection(df::AbstractDataFrame;
    obj_col::Symbol, set_col::Symbol,
    min_size::Integer=1
)
    S = nonmissingtype(eltype(df[!, set_col]))
    T = nonmissingtype(eltype(df[!, obj_col]))
    frame2collection!(Dict{S, Set{T}}(), df,
        obj_col=obj_col, set_col=set_col, min_size=min_size)
end

function collection2frame(coll::AbstractDict{<:Any, <:Set},
                          set2coll_df::Union{AbstractDataFrame, Nothing}=nothing;
                          setid_col = :set_id, objid_col = :protein_ac,
                          collid_col = :coll_id
)
    len = sum(length, values(coll))
    set_ids = sizehint!(Vector{keytype(coll)}(), len)
    obj_ids = sizehint!(Vector{eltype(valtype(coll))}(), len)
    for (set_id, set_objs) in coll
        for obj_id in set_objs
            push!(set_ids, set_id)
            push!(obj_ids, obj_id)
        end
    end
    res_df = DataFrame(setid_col => set_ids, objid_col => obj_ids)
    if set2coll_df !== nothing
        res_df = leftjoin(res_df, set2coll_df[!, [collid_col, setid_col]],
                          on = :term_id)
    end
    return res_df
end

function frame2collections(df::AbstractDataFrame;
    obj_col::Symbol = :protgroup_id,
    set_col::Symbol = :set_id,
    coll_col::Symbol = :src,
    verbose::Bool = true,
    kwargs...
)
    # convert back to collection of sets
    ObjType = eltype(df[!, obj_col])
    SetType = eltype(df[!, set_col])
    CollType = eltype(df[!, coll_col])
    colls = Dict{Symbol, Dict{SetType, Set{ObjType}}}()
    for coll_subdf in groupby(df, coll_col, sort=false);
        coll_id = Symbol(string(coll_subdf[1, coll_col]))
        verbose && @info("Processing $coll_id annotation collection...")
        coll_sets = frame2collection(coll_subdf; obj_col, set_col, kwargs...)
        colls[coll_id] = coll_sets
        verbose && @info("  $(length(coll_sets)) set(s) processed")
    end
    return colls
end

end
