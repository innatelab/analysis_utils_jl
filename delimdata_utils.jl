module DelimDataUtils

using DataFrames
using Printf: @sprintf

function expand_delim_column(list_v::AbstractVector{T},
                           key_v::Union{AbstractVector, Nothing} = nothing;
                           delim::String = ";") where T
    if key_v !== nothing
        length(key_v) == length(list_v) || throw(DimensionMismatch("Lengths of key_v and list_v differ"))
        exp_key_v = similar(key_v, 0)
        if list_v isa Vector # workaround for missing sizehint!(::WeakRefStrings.StringArray)
            sizehint!(exp_key_v, length(key_v))
        end
    else
        exp_key_v = nothing
    end

    ET = T >: Missing ? Union{String, Missing} : String
    elem_v = sizehint!(Vector{ET}(), length(list_v))
    for (rix, list_el) in enumerate(list_v)
        if !ismissing(list_el)
            list_str = string(list_el)
            elems = split(list_str, delim)
            append!(elem_v, elems)
            if key_v !== nothing
                @inbounds key_el = key_v[rix]
                for j in eachindex(elems)
                    push!(exp_key_v, key_el)
                end
            end
        else
           push!(elem_v, missing)
           (key_v !== nothing) && push!(exp_key_v, key_v[rix])
        end
    end
    #@show typeof(exp_key_v) typeof(elem_v)
    return exp_key_v, elem_v
end

function expand_delim_column(df::DataFrame; list_col::Symbol = :majority_protein_acs,
                           elem_col::Union{Nothing, Symbol} = nothing,
                           key_col::Union{Symbol, AbstractVector{Symbol}} = :protgroup_id,
                           delim = ";")
    if elem_col === nothing
        endswith(string(list_col), "s") || error("elem_col not specified, autodetection failed")
        _elem_col = Symbol(replace(string(list_col), r"s$" => ""))
    else
        _elem_col = elem_col
    end
    exp_key_v, elem_v = expand_delim_column(df[!, list_col], 1:nrow(df), delim=delim)
    res = df[exp_key_v, key_col isa AbstractVector ? key_col : [key_col]]
    res[!, _elem_col] = elem_v
    return res
end

function unique_substrings(strings::AbstractVector{T}; delim = T(";")) where {T<:AbstractString}
    unique_substr = Set{T}()
    for str in skipmissing(strings)
        union!(unique_substr, split(str, delim))
    end
    return unique_substr
end

rejoin_unique_substrings(strings::AbstractVector{T}; delim = T(";"), joindelim = T(";")) where {T<:AbstractString} =
    join(sort(collect(unique_substrings(strings, delim=delim))), joindelim)

function union_substrings(strings::AbstractVector{Union{T,Missing}}, splitter::T = T(";")) where {T<:AbstractString}
    unique_substr = mapreduce(str -> unique(split(str::T, splitter)),
                              union!, Set{T}(), skipmissing(strings))
    join(sort(collect(unique_substr)), splitter)
end

function nmatches(r::Regex, s::AbstractString)
    n = 0
    for x in eachmatch(r, s)
        n+=1
    end
    return n
end

end
