module MSImport

const DelimDataUtils = Main.DelimDataUtils
const FrameUtils = Main.FrameUtils

using CSV, DataFrames, CategoricalArrays

zero2missing(::Missing) = missing

zero2missing(x::Real) =
    ifelse(!ismissing(x) && x > 0, float.(x), missing)

zero2missing(x::AbstractVector{<:Real}, name::Union{String,Nothing} = nothing) =
    zero2missing.(x)
zero2missing(x::AbstractVector{Union{T, Missing}}, name::Union{String,Nothing} = nothing) where T=
    zero2missing.(x)

"""
Definition of a group of columns in Spectronaut output table:
the regular expression pattern, the renaming rule and
the rule to fix the data.
"""
struct ColumnGroupDef
    name::Symbol
    rename_regexs::Vector{Pair{Regex, SubstitutionString{String}}}
    update_report!::Union{Function, Nothing}
    update_coldata::Union{Function, Nothing}
end

function extractnames(df::DataFrame, colgroup::ColumnGroupDef)
    res = Vector{Pair{String, String}}()
    for (col_regex, col_subst) in colgroup.rename_regexs
        for col in names(df)
            if occursin(col_regex, col)
                push!(res, col => replace(col, col_regex => col_subst))
            end
        end
    end
    return res
end

function extract(df::DataFrame, colgroup::ColumnGroupDef;
                 verbose::Bool=false)
    cols = extractnames(df, colgroup)
    cols_ext = convert(Vector{Pair{Union{String, Nothing}, String}}, cols)
    cg_df = df[!, first.(cols)]
    if !isempty(cols)
        rename!(cg_df, cols...)
        if colgroup.update_coldata !== nothing
            for coln in names(cg_df)
                verbose && @info "Postprocessing $coln ($(eltype(cg_df[!, coln])))"
                cg_df[!, coln] = colgroup.update_coldata(cg_df[!, coln], coln)
            end
        end
        if colgroup.update_report! !== nothing
            colgroup.update_report!(cg_df, cols_ext, df)
        end
    end
    return cg_df, cols_ext
end

function update_protgroup_report!(df::DataFrame, cols::AbstractVector, orig_df::AbstractDataFrame)
    df.protgroup_id = FrameUtils.indexunique(df.protgroup_sn_id)
    # fix majority_protein_acs if only _alt version is available
    if !hasproperty(df, :majority_protein_acs) && hasproperty(df, :majority_protein_acs_alt)
        rename!(df, :majority_protein_acs_alt => :majority_protein_acs)
        deleteat!(cols, findfirst(x => last(x) == "majority_protein_acs", cols))
        ix = findfirst(x => last(x) == "majority_protein_acs_alt", cols)
        cols[ix] = first(cols[ix]) => "majority_protein_acs"
    end
    push!(cols, nothing => "protgroup_id")
    return df
end

function update_peptide_report!(df::DataFrame, cols::AbstractVector, orig_df::AbstractDataFrame)
    if !hasproperty(orig_df, :peptide_seq)
        pepmod_seqs = hasproperty(orig_df, :pepmod_seq) ? orig_df.pepmod_seq : orig_df[!, "EG.ModifiedSequence"]
        df.peptide_seq = replace.(pepmod_seqs, Ref(r"\[[^[]+\]" => ""))
        push!(cols, nothing => "peptide_seq")
    end
    df.peptide_id = FrameUtils.indexunique(df.peptide_seq)
    push!(cols, nothing => "peptide_id")
    return df
end

function update_pepmodstate_report!(df::DataFrame, cols::AbstractVector, orig_df::AbstractDataFrame)
    df.pepmod_id = FrameUtils.indexunique(df.pepmod_seq)
    df.pepmodstate_id = FrameUtils.indexunique(df.pepmodstate_seq)
    matches = match.(Ref(r"\.(\d+)$"), df.pepmodstate_seq)
    df.charge = parse.(Int, getindex.(matches, 1))
    push!(cols, nothing => "pepmodstate_id")
    push!(cols, nothing => "charge")
    return df
end

function read_wide_table(file_path::AbstractString,
        colgroups_defs::AbstractDict{Symbol, ColumnGroupDef};
        import_data::AbstractVector{Symbol} = Vector{Symbol}(),
        limit_rows::Union{Integer, Nothing} = nothing,
        verbose::Bool=false,
        csv_kwargs...,
)
    verbose && @info "Reading $file_path..."
    orig_df = CSV.read(file_path; limit=limit_rows, strict=true, csv_kwargs...)
    verbose && @info "$(size(orig_df, 1)) rows(s), $(size(orig_df, 2)) column(s)"
    res_df = orig_df[!, String[]]
    col_info = Dict{Symbol, Vector{Pair{Union{String, Nothing}, String}}}()
    for cgname in import_data
        if !haskey(colgroups_defs, cgname)
            @warn "Unknown column group: $cgname"
            continue
        end
        cgdef = colgroups_defs[cgname]
        cg_df, cg_colmap = extract(orig_df, cgdef, verbose=verbose)
        if isempty(cg_colmap)
            @warn("No $(cgdef.name) columns found in data")
        else
            verbose && @info "Appending $(cgdef.name) columns to the report: $(size(cg_df, 2)) column(s) added"
            res_df = hcat(res_df, cg_df)
            col_info[cgname] = cg_colmap
        end
    end
    for col in names(res_df)
        if nonmissingtype(eltype(res_df[!, col])) === Bool
            verbose && @info "Fixing $col boolean column"
            res_df[!, col] = coalesce.(res_df[!, col], false)
        end
    end

    return res_df, col_info
end

include(joinpath(@__DIR__, "spectronaut_import.jl"))
include(joinpath(@__DIR__, "maxquant_import.jl"))

end
