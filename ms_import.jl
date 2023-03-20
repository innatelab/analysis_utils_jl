module MSImport

using Core: Argument
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

ColumnRename = Union{Pair{String, String}, Pair{Nothing, String}}

"""
Definition of a group of columns in MS software output table:
the regular expression pattern, the renaming rule and
the rule to fix the data.
"""
struct ColumnGroupDef
    name::Symbol
    required::Bool
    rename_regexs::Vector{Pair{Regex, SubstitutionString{String}}}
    update_coldata::Union{Function, Nothing}
end

function extractnames(df::DataFrame, colgroup::ColumnGroupDef)
    res = Vector{ColumnRename}()
    for (col_regex, col_subst) in colgroup.rename_regexs
        for col in names(df)
            if occursin(col_regex, col)
                push!(res, col => replace(col, col_regex => col_subst))
            end
        end
    end
    return res
end

function update_msruns!(df::DataFrame, imported_msrun_cols::AbstractVector{<:ColumnRename},
                        msruns::AbstractDataFrame;
                        orig_column_suffix::String = "_orig",
                        verbose::Bool=false)
    verbose && @info "Updating imported data with user-supplied msrun info..."
    if !("rawfile" in names(msruns))
        error("rawfile column not found in user-supplied msruns")
    end
    if !any(rn -> rn[2] == "rawfile", imported_msrun_cols)
        error("rawfile column not found in imported data")
    end
    msruns = copy(msruns, copycols=false)
    msruns.isduplicated = nonunique(msruns, :rawfile)
    if any(msruns.isduplicated)
        error("$(sum(msruns.isduplicated)) duplicated rawfile(s) found in user-supplied msruns: $(join(msruns.rawfile[msruns.isduplicated], ", "))")
    end
    user_rawfile2ix = Dict(rawfile => i for (i, rawfile) in enumerate(msruns.rawfile))
    df2userix = Dict(rawfile => get(user_rawfile2ix, rawfile, 0) for rawfile in df.rawfile)
    missing_rawfiles = [rawfile for (rawfile, ix) in df2userix if ix==0]
    if !isempty(missing_rawfiles)
        @warn "  * $(length(missing_rawfiles)) imported rawfile(s) have no information in user-supplied msruns: $(join(missing_rawfiles, ", "))"
    end
    df_user_ixs = getindex.(Ref(df2userix), df.rawfile)
    nonmissing_mask = df_user_ixs .!= 0
    nonmissing_user_ixs = df_user_ixs[nonmissing_mask]

    for col in names(msruns)
        (col == "rawfile") && continue
        is_msrun_col = false
        for (i, rename) in enumerate(imported_msrun_cols)
            # rename the existing column with the same name
            if col == rename[2]
                verbose && @info "  * column $col exists in imported data, renaming to $(col)$(orig_column_suffix)"
                newcol = col*orig_column_suffix
                imported_msrun_cols[i] = rename[1] => newcol
                rename!(df, col => newcol)
                is_msrun_col = true
                break
            end
        end
        if !is_msrun_col && (col in names(df))
            @warn "  * user-specified msrun column $col already exists in imported data, but does not belong to msruns group, skipping"
        else
            verbose && @info "  * adding user-specified msrun column $col"
            if !isempty(missing_rawfiles)
                df[!, col] = missings(eltype(msruns[!, col]), nrow(df))
                df[nonmissing_mask, col] .= msruns[nonmissing_user_ixs, col]
            else
                df[!, col] = msruns[nonmissing_user_ixs, col]
            end
            push!(imported_msrun_cols, nothing => col)
        end
    end
end

function extract(df::DataFrame, colgroup::ColumnGroupDef; verbose::Bool=false)
    cols = extractnames(df, colgroup)
    cg_df = df[!, first.(cols)]
    if !isempty(cols)
        rename!(cg_df, cols...)
        if colgroup.update_coldata !== nothing
            for coln in names(cg_df)
                verbose && @info "  * postprocessing $coln ($(eltype(cg_df[!, coln])))"
                cg_df[!, coln] = colgroup.update_coldata(cg_df[!, coln], coln)
            end
        end
    end
    return cg_df, cols
end

function update_protgroup_report!(df::DataFrame, cols::AbstractVector, orig_df::AbstractDataFrame;
                                  verbose::Bool=false)
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
    df.pepmod_id = FrameUtils.indexunique(df.pepmod_seq, sort=false, start=0)
    df.pepmodstate_id = FrameUtils.indexunique(df.pepmodstate_seq, sort=false, start=0)
    matches = match.(Ref(r"\.(\d+)$"), df.pepmodstate_seq)
    df.charge = parse.(Int, getindex.(matches, 1))
    push(cols, nothing => "pepmod_id")
    push(cols, nothing => "pepmodstate_id")
    push(cols, nothing => "charge")
    return df
end

function update_pepmodstates!(df::DataFrame, cols::AbstractVector{<:MSImport.ColumnRename};
                              verbose::Bool=false)
    df.pepmod_id = FrameUtils.indexunique(collect(zip(df.peptide_id, df.pepmod_seq)), start=0, sort=true)
    push!(cols, nothing => "pepmod_id")

    pms_cols = ["pepmod_id"]
    if ("msfraction" in names(df)) && any(!ismissing, df.msfraction)
        verbose && @info "Using msfractions to index pepmodstates"
        push!(pms_cols, "msfraction")
        # FIXME src column
        push!(cols, nothing => "msfraction")
    end
    push!(pms_cols, "charge")
    dftmp = select(df, pms_cols, copycols=false)
    dftmp.pepmodstate_ck = collect(zip(eachcol(dftmp)...))
    df.pepmodstate_id = FrameUtils.indexunique(dftmp.pepmodstate_ck, start=0, sort=true)
    push!(cols, nothing => "pepmodstate_id")
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
    orig_df = CSV.read(file_path, DataFrame; limit=limit_rows, strict=true, csv_kwargs...)
    verbose && @info "$(size(orig_df, 1)) rows(s), $(size(orig_df, 2)) column(s)"
    res_df = orig_df[!, String[]]
    col_info = Dict{Symbol, Vector{ColumnRename}}()
    bad_cgnames = setdiff(import_data, keys(colgroups_defs))
    if !isempty(bad_cgnames)
        @warn "Unknown column groups: $(join(bad_cgnames, ","))"
    end
    for (cgname, cgdef) in colgroups_defs
        cgdef.required || (cgname ∈ import_data) || continue
        cg_df, cg_renames = extract(orig_df, cgdef, verbose=verbose)
        if isempty(cg_renames)
            @warn("  * no $(cgdef.name) columns found in data")
        else
            verbose && @info "  * appending $(cgdef.name) columns to the report: $(size(cg_df, 2)) column(s) added"
            res_df = hcat(res_df, select(cg_df, setdiff(names(cg_df), names(res_df))))
            col_info[cgname] = cg_renames
        end
    end
    if isempty(col_info)
        throw(ArgumentError("No import_data specified, expected a subset of $(join(keys(colgroups_defs) |> collect |> sort, ", "))"))
    end
    for col in names(res_df)
        if nonmissingtype(eltype(res_df[!, col])) === Bool
            verbose && @info "  * fixing $col boolean column"
            res_df[!, col] = coalesce.(res_df[!, col], false)
        end
    end

    return res_df, col_info
end

parse_msrun_colname(colname::Union{String, Symbol}) =
    match(r"^(?<measure>[^.]+)(?:\.(?<msrun>[^|]+))?(?:\|(?<mstag>.+))?$", string(colname))

function parse_msrun_colnames(colnames::AbstractVector{String};
                              measures::Union{AbstractVector{String}, Nothing}=nothing,
                              msruns::Union{AbstractVector{String}, Nothing}=nothing,
                              mstags::Union{AbstractVector{String}, Nothing}=nothing)
    parse_chunks = parse_msrun_colname.(colnames)
    df = DataFrame(colname = [Symbol(m.match) for m in parse_chunks],
                   measure = categorical(getindex.(parse_chunks, :measure)),
                   msrun = categorical([m[:msrun] !== nothing && !isempty(m[:msrun]) ?
                                       m[:msrun] : missing for m in parse_chunks]),
                   mstag = categorical([m[:mstag] !== nothing && !isempty(m[:mstag]) ?
                                       m[:mstag] : missing for m in parse_chunks]),
                   )
    (measures !== nothing) && levels!(df.measure, measures)
    (msruns !== nothing) && levels!(df.msrun, msruns)
    (mstags !== nothing) && levels!(df.mstag, mstags)
    return df
end

include(joinpath(@__DIR__, "spectronaut_import.jl"))
include(joinpath(@__DIR__, "maxquant_import.jl"))
include(joinpath(@__DIR__, "fragpipe_import.jl"))

end
