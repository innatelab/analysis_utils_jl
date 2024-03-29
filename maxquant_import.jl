module MaxQuant

using DataFrames, CategoricalArrays

const MSImport = Main.MSImport
const FrameUtils = Main.FrameUtils

const MaxQuant_NAs = ["", "NA", "NaN", "n. def.", "n.def."]
const SILAC_mstags = ["L", "M", "H", "Sum"]

function fix_ratio_column_data(x::AbstractVector, name::String)
    if !(nonmissingtype(eltype(x)) <: Number)
        if !occursin("type", string(name))
            #@show coln eltype(x) x[1:10]
            error("Ratio column $coln is of type $(eltype(x)), expected Number")
        end
        return categorical(x)
    else
        return zero2missing(x, name)
    end
end

fix_ident_type_column_data(x::AbstractVector, name::Union{String, Nothing} = nothing) =
    levels!(categorical(ifelse.(isempty.(coalesce.(x, "")), missing, x)),
            ["By matching", "By MS/MS"])

function ColumnGroupDef(name::Symbol,
    columns::Union{AbstractVector, Nothing};
    required::Bool=false,
    update_coldata::Union{Function, Nothing}=nothing
)
    rename_regexs = Vector{Pair{Regex, SubstitutionString{String}}}()
    for (col, newcol) in columns
        if col isa Regex
            subst = col => newcol
        else
            subst = Regex(string("^", col, "\$")) => newcol
        end
        push!(rename_regexs, subst)
    end
    return MSImport.ColumnGroupDef(name, required, rename_regexs, update_coldata)
end

# data column groups in proteinGroups.txt
const MQ_ProteinGroups_ColumnGroupDefs = [
    ColumnGroupDef(:protgroup, [
        "id" => "protgroup_id",
        "Protein IDs" => "protein_acs", "Majority protein IDs" => "majority_protein_acs",
        "Gene names" => "gene_names", "Protein names" => "protein_names",
        "Fasta headers" => "fasta_headers",
        "Number of proteins" => "n_proteins",
        "Score" => "ident_score", "Q-value" => "q_value",
        "Sequence length" => "seqlen", "Sequence lengths" => "seqlens",
        "Mol. weight [kDa]" => "mol_weight_kDa",
        "Peptides" => "n_peptides",
        "Razor + unique peptides" => "n_peptides_unique_razor",
        "Unique peptides" => "n_peptides_unique",
        "Sequence coverage [%]" => "seqcov",
        "Unique + razor sequence coverage [%]" => "seqcov_unique_razor",
        "Unique sequence coverage [%]" => "seqcov_unique",
        "MS/MS count" => "n_ms2",
        r"(?:Potential )?contaminant" => "is_contaminant",
        "Reverse" => "is_reverse",
        "Only identified by site" => "is_site_only_ident",
    ], required=true),
    ColumnGroupDef(:intensity, [
        r"^Intensity\s([LMH])(?:\s(.+)|)$" => s"Intensity.\2|\1",
        r"^Intensity(?:\s(.+)|)$" => s"Intensity.\1"],
        update_coldata=MSImport.zero2missing),
    ColumnGroupDef(:LFQ, [
        r"^LFQ [Ii]ntensity\s([LMH])(?:\s(.+)|)$" => s"LFQ_Intensity.\2|\1",
        r"^LFQ [Ii]ntensity(?:\s(.+)|)$" => s"LFQ_Intensity.\1"],
        update_coldata=MSImport.zero2missing),
    ColumnGroupDef(:iBAQ, [
        r"^iBAQ\s([LMH])(?:\s(.+)|)$" => s"iBAQ.\2|\1",
        r"^iBAQ(?:\s(.+)|)" => s"iBAQ.\1"],
        update_coldata=MSImport.zero2missing),
    ColumnGroupDef(:ratio, [
        r"^iBAQ\s([LMH])(?:\s(.+)|)$" => s"iBAQ.\2|\1",
        r"^iBAQ(?:\s(.+)|)" => s"iBAQ.\1"],
        update_coldata=fix_ratio_column_data),
    ColumnGroupDef(:ident_type, [
        r"^Identification\stype\s" => s"ident_type."],
        update_coldata=fix_ident_type_column_data),
    ]
const MQ_ProteinGroups_ColumnGroupDef_Map =
    Dict(cg.name => cg for cg in MQ_ProteinGroups_ColumnGroupDefs)

read_protgroups(folder_path;
                file_name::AbstractString = "proteinGroups.txt",
                kwargs...
) = MSImport.read_wide_table(joinpath(folder_path, file_name), MQ_ProteinGroups_ColumnGroupDef_Map;
        missingstring=MaxQuant_NAs, truestrings=["+"], falsestrings=[""], strict=true,
        kwargs...)

parse_ratio_colnames(colnames::AbstractVector{Symbol}) =
    match.(r"^Ratio\.(?<mstag_nom>[^ ]+)/(?<mstag_denom>[^ ]+)(?:(?:\s(?<measure>[Nn]ormalized|[Vv]ariability \[\%\]|[Cc]ount|iso-count|type))?(?:\s(?<sample>.+))?)?$",
           string.(colnames))

function ratio_colinfo(colnames::AbstractVector{Symbol};
                       measures::Union{AbstractVector{String}, Nothing}=nothing,
                       samples::Union{AbstractVector{String}, Nothing}=nothing,
                       mstags::Union{AbstractVector{String}, Nothing}=nothing)
    parse_chunks = parse_ratio_colnames(colnames)
    df = DataFrame(colname = [Symbol(m.match) for m in parse_chunks],
                   mqmeasure = categorical([m[:measure] === nothing ? missing : m[:measure] for m in parse_chunks]),
                   mstag_nom = categorical(getindex.(parse_chunks, :mstag_nom)),
                   mstag_denom = categorical(getindex.(parse_chunks, :mstag_denom)),
                   sample = categorical([m[:sample] !== nothing ? m[:sample] : missing for m in parse_chunks]),
                   )
    (measures !== nothing) && levels!(df.mqmeasure, measures)
    (samples !== nothing) && levels!(df.sample, samples)
    if mstags !== nothing
        levels!(df.mstag_nom, mstags)
        levels!(df.mstag_denom, mstags)
    end
    return df
end

# Peptide Evidence Match (PEM) types, from more to less confident
const MaxQuant_IdentTypes = ["ISO-MSMS", "MULTI-MSMS", "MSMS", "MULTI-SECPEP", "MULTI-MATCH", "MULTI-MATCH-MSMS"]

const MQ_Evidence_ColumnGroupDefs = [
    ColumnGroupDef(:mschannel, [
        "Raw file" => "rawfile",
        "Experiment" => "msexperiment",
        "Fraction" => "msfraction",
        "Labeling State" => "label_state",
    ], required=true),
    ColumnGroupDef(:evidence, [
        "id" => "evidence_id",
        "Intensity" => "intensity",
        "Score" => "psm_score", "Delta score" => "psm_score_delta",
        "Max intensity m/z 0" => "max_intensity_mz0",
        # Match == Match Between Runs (MBR)
        "Type" => "ident_type",
        "Match time difference" => "mbr_time_diff", "Match m/z difference" => "mbr_mz_diff",
        "Match q-value" => "mbr_qvalue", "Match score" => "mbr_score",
        # PSM = peptide-to-spectrum match
        "PIF" => "PIF", "PEP" => "psm_pvalue",
        "Combinatorics" => "ptm_ncombs", # belong to evidence because for MBR it's 0
    ], required=true),
    ColumnGroupDef(:ret_time, [
        "Retention time" => "rt",
        "Retention length" => "rt_len",
        "Calibrated retention time" => "calib_rt",
        "Calibrated retention time start" => "calib_rt_start",
        "Calibrated retention time finish" => "calib_rt_finish",
        "Retention time calibration" => "rt_calib",
    ]),
    ColumnGroupDef(:peptide, [
        "Peptide ID" => "peptide_id",
        "Sequence" => "peptide_seq", "Length" => "peptide_len",
        "Missed cleavages" => "n_miscleaved",
        "Proteins" => "peptide_protein_acs",
        "Leading proteins" => "lead_protein_acs",
        "Leading razor protein" => "lead_razor_protein_ac",
        "Gene [Nn]ames" => "gene_names",
        "Protein [Nn]ames" => "protein_names",
        "Fasta headers" => "fasta_headers",
        "Reverse" => "is_reverse",
        "Potential contaminant" => "is_contaminant", "Contaminant" => "is_contaminant",
        "Protein group IDs" => "protgroup_ids", # IDs between different folders do not match
    ]),
    ColumnGroupDef(:pepmodstate, [
        "Peptide ID" => "peptide_id",
        "Mod. peptide ID" => "pepmod_id_mq", # groups peptides with the same set of PTMs but different localizations
        "Modifications" => "modifs",
        "Modified sequence" => "pepmod_seq",
        "Charge" => "charge",
    ]),
    ColumnGroupDef(:ms, [
        "m/z" => "mz", "MS/MS m/z" => "ms2_mz",
        "Mass [Ee]rror \\[ppm\\]" => "mass_error_ppm", "Mass [Ee]rror \\[Da\\]" => "mass_error_da",
        "Uncalibrated [Mm]ass [Ee]rror \\[ppm\\]" => "uncalib_mass_error_ppm", "Uncalibrated [Mm]ass [Ee]rror \\[Da\\]" => "uncalib_mass_error_da",
        "Uncalibrated - Calibrated m/z \\[Da\\]" => "delta_mz_da", "Uncalibrated - Calibrated m/z \\[ppm\\]" => "delta_mz_ppm",
        "MS/MS IDs" => "ms2_ids", "Best MS/MS" => "ms2_best_id", "AIF MS/MS IDs" => "ms2_aif_ids",
        "MS/MS [Cc]ount" => "n_ms2", "MS/MS [Ss]can [Nn]umber" => "ms2_scanid",
        "Mass" => "mass_da",
        "Resolution" => "resolution",
        "Number of data points" => "n_datapoints", "Number of scans" => "n_scans", "Number of isotopic peaks" => "n_isopeaks",
        "Fraction of total spectrum" => "frac_total_spectrum", "Base peak fraction" => "frac_base_peak",
    ]),
    # FIXME move to peptides.txt groups
    ColumnGroupDef(:aacount, [
        r"^(\w) [Cc]ount$" => s"aacount_\1",
    ]),
    ColumnGroupDef(:pepmod_ptms, [
        r"(.+)\s\((.+)\)\ssite IDs" => s"ptm_ids.\1",
        r"(.+)\s\((.+)\)\sProbabilities" => s"ptm_locprob_seq.\1",
        r"(.+)\s\((.+)\)\sScore Diffs" => s"ptm_scorediff_seq.\1",
    ]),
    ColumnGroupDef(:intensity, [
        r"Intensity (.+)" => s"intensity.\1",
    ]),
]
const MQ_Evidence_ColumnGroupDef_Map =
    Dict(cg.name => cg for cg in MQ_Evidence_ColumnGroupDefs)

struct Evidence
    peptides::DataFrame
    pepmods::DataFrame
    pepmodstates::DataFrame
    peptide2protgroup::DataFrame
    peptide2protein::DataFrame
    observations::DataFrame
end

function read_evidence(folder_path;
    file_name::AbstractString = "evidence.txt",
    msruns::Union{AbstractDataFrame, Nothing}=nothing,
    verbose::Bool=false,
    kwargs...
)
    evidence_df, colgroups = MSImport.read_wide_table(joinpath(folder_path, file_name), MQ_Evidence_ColumnGroupDef_Map;
            missingstring=MaxQuant_NAs, truestrings=["+"], falsestrings=[""], strict=true,
            verbose=verbose, kwargs...)
    if !isnothing(msruns)
        MSImport.update_msruns!(evidence_df, colgroups[:mschannel], msruns, orig_column_suffix="_mq", verbose=verbose)
        if hasproperty(msruns, :is_skipped)
            verbose && @info "Skipping $(sum(msruns.is_skipped)) raw files: $(join(filter(r -> r.is_skipped, msruns).rawfile, ", "))"
            filter!(r -> !coalesce(r.is_skipped, false), evidence_df)
        end
    end
    if !any(ren -> ren[2] == "msrun", colgroups[:mschannel])
        if any(!ismissing, evidence_df.msfraction)
            verbose && @info "Setting msrun = msexperiment X msfraction"
            evidence_df.msrun = evidence_df.msexperiment .* "_" .*
                (nonmissingtype(eltype(evidence_df.msfraction)) <: Number ?
                "F" .* string.(evidence_df.msfraction) :
                evidence_df.msfraction)
        else
            verbose && @info "No MS fractions detected. Setting msrun = msexperiment"
            evidence_df.msrun = evidence_df.msexperiment
        end
        push!(colgroups[:mschannel], nothing => "msrun")
    end
    FrameUtils.categorical!(evidence_df, :ident_type, levels=MaxQuant_IdentTypes, ordered=true)
    catcols = [newcol for (_, newcol) in colgroups[:mschannel]
                if !(nonmissingtype(eltype(evidence_df[!, newcol])) <: Number)]
    if !isempty(catcols)
        verbose && @info "Converting MS run columns to categorical type: $(join(catcols, ", "))"
        FrameUtils.categorical!(evidence_df, catcols)
    end
    if haskey(colgroups, :pepmodstate)
        MSImport.update_pepmodstates!(evidence_df, colgroups[:pepmodstate], verbose=verbose)
    end
    return evidence_df, colgroups
end

function pivot_ptms_longer(df::AbstractDataFrame, idcols::Union{AbstractVector, Symbol})
    res = FrameUtils.pivot_longer(df, idcols,
            measure_vars_regex=r"^(?<value>[^.]+)\.(?<var>.+)$",
            var_col=:ptm_type)
    return res
end

#= old Evidence code
    file_path = joinpath(folder_path, file_name)
    @info "Reading $file_path..."
    evid_df = CSV.read(file_path, delim='\t', missingstring=MaxQuant_NAs,
                       limit=limit_rows,
                       truestrings=["+"], falsestrings=[""], strict=true)
    @info "$(size(evid_df, 1)) evidence record(s), $(size(evid_df, 2)) column(s)"
    cols_map = copy(MQ_Evidence_ColumnMap)
    cols_mask = first.(cols_map) .∈ Ref(Set(names(evid_df)))
    col_info = Dict{Symbol, Vector{Symbol}}()
    if !all(cols_mask)
        @warn "Column(s) $(join(string.("\"", first.(cols_map[.!cols_mask]), "\""), ", ")) not found in the file"
        cols_map = cols_map[cols_mask]
    end
    evid_cols_mask = names(evid_df) .∈ Ref(first.(cols_map))
    if !all(evid_cols_mask)
        @warn "File column(s) $(join(string.("\"", names(evid_df)[.!evid_cols_mask], "\""), ", ")) ignored"
    end
    @info "Importing evidences..."
    obs_df = select(evid_df, first.(cols_map), copycols=false)
    rename!(obs_df, cols_map...)
    categorical!(obs_df, [:ident_type, :msrun, :raw_file])

    peptide_cols = [:peptide_id, :peptide_seq, :peptide_len, :n_miscleaved,
                    :protein_acs, :lead_protein_acs, :lead_razor_protein_ac, :protgroup_ids,
                    :is_reverse, :is_contaminant]
    pepmod_cols = [:pepmod_id, :modifs, :pepmod_seq, :peptide_id]
    pepmodstate_cols = [:pepmod_id, :peptide_id, :charge]
    peptides_df = sort!(unique!(select(obs_df, peptide_cols)), :peptide_id)
    peptides_df.is_reverse = coalesce.(peptides_df.is_reverse, false)
    peptides_df.is_contaminant = coalesce.(peptides_df.is_contaminant, false)
    @info "$(size(peptides_df, 1)) peptide(s) imported..."
    pepmods_df = sort!(unique!(select(obs_df, pepmod_cols)), :pepmod_id)
    @info "$(size(pepmods_df, 1)) pepmods(s) imported..."
    pepmodstates_df = sort!(unique!(select(obs_df, pepmodstate_cols)), [:pepmod_id, :charge])
    pepmodstates_df.pepmodstate_id = collect(1:size(pepmodstates_df, 1))
    @info "$(size(pepmodstates_df, 1)) pepmod state(s) imported..."

    # add pepmodstate_id to observations
    obs_df = join(obs_df, select(pepmodstates_df, Not(:peptide_id)),
                  kind=:inner, on=[:pepmod_id, :charge])
    # remove peptide/pepmod/pepmodstate columns from observations
    select!(obs_df, Not(setdiff(union(peptide_cols, pepmod_cols, pepmodstate_cols),
                                [:pepmodstate_id, :pepmod_id, :peptide_id])))

    @info "Assembling peptide-to-protgroup map..."
    pep2pg_df = DelimDataUtils.expand_delim_column(peptides_df, list_col=:protgroup_ids, key_col=:peptide_id)
    pep2pg_df.protgroup_id = parse.(Int, pep2pg_df.protgroup_id)
    @info "Assembling peptide-to-protein map..."
    pep2protac_df = DelimDataUtils.expand_delim_column(peptides_df, list_col=:protein_acs, key_col=:peptide_id)
    pep2protac_df[!, :is_lead] .= false
    pep2protac_df[!, :is_razor] .= false
    pep2leadprotac_df = DelimDataUtils.expand_delim_column(peptides_df, list_col=:lead_protein_acs, key_col=:peptide_id, elem_col=:protein_ac)
    pep2leadprotac_df[!, :is_lead] .= true
    pep2leadprotac_df[!, :is_razor] .= false
    pep2leadrazorprotac_df = DelimDataUtils.expand_delim_column(peptides_df, list_col=:lead_razor_protein_ac, key_col=:peptide_id, elem_col=:protein_ac)
    pep2leadrazorprotac_df[!, :is_lead] .= true
    pep2leadrazorprotac_df[!, :is_razor] .= true
    pep2ac_df = vcat(pep2protac_df, pep2leadprotac_df, pep2leadrazorprotac_df)
    pep2ac_df = by(pep2ac_df, [:peptide_id, :protein_ac],
                   is_lead = :is_lead => any,
                   is_razor = :is_razor => any)
    return Evidence(peptides_df, pepmods_df, pepmodstates_df,
                    pep2pg_df, pep2ac_df, obs_df)
end
=#
end
