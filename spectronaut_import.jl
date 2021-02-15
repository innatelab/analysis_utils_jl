module Spectronaut

using DataFrames

const MSImport = Main.MSImport
const FrameUtils = Main.FrameUtils

const Spectronaut_NAs = ["", "Filtered"]
const Spectronaut_Scopes = ["PG", "PEP", "EG", "R"]

function ColumnGroupDef(name::Symbol,
    columns::Union{AbstractVector, Nothing};
    update_report!::Union{Function, Nothing}=nothing,
    update_coldata::Union{Function, Nothing}=nothing
)
    rename_regexs = Vector{Pair{Regex, SubstitutionString{String}}}()
    for (col, newcol) in columns
        ismetric = occursin(r"^\\\.", col)
        hasscope = occursin(Regex(string("^\\.?(?:", join(Spectronaut_Scopes, "|"), ")\\.")), col)
        if hasscope
            if ismetric
                subst = Regex(string("^(.+)", col, "\$")) =>
                        SubstitutionString(string(newcol, ".\\1"))
            else
                subst = Regex(string("^", col, "\$")) => newcol
            end
        else
            if ismetric
            subst = Regex(string("^(.+)\\.(", join(Spectronaut_Scopes, "|"),")", col, "\$")) =>
                    SubstitutionString(string("\\2_", newcol, ".\\1"))
            else
            subst = Regex(string("^(", join(Spectronaut_Scopes, "|"),")\\.", col, "\$")) =>
                    SubstitutionString(string("\\1_", newcol))
            end
        end
        push!(rename_regexs, subst)
    end
    return MSImport.ColumnGroupDef(name, rename_regexs, update_report!, update_coldata)
end

# data column groups in proteins report
const Spectronaut_ColumnGroupDefs = [
    ColumnGroupDef(:quantity, ["\\.Quantity" => "intensity",
                               "\\.TotalQuantity \\(Settings\\)" => "intensity",
                               "\\.NormalizationFactor" => "normfactor",
                               "\\.Qvalue" => "qvalue"],
                   update_coldata=MSImport.zero2missing),
    ColumnGroupDef(:properties, ["IsSingleHit" => "is_singlehit"]),
    ColumnGroupDef(:ident_stats, ["\\.RunEvidenceCount" => "nevidences",
                                "\\.NrOfStrippedSequencesMeasured" => "npeptides_quanted",
                                "\\.NrOfStrippedSequencesIdentified" => "npeptides_idented",
                                "\\.NrOfModifiedSequencesMeasured" => "npepmods_quanted",
                                "\\.NrOfModifiedSequencesIdentified" => "npepmods_idented",
                                "\\.NrOfPrecursorsIdentified" => "npepmodstates_idented"],
                   update_coldata=MSImport.zero2missing),
    ColumnGroupDef(:protgroup, ["PG.ProteinGroups" => "protgroup_sn_id",
                                "PG.ProteinAccessions" => "majority_protein_acs",
                                "PG.UniProtIds" => "majority_protein_acs_alt",
                                "PG.Genes" => "gene_names", "PG.ProteinNames" => "protein_names",
                                "PG.ProteinDescription" => "protein_descriptions",
                                "PG.Organisms" => "organisms",
                                "PG.Qvalue" => "qvalue"],
                   update_report! = MSImport.update_protgroup_report!),
    ColumnGroupDef(:pepmodstate, ["EG.PrecursorId" => "pepmodstate_seq",
                                  "EG.ModifiedSequence" => "pepmod_seq",
                                  "EG.IsDecoy" => "is_decoy"],
                   update_report! = MSImport.update_pepmodstate_report!),
    ColumnGroupDef(:peptide, ["PEP.PeptidePosition" => "peptide_poses",
                              "PEP.IsProteinGroupSpecific" => "is_specific_peptide",
                              "PEP.AllOccurringProteinAccessions" => "peptide_protein_acs"],
                   update_report! = MSImport.update_peptide_report!),
    ColumnGroupDef(:msrun, ["R.FileName" => "rawfile",
                            "R.Replicate" => "msrun_ix"]),
    ColumnGroupDef(:metrics, ["\\.PTMLocalizationProbabilities" => "locprob_seq",
                              ])
]

const Spectronaut_ColumnGroupDef_Map =
    Dict(cg.name => cg for cg in Spectronaut_ColumnGroupDefs)

read_report(file_path::AbstractString; kwargs...) =
    MSImport.read_wide_table(file_path, Spectronaut_ColumnGroupDef_Map;
        missingstrings = Spectronaut_NAs,
        truestrings=["+", "True"], falsestrings=["", "False"],
        kwargs...)

const Metric_Colname_Regex = Regex(string("^(?:(?<scope>",
    join(Spectronaut_Scopes, "|"), ")_)?",
    "(?<measure>[^.]+)\\.\\[(?<rawfile_ix>\\d+)\\]\\s(?<rawfile>.+)?\$"))

const METRIC4PIVOT_REGEX = r"^(?<value>[^.]+)\.\[(?<var>\d+)\]\s(?<rawfile>.+)?$"

parse_metrics_colnames(colnames::AbstractVector{String}) =
    match.(Metric_Colname_Regex, colnames)

function metrics_colinfo(colnames::AbstractVector{String};
                         measures::Union{AbstractVector{String}, Nothing}=nothing,
                         rawfiles::Union{AbstractVector{String}, Nothing}=nothing)
    parse_chunks = parse_metrics_colnames(colnames)
    df = DataFrame(colname = [String(m.match) for m in parse_chunks],
                scope = levels!(categorical(getindex.(parse_chunks, :scope)), Spectronaut_Scopes),
                measure = categorical(getindex.(parse_chunks, :measure)),
                rawfile_ix = parse.(Int, getindex.(parse_chunks, :rawfile_ix)),
                rawfile = categorical(getindex.(parse_chunks, :rawfile)))
    (measures !== nothing) && levels!(df.measure, measures)
    (rawfiles !== nothing) && levels!(df.rawfile, rawfiles)
    return df
end

function pivot_longer(df::AbstractDataFrame, idcols::Union{AbstractVector, Symbol})
    res = FrameUtils.pivot_longer(df, idcols,
            measure_vars_regex=METRIC4PIVOT_REGEX,
            var_col=:rawfile_ix, value_col=:metric)
    res.rawfile_ix = parse.(Int, get.(res.rawfile_ix))
    return res
end

end
