module PhosphoSitePlus

using DataFrames, CSV, CodecZlib

const AnnotationFiles = [
    :KinaseSubstrate => "Kinase_Substrate_Dataset",
    :DiseaseAssoc => "Disease-associated_sites",
    :PhosphoSites => "Phosphorylation_site_dataset",
    :UbiSites => "Ubiquitination_site_dataset",
    :RegSites => "Regulatory_sites"
]

const GenericColumnMap = Dict(
    :GENE => :gene_name, :SUB_GENE => :gene_name,
    :PROTEIN => :protein_name, :KINASE => :kinase_protein_name, :SUBSTRATE => :protein_name,
    :ACC_ID => :protein_ac, :KIN_ACC_ID => :kinase_protein_ac, :SUB_ACC_ID => :protein_ac,
    :GENE_ID => :entrez_id, :SUB_GENE_ID=> :entrez_id,
    :ORGANISM => :organism, :KIN_ORGANISM => :kinase_organism, :SUB_ORGANISM => :organism,
    Symbol("SITE_+/-7_AA") => :flanking_AAs,
    :DOMAIN => :domain,
    :DISEASE => :disease
)

const SpecificColumnMaps = Dict(
    :KinaseSubstrate => Dict(
        :GENE => :kinase_gene_name,
    ),
    :RegSites => Dict(
        :ON_FUNCTION => :reg_function,
        :ON_PROT_INTERACT => :reg_prot_iactions,
        :ON_OTHER_INTERACT => :reg_other_iactions,
        :PMIDs => :reg_pubmed_ids,
    ),
    :DiseaseAssoc => Dict(
        :PMIDs => :disease_pubmed_ids,
    )
)

collapsevals(vals; delim="; ") =
    all(ismissing, vals) ? missing : join(skipmissing(vals), delim)

function read_annotations(folder::AbstractString; verbose::Bool=false)
    res = Dict(begin
        verbose && @info "Reading $dsname from $filename..."
        if isfile(joinpath(folder, filename))
            df = CSV.read(joinpath(folder, filename), header=4, threaded=false)
        elseif isfile(joinpath(folder, filename * ".gz"))
            df = open(GzipDecompressorStream, joinpath(folder, filename * ".gz"), "r") do io
                CSV.read(io, header=4, threaded=false)
            end
        end
        verbose && @info "  $(nrow(df)) $dsname annotations read"
        custom_colmap = get(SpecificColumnMaps, dsname, nothing)
        for col in propertynames(df)
            newcol = nothing
            isnothing(custom_colmap) || (newcol = get(custom_colmap, col, nothing))
            isnothing(newcol) && (newcol = get(GenericColumnMap, col, nothing))
            if !isnothing(newcol)
                verbose && @info "    renaming $col => $newcol"
                rename!(df, col => newcol)
            end
        end
        if hasproperty(df, :MOD_RSD)
            ptm_matches = match.(Ref(r"(\w)(\d+)-(\w+)"), df.MOD_RSD)
            df.ptm_AA = string.(getindex.(ptm_matches, 1))
            df.ptm_pos = parse.(Int, getindex.(ptm_matches, 2))
            df.ptm_code = string.(getindex.(ptm_matches, 3))
        elseif hasproperty(df, :SUB_MOD_RSD)
            ptm_matches = match.(Ref(r"(\w)(\d+)"), df.SUB_MOD_RSD)
            df.ptm_AA = string.(getindex.(ptm_matches, 1))
            df.ptm_pos = parse.(Int, getindex.(ptm_matches, 2))
        end
        dsname => df
    end for (dsname, filename) in AnnotationFiles)
    return res
end

function combine_annotations(annot_dfs::AbstractDict)
    phospho_annots_df = select(annot_dfs[:PhosphoSites], [:protein_ac, :ptm_pos, :ptm_AA, :ptm_code, :domain])
    ubi_annots_df = select(annot_dfs[:UbiSites], [:protein_ac, :ptm_pos, :ptm_AA, :ptm_code, :domain])
    phubi_annots_df = vcat(phospho_annots_df, ubi_annots_df)
    phubi_annots_df[!, :ptm_is_known] .= true
    ks_annots_df = combine(groupby(annot_dfs[:KinaseSubstrate], [:protein_ac, :ptm_pos, :ptm_AA]),
                           :kinase_gene_name => collapsevals => :kinase_gene_names)
    regsite_annots_df = combine(groupby(annot_dfs[:RegSites], [:protein_ac, :ptm_pos, :ptm_AA, :ptm_code]),
        :reg_function => collapsevals => :reg_function,
        :reg_prot_iactions => collapsevals => :reg_prot_iactions,
        :reg_other_iactions => collapsevals => :reg_other_iactions,
        :reg_pubmed_ids => collapsevals => :reg_pubmed_ids)
    disease_annots_df = combine(groupby(annot_dfs[:DiseaseAssoc], [:protein_ac, :ptm_pos, :ptm_AA, :ptm_code]),
                                :disease => collapsevals => :diseases,
                                :disease_pubmed_ids => collapsevals => :diseases_pubmed_ids)
    annots_df = reduce((df1, df2) -> outerjoin(df1, df2, on=intersect(intersect([:protein_ac, :ptm_pos, :ptm_AA, :ptm_code],
                                                                      propertynames(df1), propertynames(df2)))),
                       [phubi_annots_df, ks_annots_df,
                        regsite_annots_df, disease_annots_df])
    return annots_df
end

end
