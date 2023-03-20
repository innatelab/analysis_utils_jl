module GMT

using DataFrames, CategoricalArrays

"""
    read(ac_type::Type{T}, filename) where T -> Dict{String, Set{T}}

Read GMT (Gene Matrix Transposed) file.

### See Also
[GMT format description](http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29).
"""
function read(ac_type::Type{T}, filename;
              id_col::Symbol = :id,
              src_col::Symbol = :src,
              name_col::Symbol = :name,
              descr_col::Symbol = :descr,
              fix_src::Bool = true) where T <: Union{String, Integer}
    gmt_txt = readlines(filename)
    names = sizehint!(Vector{String}(), length(gmt_txt))
    srcs = similar(names)
    ids = similar(names)
    descrs = similar(names)
    coll = Dict{String, Set{T}}()
    for gmt_line in gmt_txt
        chunks = split(rstrip(gmt_line), '\t')
        info = chunks[1]
        descr = chunks[2]
        name, src, id = split(info, '%')
        if fix_src
            src = fix_annotation_src(src)
        end
        if src == "Reactome" && occursin(r"^\d+$", id)
            # fix Reactome IDs, otherwise plotly thinks they are numbers and messes heatmap
            id = string("R-HSA-", id)
        end
        push!(names, name)
        push!(descrs, descr)
        push!(srcs, src)
        push!(ids, id)
        if T <: Integer
            coll[id] = Set(parse(Int, chunks[i]) for i in 3:length(chunks))
        else
            coll[id] = Set(chunks[3:end])
        end
    end

    return DataFrame(id_col => ids,
                     src_col => srcs,
                     name_col => names,
                     descr_col => descrs),
           coll
end

"""
Fix annotation source names.
"""
function fix_annotation_src(src::AbstractString)
    for (regex, fixedname) in [
        r"GObp"i => "GO_BP", r"GOcc"i => "GO_CC", r"GOmf"i => "GO_MF",
        r"Reactome Database ID Release|REACTOME"i => "Reactome",
        r"PANTHER PATHWAY"i => "Panther_Pathway",
        r"PATHWAY INTERACTION"i => "NCI_Nature",
        r"WikiPathways"i => "WikiPathways"]
        occursin(regex, src) && return fixedname
    end
    return src
end

fix_annotation_src(src::AbstractCategoricalVector) =
    levels!(src, levels(fix_annotation_src(levels(src))))

"""
    group_by_src(annot_info_df, annot_sets;
                 src_col::Symbol = :src) -> Dict{Symbol, Dict{String, <:Set}}

Split `annot_sets` into set collections by the source of annotation,
which is taken from `src_col` of `annot_info_df`.
"""
function group_by_src(annot_info_df::AbstractDataFrame,
                      annot_sets::AbstractDict{String, <:Set};
                      src_col::Union{Symbol, String} = :src,
                      verbose::Bool = false)
    colls = Dict{Symbol, typeof(annot_sets)}()
    for coll_subdf in groupby(annot_info_df, src_col, sort=false)
        coll_id = Symbol(get(coll_subdf[1, src_col]))
        verbose && @info "Processing $coll_id annotation collection..."
        colls[coll_id] = Dict(set_id => annot_sets[set_id]
                              for set_id in coll_subdf.id)
        verbose && @info "  $(length(colls)) set(s) processed"
    end
    return colls
end

end
