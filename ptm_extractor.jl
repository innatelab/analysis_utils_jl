module PTMExtractor

using DataFrames, Printf
using BioSequences, BioAlignments

const DelimDataUtils = Main.DelimDataUtils

# spectronaut-specific extraction of protein(ac)-to-peptide matches
function peptide_matches(peptides_df::AbstractDataFrame;
                         protein_acs_col::Symbol=:majority_protein_acs)
    protac_counts = length.(findall.(Ref(";"), coalesce.(peptides_df[!, protein_acs_col], "")))
    pos_counts = length.(findall.(Ref(";"), coalesce.(peptides_df.peptide_poses, "")))
    if protac_counts != pos_counts
        @show peptides_df[protac_counts .!= pos_counts, [:peptide_id, protein_acs_col, :peptide_poses]]
        error("Mismatches in protein and peptide positions")
    end
    poses_df = DelimDataUtils.expand_delim_column(peptides_df, list_col=:peptide_poses, elem_col=:peptide_poses2, key_col=:peptide_id)
    acs_df = DelimDataUtils.expand_delim_column(peptides_df, list_col=protein_acs_col, elem_col=:protein_ac, key_col=:peptide_id)
    @assert poses_df.peptide_id == acs_df.peptide_id
    # one peptide may have multiple matches in protein sequence (comma-separated)
    poses_df.rowix = 1:nrow(poses_df)
    pos_df = DelimDataUtils.expand_delim_column(poses_df, list_col=:peptide_poses2, elem_col=:peptide_pos, key_col=:rowix, delim=",")
    pos_df.peptide_id = poses_df[pos_df.rowix, :peptide_id]
    pos_df.protein_ac = acs_df[pos_df.rowix, :protein_ac]
    pos_df.peptide_pos = parse.(Int, pos_df.peptide_pos)
    return select!(pos_df, Not(:rowix))
end

# find peptide occurrenes in the sequences
function match_peptides(peptides_df::AbstractDataFrame,
                        seqs::Union{Nothing, AbstractDict}=nothing;
                        protein_acs_col::Symbol=:peptide_protein_acs)
    res = DataFrame()
    for r in eachrow(peptides_df)
        pep = replace(r.peptide_seq, r"^_|_$" => "") # remove _ guards
        acs = r[protein_acs_col]
        for ac in split(acs, ';')
            seq = get(seqs, ac, nothing)
            if isnothing(seq)
                @warn "No sequence found for $(ac)"
                continue
            end
            matches = findall(pep, seq)
            isempty(matches) && @warn "No occurences of $(pep) found in $(ac)"
            # one peptide may have multiple matches in protein sequence (comma-separated)
            for m in matches
                push!(res, (peptide_id = r.peptide_id,
                            protein_ac = String(ac),
                            peptide_pos = first(m)))
            end
        end
    end
    return res
end

# generate the dataframe with all potential ptm sites
# for given protein sequences
function ptmsites(proteins_df::AbstractDataFrame,
                  ptmtype2aas::Union{AbstractVector, Nothing}=nothing)
    res = DataFrame()
    for r in eachrow(proteins_df)
        for (pos, aa) in enumerate(r.seq)
            if isnothing(ptmtype2aas)
                ptmr = (ptm_AA_seq = aa, ptm_pos = pos,
                        protein_ac = r.protein_ac)
                push!(res, ptmr)
            else
                for (ptmtype, aas) in ptmtype2aas
                    if occursin(aa, aas)
                        ptmr = (ptm_type = ptmtype, ptm_AAs = aas, ptm_AA_seq = aa, ptm_pos = pos,
                                protein_ac = r.protein_ac)
                        push!(res, ptmr)
                    end
                end
            end
        end
    end
    return res
end

regex2string(r::Regex) = replace(string(r), r"^r\"|\"$" => "")

abstract type AbstractPTMAnnotationFormat end

struct MaxQuantPTMAnnotationFormat <: AbstractPTMAnnotationFormat end
struct SpectronautPTMAnnotationFormat <: AbstractPTMAnnotationFormat end

const PTM_ANNOTATION_FORMATS = Dict(
    :maxquant => MaxQuantPTMAnnotationFormat(),
    :spectronaut => SpectronautPTMAnnotationFormat()
)

PTMAnnotationFormat(format::Symbol) = get!(PTM_ANNOTATION_FORMATS, format) do f
    throw(ArgumentError("Unsupported PTM annotation format $f"))
end

const PTM_REGEX  = r"(?<type>\w+)\s\((?<AA>\w+|Protein [NC]-term)\)" # PTM annotation regex without brackets
const SPECTRONAUT_PTM_REGEX = Regex(string("\\[", regex2string(PTM_REGEX), "\\]"))
const MAXQUANT_PTM_REGEX = Regex(string("\\(", regex2string(PTM_REGEX), "\\)"))

# regex for PTM annotation within modified peptide sequence
annotation_regex(::MaxQuantPTMAnnotationFormat) = MAXQUANT_PTM_REGEX
annotation_regex(::SpectronautPTMAnnotationFormat) = SPECTRONAUT_PTM_REGEX

# regex for annotation of PTM localization probability within modified peptide sequence
locprob_regex(::MaxQuantPTMAnnotationFormat) = r"\((?<locprob>\d+(?:\.\d+)?)\)"
locprob_regex(::SpectronautPTMAnnotationFormat) = r"\[(?<type>\w+)\s\((?<AA>\w+|Protein [NC]-term)\):\s(?<locprob>\d+(?:\.\d+)?)%\]"

# extract inidividual PTMs and their offsets from the modified sequence
function append_ptms!(out::AbstractDataFrame,
                      format::AbstractPTMAnnotationFormat,
                      object_id::Integer,
                      modseq::AbstractString;
                      selptms::Union{Regex, Nothing}=nothing)
    ptm_regex = annotation_regex(format)
    aaseq = replace(replace(modseq, ptm_regex => ""), "_"=>"")
    curoffset = 0 # offset of AA between aaseq (stripped of PTM annots) and seq
    nptms = nselptms = 0
    for m in eachmatch(ptm_regex, modseq)
        is_selptm = isnothing(selptms) || occursin(selptms, m[:type])
        push!(out, (object_id=object_id, peptide_seq=aaseq,
              ptm_type=m[:type],
              ptm_AAs=m[:AA], ptm_AA_seq=modseq[m.offset-1],
              ptm_offset=m.offset-1 - curoffset-2,
              is_selptm = is_selptm))
        curoffset += length(m.match)
        nptms += 1
        is_selptm && (nselptms += 1)
    end
    return out, nptms, nselptms
end

function extract_ptms(df::AbstractDataFrame,
    format::AbstractPTMAnnotationFormat;
    objid_col::Symbol=:pepmod_id,
    modseq_col::Symbol=:pepmod_seq,
    selptms::Union{Regex, Nothing}=nothing
)
    ptms_df = DataFrame()
    obj_ptmstats_df = select(df, [objid_col])
    obj_ptmstats_df[!, :nptms] .= 0
    obj_ptmstats_df[!, :nselptms] .= 0
    for (i, r) in enumerate(eachrow(df))
        _, nptms, nselptms = append_ptms!(ptms_df, format, r[objid_col], r[modseq_col], selptms=selptms)
        obj_ptmstats_df[i, :nptms] = nptms
        obj_ptmstats_df[i, :nselptms] = nselptms
    end
    categorical!(ptms_df, [:ptm_type, :ptm_AAs])
    rename!(ptms_df, :object_id=>objid_col), obj_ptmstats_df
end

extract_ptms(df::AbstractDataFrame; format::Symbol=:maxquant, kwargs...) =
    extract_ptms(df, PTMAnnotationFormat(format); kwargs...)

function elements2sets!(els2sets::Dict{K, VS}, df::AbstractDataFrame, setid_col::Symbol, elem_col::Symbol,
                        group_cols::Union{AbstractVector{Symbol}, Nothing}) where {K <: Tuple, VS <: AbstractVector}
    N = length(K.types)
    groupby_cols = [setid_col]
    isnothing(group_cols) || append!(groupby_cols, group_cols)
    for set_df in groupby(df, groupby_cols)
        key = ntuple(i -> i < N ? set_df[1, group_cols[i]] : Set(set_df[:, elem_col]), N)
        push!(get!(() -> VS(), els2sets, key), set_df[1, setid_col])
    end
    setgroups = sort([kv[1] => sort!(kv[2]) for kv in pairs(els2sets)],
                     by=kv -> (Base.front(kv[1])..., kv[2]))
    return setgroups
end

const PTM_AAs_TERM = ["Protein N-term", "Protein C-term"] # ptm_AAs column value for terminal PTMs

@inline default_locprob_sequence_column(::MaxQuantPTMAnnotationFormat) = :ptm_locprob_seq
@inline default_locprob_sequence_column(::SpectronautPTMAnnotationFormat) = :EG_locprob_seq

@inline default_locprob_ptmtype_column(::MaxQuantPTMAnnotationFormat) = :ptm_type
@inline default_locprob_ptmtype_column(::SpectronautPTMAnnotationFormat) = nothing

@inline locprob_regex(::MaxQuantPTMAnnotationFormat) = r"\((?<locprob>\d+(?:\.\d+)?)\)"
@inline locprob_regex(::SpectronautPTMAnnotationFormat) = r"\[(?<type>\w+)\s\((?<AAs>\w+|Protein [NC]-term)\):\s(?<locprob>\d+(?:\.\d+)?)%\]"

@inline parse_locprob(probstr::AbstractString, ::MaxQuantPTMAnnotationFormat) = parse(Float64, probstr)
@inline parse_locprob(probstr::AbstractString, ::SpectronautPTMAnnotationFormat) = parse(Float64, probstr) / 100

function append_ptm_locprobs!(out::AbstractDataFrame,
    format::AbstractPTMAnnotationFormat,
    ptm_type::Union{AbstractString, Nothing},
    report_ix::Integer,
    object_id::Integer,
    msrun_id::Any,
    locprobseq::AbstractString,
    modseq::Union{AbstractString, Nothing}=nothing,
    ptmsdf_buf::DataFrame=nothing
)
    ptms_df = isnothing(modseq) ? nothing : append_ptms!(isempty(ptmsdf_buf) ? ptmsdf_buf : mapcols!(empty!, ptmsdf_buf), format, object_id, modseq)[1]
    curoffset = 0 # offset of AA between aaseq (stripped of PTM annots) and seq
    ptmix = 0
    locprobstart = !isempty(locprobseq) && locprobseq[1] == '_' ? 1 : 0
    for m in eachmatch(locprob_regex(format), locprobseq)
        locprob = parse_locprob(m[:locprob], format)
        offset = m.offset - locprobstart - curoffset-2
        aa = ptm_AA_seq = locprobseq[m.offset-1]
        curoffset += length(m.match)
        if !isnothing(ptms_df)
            ptm_ix = searchsortedfirst(ptms_df.ptm_offset, offset)
            if (ptm_ix > nrow(ptms_df) || ptms_df.ptm_offset[ptm_ix] != offset)
                continue # skip ptm if it's not in the modified sequence (locprob too low)
            else
                cur_ptm_type = isnothing(ptm_type) ? m[:type] : ptm_type
                #if (aa != ptms_df.ptm_AA_seq[ptm_ix]) || (ptms_df.ptm_type[ptm_ix] != cur_ptm_type)
                #    @show modseq locprobseq ptms_df offset ptm_ix aa ptms_df.ptm_AA_seq[ptm_ix]
                #end
                @assert aa == ptms_df.ptm_AA_seq[ptm_ix] "object_id=$(object_id), pos=$(ptm_ix): Found AA=$aa, expected $(ptms_df.ptm_AA_seq[ptm_ix])"
                # allow PTM type mismatch if pepmodseq one is terminal modification
                @assert (ptms_df.ptm_type[ptm_ix] == cur_ptm_type) ||
                    in(ptms_df.ptm_AAs[ptm_ix], PTM_AAs_TERM) "object_id=$(object_id), pos=$(ptm_ix): Found ptm_type=$(cur_ptm_type), expected $(ptms_df.ptm_type[ptm_ix])"
            end
        end
        push!(out, (report_ix=report_ix, object_id=object_id, msrun_id=msrun_id,
            ptm_type=cur_ptm_type,
            ptm_AA_seq=aa,
            ptm_offset=offset,
            ptm_locprob=locprob))
    end
    return out
end

function extract_ptm_locprobs(df::AbstractDataFrame,
    format::AbstractPTMAnnotationFormat;
    objid_col::Symbol=:pepmodstate_id,
    msrun_col::Union{AbstractVector, Symbol}=:rawfile_ix,
    ptmtype_col::Union{Symbol, Nothing}=default_locprob_ptmtype_column(format),
    locprobseq_col::Symbol=default_locprob_sequence_column(format),
    modseq_col::Union{Symbol, Nothing}=nothing
)
    res = DataFrame()
    ptmsbuf_df = isnothing(modseq_col) ? nothing : DataFrame()
    for (i, r) in enumerate(eachrow(df))
        locprobseq = r[locprobseq_col]
        ismissing(locprobseq) ||
            append_ptm_locprobs!(res, format,
                                 isnothing(ptmtype_col) ? nothing : convert(String, r[ptmtype_col]),
                                 i, r[objid_col], r[msrun_col], locprobseq,
                                 isnothing(modseq_col) ? nothing : r[modseq_col],
                                 ptmsbuf_df)
    end
    categorical!(res, [:ptm_type])
    rename!(res, :object_id=>objid_col)
    if msrun_col isa Symbol
        rename!(res, :msrun_id => msrun_col)
    else
        for (i, col) in enumerate(msrun_col)
            res[!, col] = getindex.(res.msrun_id, i)
        end
        select!(res, Not(:msrun_id))
    end
    return res
end

extract_ptm_locprobs(df::AbstractDataFrame; format::Symbol=:maxquant, kwargs...) =
    extract_ptm_locprobs(df, PTMAnnotationFormat(format); kwargs...)

function flanking_sequence(seq::AbstractString, pos::Integer; flanklen::Integer=15)
    start_flank = pos - flanklen
    res = repeat('_', length(start_flank:0))
    res *= uppercase(seq[max(start_flank, 1):pos-1])
    res *= uppercase(seq[pos])
    end_flank = pos + flanklen
    res *= uppercase(seq[pos+1:min(end_flank, length(seq))])
    res *= repeat('_', length(length(seq)+1:end_flank))
    return res
end

# map the positions of aaobjs from srcseqs to destdeqs of the same group using global pairwise alignment
function map_aapos(aaobjs_df::AbstractDataFrame,
                   srcseqs_df::AbstractDataFrame,
                   destseqs_df::AbstractDataFrame;
                   seqgroup_col::Union{AbstractVector, Symbol}=:protein_ac,
                   srcseqid_col::Union{AbstractVector, Symbol}=seqgroup_col,
                   destseqid_col::Union{AbstractVector, Symbol}=seqgroup_col,
                   obj_prefix::Symbol=:ptm_,
                   objpos_col::Symbol=Symbol(obj_prefix, "pos"),
                   objid_col::Union{AbstractVector, Symbol}=Symbol(obj_prefix, "id"),
                   destmap_prefix::Symbol=:dest_,
                   scoremodel::AbstractScoreModel=AffineGapScoreModel(BLOSUM90, gap_open=-10, gap_extend=-1),
                   selfmatches::Bool=true,
                   verbose::Bool=false
)
    objs_df2 = select(aaobjs_df, unique!([objid_col; srcseqid_col; objpos_col]), copycols=false)
    srcseqs_df2 = select(srcseqs_df, unique!([seqgroup_col; srcseqid_col; :seq]), copycols=false)
    srcseqs_df2.srcseq_ix = 1:nrow(srcseqs_df2)
    srcseqs_df2.srcseq_len = length.(srcseqs_df.seq)
    destseqs_df2 = select(destseqs_df, unique!([seqgroup_col; destseqid_col; :seq]), copycols=false)
    if (destseqid_col == srcseqid_col) && (destseqid_col != seqgroup_col)
        if destseqid_col isa Symbol
            used_destseqid_col = Symbol(destmap_prefix, destseqid_col)
            rename!(destseqs_df2, destseqid_col => used_destseqid_col)
        else
            used_destseqid_col = [in(col, seqgroup_col) ? col : Symbol(destmap_prefix, col) for col in destseqid_col]
            rename!(destseqs_df2, [col => newcol for (col, newcol) in zip(destseqid_col, used_destseqid_col) if col != newcol]...)
        end
    else
        used_destseqid_col = destseqid_col
    end
    destseqs_df2.destseq_ix = 1:nrow(destseqs_df2)
    used_srcseqs_df2 = semijoin(srcseqs_df2, aaobjs_df, on=srcseqid_col)
    src2dest_df = innerjoin(select(used_srcseqs_df2, unique!([:srcseq_ix; seqgroup_col; srcseqid_col; :srcseq_len]), copycols=false),
                            select(destseqs_df2, unique!([:destseq_ix; seqgroup_col; used_destseqid_col]), copycols=false),
                            on=seqgroup_col)
    selfmatches || filter!(r -> r[srcseqid_col] != r[used_destseqid_col], src2dest_df)
    src2dest_df.src2dest_ix = 1:nrow(src2dest_df)
    verbose && @info "Calculating $(nrow(src2dest_df)) pairwise alignments of $(nrow(used_srcseqs_df2)) source $(obj_prefix) seqs to $(nrow(destseqs_df2)) $(destmap_prefix) seqs..."
    src2dest_agns = [pairalign(GlobalAlignment(),
                        LongAminoAcidSeq(srcseqs_df.seq[r.srcseq_ix]), LongAminoAcidSeq(destseqs_df.seq[r.destseq_ix]),
                        scoremodel)
                     for r in eachrow(src2dest_df)]
    obj2dest_df = leftjoin(objs_df2, src2dest_df, on=srcseqid_col)
    destpos_col, destmatch_col, destaa_col = Symbol.(destmap_prefix, obj_prefix, ("pos", "match", "AA"))
    agnseqlen_col, agnscore_col, agnmatchratio_col, agnpos_col = Symbol.("agn_", ("len", "score", "match_ratio", Symbol(obj_prefix, "pos")))
    obj2dest_df[!, destpos_col] = missings(Int, nrow(obj2dest_df))
    obj2dest_df[!, destmatch_col] = missings(Char, nrow(obj2dest_df))
    obj2dest_df[!, destaa_col] = missings(Char, nrow(obj2dest_df))
    obj2dest_df[!, agnseqlen_col] = missings(Int, nrow(obj2dest_df))
    obj2dest_df[!, agnscore_col] = missings(Float64, nrow(obj2dest_df))
    obj2dest_df[!, agnmatchratio_col] = missings(Float64, nrow(obj2dest_df))
    obj2dest_df[!, agnpos_col] = missings(Int, nrow(obj2dest_df))

    for r in eachrow(obj2dest_df)
        (!ismissing(r.src2dest_ix) && (r[objpos_col] <= r.srcseq_len)) || continue
        agn = alignment(src2dest_agns[r.src2dest_ix])
        aln = agn.a.aln
        s2a = seq2agn(agn, r[objpos_col])
        s2r = seq2ref(agn, r[objpos_col])
        r[agnpos_col] = s2a[1]
        r[agnmatchratio_col] = count_matches(agn) / count_aligned(agn)
        r[agnscore_col] = score(src2dest_agns[r.src2dest_ix])
        r[agnseqlen_col] = count_aligned(agn)
        r[destmatch_col] = s2r[2]
        destseq = destseqs_df.seq[r.destseq_ix]
        if 1 <= s2r[1] <= length(destseq)
            r[destpos_col] = s2r[1]
            r[destaa_col] = destseq[s2r[1]]
        end
    end
    return select(obj2dest_df, Not([:srcseq_len, :src2dest_ix]))
end

default_ptm_labeler(r::Any; groupid_col::Symbol=:genename, obj_prefix::Symbol=:ptm_, idextra::Any=nothing) =
    coalesce(r[Symbol(obj_prefix, "type")], "Unknown") * "_" * string(coalesce(r[groupid_col], "Unknown")) *
    (isnothing(idextra) ? "" : string("-", idextra)) * "_" *
    r[Symbol(obj_prefix, "AA_seq")] * string(r[Symbol(obj_prefix, "pos")])

# group aaobjs that map to the same positions on the sequences of the same group
function group_aaobjs(
    aaobjs_df::AbstractDataFrame,
    seqs_df::AbstractDataFrame;
    seqgroup_col::Union{AbstractVector, Symbol}=:protgroup_id,
    seqid_col::Union{AbstractVector, Symbol}=:protein_ac,
    seqrank_col::Union{Symbol, Nothing}=nothing,
    obj_prefix::Symbol=:ptm_,
    objid_col::Union{AbstractVector, Symbol}=Symbol(obj_prefix, "id"),
    objpos_col::Symbol=Symbol(obj_prefix, "pos"),
    objgroupid_col::Symbol=Symbol(obj_prefix, :id),
    objgrouplabel_col::Symbol=Symbol(obj_prefix, :label),
    force_refseqs::Bool=true,
    verbose::Bool=false,
    labeler::Function=(r, idextra) -> default_ptm_labeler(r, groupid_col=isa(seqgroup_col, Symbol) ? seqgroup_col : first(seqgroup_col),
                                                          obj_prefix=obj_prefix, idextra=idextra)
)
    # select the required columns from aaobjs_df
    aaobjs_df2 = select(aaobjs_df, unique!([objid_col; seqid_col; objpos_col]), copycols=false)
    seq_cols = [seqgroup_col; isnothing(seqrank_col) ? seqid_col : seqrank_col; seqid_col;] |> unique!
    # add sequence group and rank columns
    aaobjs_df2 = leftjoin(aaobjs_df2, select(seqs_df, seq_cols, copycols=false), on=seqid_col)
    @assert nrow(aaobjs_df2) == nrow(aaobjs_df)
    sort!(aaobjs_df2, seq_cols)
    aaobjs_df2.aaobj_ix = 1:nrow(aaobjs_df) # simple primary key
    # flag whether the expected aaobj was observed
    objobs_col = Symbol(obj_prefix, "is_observed")
    aaobjs_df2[!, objobs_col] .= true
    # align aaobjs onto the sequences of the same group
    destseqs_df = semijoin(seqs_df, aaobjs_df2, on=seqgroup_col)
    aaobjs_agn_df = map_aapos(aaobjs_df2, seqs_df, destseqs_df,
                        seqgroup_col=seqgroup_col,
                        srcseqid_col=seqid_col,
                        destseqid_col=seqid_col,
                        objid_col=unique!([:aaobj_ix; objid_col]),
                        obj_prefix=obj_prefix,
                        destmap_prefix=:agn_, selfmatches=false, verbose=verbose)
    # remove non-exact matches
    filter!(r -> coalesce(r.agn_ptm_match, ' ') == '=', aaobjs_agn_df)
    agnobj_colmap = [col => Symbol(:agn_, col) for col in [:aaobj_ix; objpos_col; seqid_col] if !in(col, [seqgroup_col;])]
    expected_agnobjs_df = rename(copy(aaobjs_df2, copycols=false), agnobj_colmap...)
    agnjoin_cols = Symbol.(intersect(names(expected_agnobjs_df), names(aaobjs_agn_df)))
    @debug "Joining aaobjs alignment and expected aaobjs by: $(join(agnjoin_cols, ", "))"
    # if after the alignment some aaobj matches another aaobj,
    # the latter would be in expected_agnobjs_df
    objagns_df = leftjoin(aaobjs_agn_df, expected_agnobjs_df, on=agnjoin_cols)
    # remove unobserved alignment predictions or keep them,
    # so that we can use better ranked sequence as the reference
    # e.g. some aaobjs are only observed with non-canonical isoforms,
    # although they could be matched to the canonical one as well.
    # using canonical (although unobserved) might be preferred, because
    # annotation databases report annotation positions for canonical sequences
    if force_refseqs
        # extract unique aligned objects that didn't match observed ones
        unobjs_df = select(filter(r -> ismissing(r.agn_aaobj_ix), objagns_df),
                            [agnjoin_cols; :destseq_ix]) |> unique!
        rename!(unobjs_df, [agncol => col for (col, agncol) in agnobj_colmap if col != :aaobj_ix]...)
        unobjs_df[!, seqrank_col] = destseqs_df[unobjs_df.destseq_ix, seqrank_col]
        select!(unobjs_df, Not(:destseq_ix))
        sort!(unobjs_df, unique!([seqgroup_col; seqid_col; objpos_col]))
        unobjs_df[!, objobs_col] .= false
        # assign simple primary key to unobserved aaobjs
        unobjs_df[!, :aaobj_ix] .= nrow(aaobjs_df2) .+ (1:nrow(unobjs_df))
        # add unobserved objects
        aaobjs_df2 = vcat(aaobjs_df2, unobjs_df) |> unique!
        @assert 1:nrow(aaobjs_df2) == aaobjs_df2.aaobj_ix
        expected_agnobjs_df = rename(copy(aaobjs_df2, copycols=false), agnobj_colmap...)
        objagns_df = leftjoin(aaobjs_agn_df, expected_agnobjs_df, on=agnjoin_cols)
    else
        filter!(r -> !ismissing(r.agn_aaobj_ix), objagns_df)
    end
    # fill in the newly assigned simple primary key in objagns_df
    obj2obj_df = select(objagns_df, [:agn_aaobj_ix, :aaobj_ix]) |> unique!
    obj2group = copy(aaobjs_df2.aaobj_ix)
    group2objs = [Set((i,)) for i in obj2group]
    @inbounds for r in eachrow(obj2obj_df)
        gr = obj2group[r.aaobj_ix]
        agngr = obj2group[r.agn_aaobj_ix]
        if gr != agngr
            destgr, srcgr = minmax(gr, agngr)
            srcgr_objs = group2objs[srcgr]
            for ix in srcgr_objs
                obj2group[ix] = destgr
            end
            union!(group2objs[destgr], srcgr_objs)
            empty!(srcgr_objs)
        end
    end
    nextgroupid = 0
    groupid = fill(0, length(group2objs))
    for (i, gr) in enumerate(group2objs)
        isempty(gr) && continue
        groupid[i] = (nextgroupid += 1)
    end
    aaobjs_df2[!, objgroupid_col] = groupid[obj2group]
    if !isnothing(seqrank_col)
        # assign group id using the position in highest-ranking sequence
        objisref_col = Symbol(obj_prefix, :is_reference) 
        aaobjs_df2[!, objisref_col] .= false
        aaobjs_df2[!, objgrouplabel_col] .= ""
        label2id = Dict{String, Int}()
        for gr in groupby(aaobjs_df2, objgroupid_col)
            refix = findmin(gr[!, seqrank_col])[2]
            gr[refix, objisref_col] = true
            id = gr[1, objgroupid_col]
            # generate unique label
            idextra = coalesce(gr[refix, seqrank_col], 1)
            label = ""
            while true
                label = labeler(gr[refix, :], idextra > 1 ? idextra : nothing)
                previd = get(label2id, label, nothing)
                if isnothing(previd) || (previd == id)
                    break
                end
                #@warn "Duplicate $objgroupid_col=$groupid ($objgroupix_col=$prevgroupix and $objgroupix_col=$groupix)"
                idextra += 1
            end
            @assert !haskey(label2id, label)
            label2id[label] = id
            gr[:, objgrouplabel_col] .= label
        end
        aaobjs_df2[!, objgrouplabel_col] = categorical(aaobjs_df2[!, objgrouplabel_col])
    end
    # keep either the reference or observed
    filter!(r -> r[objobs_col] || r[objisref_col], aaobjs_df2)
    return select!(aaobjs_df2, Not([:aaobj_ix]))
end

function setgroup_frame(df::AbstractDataFrame, setgroup_col::Symbol, setid_col::Symbol, elem_col::Symbol,
                        group_cols::Union{AbstractVector{Symbol}, Nothing})
    KTypes = Type[]
    res_cols = Vector{Any}()
    res_names = Vector{Symbol}()
    if !isnothing(group_cols)
        append!(res_names, group_cols)
        for col in group_cols
            T = eltype(df[col])
            push!(KTypes, T)
            push!(res_cols, Vector{T}())
        end
    end
    E = eltype(df[elem_col])
    S = eltype(df[setid_col])
    push!(res_cols, Int[]); push!(res_names, setgroup_col)
    push!(res_cols, S[]); push!(res_names, setid_col)
    push!(KTypes, Set{E})
    K = Tuple{KTypes...}
    T = eltype(df[setid_col])
    els2sets = elements2sets!(Dict{K, Vector{S}}(), df, setid_col, elem_col, group_cols)
    for (i, (key, setids)) in enumerate(els2sets)
        for setid in setids
            for (j, keypart) in enumerate(Base.front(key))
                push!(res_cols[j], keypart)
            end
            push!(res_cols[end-1], i)
            push!(res_cols[end], setid)
        end
    end
    return DataFrame(res_cols, res_names)
end

# group PTMs sharing the same data (same offsets of the same peptides/pepmodstates)
function ptmgroups(ptmn2pms_df::AbstractDataFrame; verbose::Bool=false)
    ptmn2pms_df = copy(ptmn2pms_df, copycols=false)

    pepXoffs = [(r.peptide_id, r.ptm_offset) for r in eachrow(ptmn2pms_df)]
    pepXoffs_unique = sort!(unique(pepXoffs))
    ptmn2pms_df.peptideXoffset = searchsortedfirst.(Ref(pepXoffs_unique), pepXoffs)

    pmsXoffs = [(r.pepmodstate_id, r.ptm_offset) for r in eachrow(ptmn2pms_df)]
    pmsXoffs_unique = sort!(unique(pmsXoffs))
    ptmn2pms_df.pepmodstateXoffset = searchsortedfirst.(Ref(pmsXoffs_unique), pmsXoffs)

    ptm2ptmgroup_df = setgroup_frame(ptmn2pms_df, :ptmgroup_id, :ptm_id, :peptideXoffset, [:ptm_type])
    ptmn2ptmngroup_df = setgroup_frame(ptmn2pms_df, :ptmngroup_id, :ptmn_id, :pepmodstateXoffset, [:ptm_type, :nselptms])
    ptm2ptmn_df = unique!(select(ptmn2pms_df, [:ptm_id, :ptmn_id]))
    ptmgroup2ptmngroup_df = DataFrame()
    # TODO: use protgroup_asembly()?
    while true
        ptmgroup2ptmngroup_df = unique!(select!(innerjoin(innerjoin(ptmn2ptmngroup_df, ptm2ptmn_df, on=:ptmn_id), ptm2ptmgroup_df, on=[:ptm_type, :ptm_id]),
                                        [:ptm_type, :nselptms, :ptmgroup_id, :ptmngroup_id]))
        # we are not done yet: there could be ptmgroups that are strict subsets of the other ptmgroups (i.e. contain subsets of the ptmngroups)
        # that's because at this modent some ptms could belong to multiple ptmngroups (because of peptides with modifications)
        # 1. identify ptmgroups sharing ptmns
        pgn2pgs = Dict{Int, Set{Int}}()
        for r in eachrow(ptmgroup2ptmngroup_df)
            push!(get!(() -> Set{Int}(), pgn2pgs, r.ptmngroup_id), r.ptmgroup_id)
        end
        filter!(kv -> length(kv[2]) > 1, pgn2pgs)
        if !isempty(pgn2pgs)
            verbose && @info "$(length(pgn2pgs)) PTMn groups are shared by PTM groups"
            # 2. merge them
            pg2merged = Dict{Int, Int}()
            for pgs in values(pgn2pgs)
                newpg = minimum(pgs)
                for pg in pgs
                    pg2merged[pg] = newpg
                end
            end
            ptm2ptmgroup_df.ptmgroup_id = [get(pg2merged, r.ptmgroup_id, r.ptmgroup_id) for r in eachrow(ptm2ptmgroup_df)]
            # 3. remove indexing gaps
            ptmgroup_ids = sort!(unique(ptm2ptmgroup_df.ptmgroup_id))
            ptm2ptmgroup_df.ptmgroup_id = searchsortedfirst.(Ref(ptmgroup_ids), ptm2ptmgroup_df.ptmgroup_id)
            # 4. repeat: update ptmgroup2ptmngroup and check for sharing again
        else
            break # we are done
        end
    end

    return ptm2ptmgroup_df, ptmn2ptmngroup_df, ptmgroup2ptmngroup_df
end

end
