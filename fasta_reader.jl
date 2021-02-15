module Fasta

using DataFrames, CategoricalArrays, FastaIO

#=
Regular expression for extracting fields from the UniProt Fasta header
=#
const UNIPROT_REGEX = r"^(?<src>sp|tr)\|(?<ac>[^|]+)\|(?<code>\S+)\s(?<name>.+)\sOS=(?<organism>.+?(?=\s(?:GN|PE|SV)\=|$))?(?:\sGN=(?<genename>[^=]+))?(?:\sPE=(?<existence>\S+))?(?:\sSV=(?<seqver>\S+))?$";

function read_uniprot(filename; verbose=false)
    fr = FastaReader(filename)

    headers = Vector{String}()
    seqs = Vector{String}()
    srcs = Vector{Union{String, Missing}}()
    acs = Vector{Union{String, Missing}}()
    codes = Vector{Union{String, Missing}}()
    names = Vector{Union{String, Missing}}()
    organisms = Vector{Union{String, Missing}}()
    genenames = Vector{Union{String, Missing}}()
    existences = Vector{Union{Int, Missing}}()
    seqvers = Vector{Union{Int, Missing}}()
    while !eof(fr)
        header, seq = readentry(fr)
        push!(headers, header)
        push!(seqs, seq)
        hm = match(UNIPROT_REGEX, header)
        if hm === nothing
            @warn "Fasta header parsing failed: $header"
            push!(srcs, missing)
            push!(acs, missing)
            push!(codes, missing)
            push!(names, missing)
            push!(organisms, missing)
            push!(genenames, missing)
            push!(existences, missing)
            push!(seqvers, missing)
        else
            push!(srcs, hm[:src])
            push!(acs, hm[:ac])
            push!(codes, hm[:code])
            push!(names, hm[:name])
            push!(organisms, hm[:organism])
            push!(genenames, hm[:genename] !== nothing ? hm[:genename] : missing)
            push!(existences, hm[:existence] !== nothing ? parse(Int, hm[:existence]) : missing)
            push!(seqvers, hm[:seqver] !== nothing ? parse(Int, hm[:seqver]) : missing)
        end
    end
    return DataFrame(
        src_db = categorical(srcs),
        protein_ac = acs,
        protein_code = codes,
        protein_name = names,
        organism = organisms,
        genename = genenames,
        protein_existence = existences,
        seq_version = seqvers,
        seq = seqs,
        fasta_header = headers,
    )
end

function strip_uniprot_isoform(ac::AbstractString)
    acmatch = match(r"^(\w+)-(?:PRO_)?\d+$", ac)
    acmatch === nothing ? ac : convert(String, acmatch[1])
end

# FIXME based on GRCm38 release 79, later releases extend the format
const ENSEMBL_AC_PART = "(?<ac>ENS\\w+\\d+)"
const ENSEMBL_ANNOTS_PART = "\\spep:(?<pep>.+)\\s(?<loctype>chromosome|supercontig|scaffold):(?<location>.+)\\sgene:(?<geneac>.+)\\stranscript:(?<transcript>.+)\\sgene_biotype:(?<gene_biotype>.+)\\stranscript_biotype:(?<transcript_biotype>.+)"
const ENSEMBL_SAAVs_PART = "(?:\\s+(?<SAAVs>\\S+))?"
const ENSEMBL_AC_REGEX = Regex(ENSEMBL_AC_PART)
const ENSEMBL_ANNOTS_REGEX = Regex(ENSEMBL_AC_PART * ENSEMBL_ANNOTS_PART)
const ENSEMBL_SAAVs_REGEX = Regex(ENSEMBL_AC_PART * ENSEMBL_SAAVs_PART)

function read_ensembl(filename; 
                      has_annotations::Bool=true,
                      has_SAAVs::Bool=false,
                      verbose::Bool=false)
    fr = FastaReader(filename)

    headers = Vector{String}()
    seqs = Vector{String}()
    acs = Vector{Union{String, Missing}}()
    if has_annotations
        @assert !has_SAAVs
        header_regex = ENSEMBL_ANNOTS_REGEX
        peps = Vector{Union{String, Missing}}()
        loctypes = Vector{Union{String, Missing}}()
        locations = Vector{Union{String, Missing}}()
        geneacs = Vector{Union{String, Missing}}()
        transcripts = Vector{Union{String, Missing}}()
        gene_biotypes = Vector{Union{String, Missing}}()
        transcript_biotypes = Vector{Union{String, Missing}}()
    elseif has_SAAVs
        header_regex = ENSEMBL_SAAVs_REGEX
        SAAVs = Vector{Union{String, Missing}}()
    else
        header_regex = ENSEMBL_AC_REGEX
    end
    while !eof(fr)
        header, seq = readentry(fr)
        push!(headers, header)
        push!(seqs, seq)
        hm = match(header_regex, header)
        if isnothing(hm)
            @warn "Fasta header parsing failed: $header"
            push!(acs, missing)
            if has_annotations
                push!(peps, missing)
                push!(loctypes, missing)
                push!(locations, missing)
                push!(geneacs, missing)
                push!(transcripts, missing)
                push!(gene_biotypes, missing)
                push!(transcript_biotypes, missing)
            end
            if has_SAAVs
                push!(SAAVs, missing)
            end
        else
            push!(acs, hm[:ac])
            if has_annotations
                push!(peps, hm[:pep])
                push!(loctypes, hm[:loctype])
                push!(locations, hm[:location])
                push!(geneacs, hm[:geneac])
                push!(transcripts, hm[:transcript])
                push!(gene_biotypes, hm[:gene_biotype])
                push!(transcript_biotypes, hm[:transcript_biotype])
            end
            if has_SAAVs
                push!(SAAVs, !isnothing(hm[:SAAVs]) ? hm[:SAAVs] : missing)
            end
        end
    end
    res = DataFrame(fasta_header = headers,
                    protein_ac = acs,
                    seq = seqs)
    if has_annotations
        res.status = categorical(peps)
        res.loctype = categorical(loctypes)
        res.location = locations
        res.gene_ac = geneacs
        res.transcript_ac = transcripts
        res.gene_biotype = categorical(gene_biotypes)
        res.transcript_biotypes = categorical(transcript_biotypes)
    end
    if has_SAAVs
        res.SAAVs = SAAVs
    end
    return res
end

const CONTAMINANT_REGEX = r"(?<ac>CON__\w+(?:-\d+)?)\s?(?<src>SWISS-PROT|TREMBL|REFSEQ|ENSEMBL|H-INV):(?<ac2>\w+(?:-\d+)?)(?:\||\s|;)(?:\((?<organism>.+)\)\s)?(?<info>.+)"
const CONTAMINANT2_REGEX = r"CON__(?<src>SWISS-PROT|TREMBL|REFSEQ|ENSEMBL|H-INV):(?<ac>\w+(?:-\d+)?)(?:\||\s|;)(?:\((?<organism>.+)\)\s)?(?<info>.+)"

function read_contaminants(filename; verbose=false)
    fr = FastaReader(filename)

    headers = Vector{String}()
    seqs = Vector{String}()
    srcs = Vector{Union{String, Missing}}()
    acs = Vector{Union{String, Missing}}()
    protein_names = Vector{Union{String, Missing}}()
    organisms = Vector{Union{String, Missing}}()
    while !eof(fr)
        header, seq = readentry(fr)
        push!(headers, header)
        push!(seqs, seq)
        hm = match(CONTAMINANT_REGEX, header)
        if hm === nothing
            hm = match(CONTAMINANT2_REGEX, header)
        end
        if hm === nothing
            @warn "Fasta header parsing failed: $header"
            push!(srcs, missing)
            push!(acs, missing)
            push!(protein_names, missing)
            push!(organisms, missing)
        else
            push!(srcs, hm[:src])
            push!(acs, hm[:ac])
            push!(protein_names, hm[:info])
            push!(organisms, isnothing(hm[:organism]) ? missing : hm[:organism])
        end
    end
    return DataFrame(
        src_db = categorical(srcs),
        protein_ac = acs,
        protein_name = protein_names,
        organism = organisms,
        seq = seqs,
        fasta_header = headers,
    )
end

const PHOSPHOSITEPLUS_REGEX = r"GN:(?<genename>[^|]+)\|(?<name>[^|]+)\|(?<organism>[^|]+)\|(?<ac>[^|]*)$"

function read_phosphositeplus(filename)
    fr = FastaReader(filename)

    headers = Vector{String}()
    seqs = Vector{String}()
    genenames = Vector{Union{String, Missing}}()
    acs = Vector{Union{String, Missing}}()
    protein_names = Vector{Union{String, Missing}}()
    organisms = Vector{Union{String, Missing}}()

    minus2missing(x) = x == "-" ? missing : x

    while !eof(fr)
        header, seq = readentry(fr)
        push!(headers, header)
        push!(seqs, seq)
        hm = match(PHOSPHOSITEPLUS_REGEX, header)
        if hm === nothing
            @warn "Fasta header parsing failed: $header"
            push!(genenames, missing)
            push!(acs, missing)
            push!(protein_names, missing)
            push!(organisms, missing)
        else
            push!(genenames, minus2missing(hm[:genename]))
            push!(acs, minus2missing(hm[:ac]))
            push!(protein_names, minus2missing(hm[:name]))
            push!(organisms, minus2missing(hm[:organism]))
        end
    end
    return DataFrame(
        protein_ac = acs,
        protein_name = protein_names,
        genename = genenames,
        organism = organisms,
        seq = seqs,
        fasta_header = headers,
    )
end

end
