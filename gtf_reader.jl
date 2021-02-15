module GTF

using CSV, DataFrames

"""
    read(filename) -> DataFrame

Reads Ensembl GTF files.

### See Also
[GTF format description](https://www.ensembl.org/info/website/upload/gff.html).
"""
read(filename::AbstractString) =
    CSV.read(filename, delim='\t', limit=10000000,# missingstrings = ".",
             comment = "#!", #allowmissing=:auto,
             header=[:seqname, :src_db, :feature, :start, :end, :score, :strand, :frame, :attribute])

"""
    parse_attributes(gtf_df::AbstractDataFrame) -> DataFrame

Parses the attribute string of each GTF entry and returns the results as a data frame.

### See Also
[GTF format description](https://www.ensembl.org/info/website/upload/gff.html).
"""
function parse_attributes(attrs::AbstractVector{String})
    taggedvals = Dict{Symbol, Vector{Tuple{Int, String}}}()
    for (row_ix, attr_line) in enumerate(attrs)
        attr_line !== missing || continue
        chunks = split(attr_line, ';', keepempty=false)
        for chunk in chunks
            cm = match(r"^\s*(?<tag>.+)\s\"(?<val>.+)\"\s*$", chunk)
            if cm === nothing
                @warn "Tag value parsing failed: $chunk"
            else
                vals = get!(() -> sizehint!(valtype(taggedvals)(), length(attrs)), taggedvals, Symbol(cm[:tag]::SubString{String}))
                push!(vals, (row_ix, string(cm[:val]::SubString{String})))
            end
        end
    end
    df = DataFrame(Vector{Union{String, Missing}}[begin
        x = Vector{Union{String, Missing}}(missing, length(attrs))
        setindex!(x, last.(vals), first.(vals))
        x end for (tag, vals) in taggedvals], collect(keys(taggedvals)))
    df.row_ix = eachindex(attrs)
    return df
end

end
