module ProSite

using DataFrames, CategoricalArrays, FastaIO

function read_prosite_alignments(folder; verbose=false)
    agnfiles = readdir(folder)
    verbose && @info("Found $(length(agnfiles)) in $folder, processing...")
    prot_names = Vector{String}()
    prot_acs = Vector{String}()
    seq_ranges = Vector{UnitRange}()
    motif_names = Vector{String}()
    motif_ids = Vector{String}()
    motif_ranges = Vector{Union{UnitRange, Missing}}()
    for agnfile in agnfiles
        endswith(agnfile, ".msa") || continue
        verbose && @info("Parsing $(agnfile)")
        fr = FastaReader(joinpath(folder, agnfile))
        while !eof(fr)
            header, seq = readentry(fr)
            hmatch = match(r"^(.+)\|(.+)\/(\d+)-(\d+)\: (.+)\|(.+)$", header)
            push!(prot_names, hmatch[1])
            push!(prot_acs, hmatch[2])
            push!(seq_ranges, parse(Int, hmatch[3]):parse(Int, hmatch[4]))
            push!(motif_names, hmatch[5])
            mmatch = match(r"^(.+)\/(\d+)\.(\d+)$", hmatch[6])
            if mmatch !== nothing
                push!(motif_ids, mmatch[1])
                push!(motif_ranges, parse(Int, mmatch[2]):parse(Int, mmatch[3]))
            else
                push!(motif_ids, hmatch[6])
                push!(motif_ranges, missing)
            end
        end
    end
    return DataFrame(
        protein_ac = prot_acs,
        protein_name = prot_names,
        seq_range = seq_ranges,
        motif_id = categorical(motif_ids),
        motif_name = categorical(motif_names),
        motif_range = motif_ranges,
    )
end

end
