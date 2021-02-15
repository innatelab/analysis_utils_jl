module SigDB

using FileIO, RData

# SigDB mapping from http://bioinf.wehi.edu.au/software/MSigDB/

function r2setcoll(::Type{I}, rdict) where I
    coll = sizehint!(Dict{String, Set{I}}(), length(rdict))
    for set_id in keys(rdict)
        coll[set_id] = Set(I[parse(I, e_id) for e_id in rdict[set_id]])
    end
    return coll
end

# FIXME add more collections
const SigDBs = Dict("c2" => "C2", "c6" => "C6", "c7" => "C7", "H" => "Hallmark")
const SigDBOrganisms = Dict("Mm" => "mouse", "Hs" => "human")

"""
Read SigDB data from `sigdbpath` into existing `colls` collection.
"""
function readto!(colls, sigdbpath; organism = "Hs", version = "v5p2")
    I = eltype(valtype(valtype(colls))) # type of gene/protein IDs
    for sigdb in keys(SigDBs)
        coll_file = joinpath(sigdbpath, SigDBOrganisms[organism]*"_"*sigdb*"_"*version*".rdata")
        @info "Loading $sigdb from $coll_file"
        if !isfile(coll_file)
            @warn "SigDB collection $sigdb not found at $coll_file, skipping"
            continue
        end
        coll = r2setcoll(I, load(coll_file, convert=true)[organism*"."*sigdb])
        colls[Symbol(organism * "_SigDB_", SigDBs[sigdb])] = coll
    end
    return colls
end

end
