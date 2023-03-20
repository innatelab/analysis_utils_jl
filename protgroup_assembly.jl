module ProtgroupAssembly

using DataFrames, Graphs
using Printf: @sprintf

# the relationship between the two sets of T elements
struct SetXSet{T}
    common::Set{T}  # elements common to both sets (left & right)
    left::Set{T}    # left set-only elements
    right::Set{T}   # right set-only elements

    SetXSet{T}() where T =
        new{T}(Set{T}(), Set{T}(), Set{T}())
end

Base.eltype(::Type{SetXSet{T}}) where T = T
Base.eltype(x::SetXSet) = eltype(typeof(x))

struct Protgroup{T,P}
    major_prots::Set{T}
    prots::Set{T}

    spec_peps::Set{P}
    peps::Set{P}
end

# single linkage clustering of proteins sharing peptides
# one protein in the cluster is representative (id2sg is positive)
function peptide_sharing_groups(idXid::Dict{Tuple{Int, Int}, SetXSet{P}}, nprots::Integer) where P
    maxsg = 0
    id2sg = fill(0, nprots) # sg < 0 is for non-representative prot id of the class
    sg2ids = Dict{Int, Set{Int}}()
    for ((id1, id2), pep_info) in idXid
        @assert (id1 > 0) && (id2 > 0) && (id1 <= nprots) && (id2 <= nprots)
        if isempty(pep_info.left) && isempty(pep_info.right)
            sg1 = id2sg[id1]
            sg2 = id2sg[id2]
            if sg1 == 0 # id1 is not assigned to the sharing group
                if sg2 == 0 # create the new group and make id1 its representative
                    sg = (maxsg += 1)
                    sg2ids[sg] = Set{Int}((id1, id2))
                    id2sg[id1] = id1 <= id2 ? sg : -sg
                    id2sg[id2] = id1 < id2 ? -sg : sg
                else
                    # assign id1 as non-representative of the group of id2
                    sg = abs(sg2)
                    id2sg[id1] = -sg
                    push!(sg2ids[sg], id1)
                end
            else # id1 is assigned to a group
                sg = abs(sg1)
                if sg2 == 0 # add id2 to the group of id1 as non-representative
                    id2sg[id2] = -sg
                    push!(sg2ids[sg], id2)
                elseif sg != abs(sg2) # merge sg1 and sg2
                    sg, sgold = minmax(sg, abs(sg2))
                    union!(sg2ids[sg], sg2ids[sgold])
                    @inbounds for id in sg2ids[sgold]
                        id2sg[id] = -sg
                    end
                    delete!(sg2ids, sgold)
                end
            end
        end
    end
    # recode sharing group indices to be 1:n
    sgix = sort!(collect(keys(sg2ids)))
    sg_ids = [sg2ids[sg] for sg in sgix]
    @inbounds for (sg, ids) in enumerate(sg_ids)
        for id in ids
            id2sg[id] = id2sg[id] > 0 ? sg : -sg
        end
    end
    return id2sg, sg_ids
end

function assemble_protgroups(pep2prots::Dict{P, Tuple{Set{T}, I}};
                             minspec_peptides::Integer=1,
                             nspec_peptides::Integer=minspec_peptides,
                             rspec_peptides::Real=0.0,
                             verbose::Bool=false) where {T, P, I<:Integer}
    # enumerate proteins
    verbose && @info("Enumerating proteins...")
    protset = Set{T}()
    for (prots, rank) in values(pep2prots)
        (rank > 0) && union!(protset, prots)
    end
    allprots = sort!(collect(protset))
    prot2id = Dict(prot => i for (i, prot) in enumerate(allprots))
    verbose && @info("$(length(prot2id)) protein(s) referenced by $(sum(x -> x[2] > 0, values(pep2prots))) decisive peptide(s) of $(length(pep2prots))")
    # build protein -> peptides map
    verbose && @info("Building protein → peptides maps...")
    # map all peptides
    id2allpeps = Dict{Int, Set{P}}()
    pep_ranks = Vector{Int}()
    for (pep, (prots, pep_rank)) in pep2prots
        if pep_rank > 0
            pos = searchsortedfirst(pep_ranks, pep_rank)
            if pos > length(pep_ranks)
                push!(pep_ranks, pep_rank)
            elseif pep_ranks[pos] != pep_rank
                insert!(pep_ranks, pos, pep_rank)
            end
        end
        for prot in prots
            id = get(prot2id, prot, 0)
            (id == 0) && continue # not referenced by used peptides
            idallpeps = get!(() -> Set{P}(), id2allpeps, id)
            push!(idallpeps, pep)
        end
    end

    # map the peptides used to define the groups
    id2peps_ranked = Dict{Int, Tuple{Int, Set{P}}}()
    for cur_rank in pep_ranks
        for (pep, (prots, pep_rank)) in pep2prots
            if pep_rank == cur_rank
                for prot in prots
                    id = get(prot2id, prot, 0)
                    (id == 0) && continue # not referenced by used peptides
                    prot_pep_rank, idpeps = get!(() -> (cur_rank, Set{P}()), id2peps_ranked, id)
                    # use the mapping only if it's the same as protein rank
                    if prot_pep_rank >= pep_rank
                        push!(idpeps, pep)
                    end
                end
            end
        end
        @info "$(length(id2peps_ranked)) proteins referenced by $cur_rank level peptides"
    end
    # strip ranks
    id2peps = Dict(id => peps for (id, (rank, peps)) in pairs(id2peps_ranked))

    # build all pairs of proteins that share used peptide (including self-self)
    verbose && @info("Collecting protein pairs that share decisive peptide(s)...")
    idXid = Dict{Tuple{Int, Int}, SetXSet{P}}()
    for (pep, (prots, pep_rank)) in pep2prots
        for prot1 in prots
            id1 = get(prot2id, prot1, 0)
            (id1 == 0) && continue
            for prot2 in prots
                id2 = get(prot2id, prot2, 0)
                ((id2 == 0) || (id1 > id2)) && continue
                id_pair = get!(() -> SetXSet{P}(), idXid, (id1, id2))
                push!(id_pair.common, pep)
            end
        end
    end
    verbose && @info("$(length(idXid)) protein pair(s) sharing decisive peptide(s) collected")
    # update protein pairs with all peptides that belong to the left or the right one
    for ((id1, id2), pep_info) in idXid
        (id1 == id2) && continue
        empty!(pep_info.left)
        for pep in id2peps[id1]
            in(pep, pep_info.common) || push!(pep_info.left, pep)
        end
        empty!(pep_info.right)
        for pep in id2peps[id2]
            in(pep, pep_info.common) || push!(pep_info.right, pep)
        end
    end
    # build equality classes (maj.protein.acs) = proteins that have identical peptides
    verbose && @info("Detecting protein(s) with identical decisive peptides...")
    id2sg, sg2ids = peptide_sharing_groups(idXid, length(allprots))

    verbose && @info("$(length(sg2ids)) distinct peptide-sharing group(s) detected...")
    # set sub/super equivalence classes
    # note: no need to include subclasses of subclasses into superclass as
    #       idXid already transitively closed
    verbose && @info("Initializing set-theoretical relationships between peptide-sharing groups...")
    sg2subs = fill!(Vector{Union{Set{Int}, Nothing}}(undef, length(sg2ids)), nothing)
    sg2supers = fill!(Vector{Union{Set{Int}, Nothing}}(undef, length(sg2ids)), nothing)
    sgXsg = Vector{Tuple{Int, Int}}()
    for ((id1, id2), pep_info) in idXid
        sg1 = id2sg[id1]
        sg2 = id2sg[id2]
        ((sg1 > 0) && (sg2 > 0)) || continue # skip non-representative edges
        n1 = length(pep_info.left)
        n2 = length(pep_info.right)
        (n1 == 0) && (n2 == 0) && continue # skip equivalence relation
        if n1 == 0 # sg1 subset of sg2
            subs = isnothing(sg2subs[sg2]) ? (sg2subs[sg2] = Set{Int}()) : sg2subs[sg2]
            push!(subs, sg1)
            sups = isnothing(sg2supers[sg1]) ? (sg2supers[sg1] = Set{Int}()) : sg2supers[sg1]
            push!(sups, sg2)
        elseif n2 == 0 # sg2 subset of sg1
            subs = isnothing(sg2subs[sg1]) ? (sg2subs[sg1] = Set{Int}()) : sg2subs[sg1]
            push!(subs, sg2)
            sups = isnothing(sg2supers[sg2]) ? (sg2supers[sg2] = Set{Int}()) : sg2supers[sg2]
            push!(sups, sg1)
        elseif (rspec_peptides > 0) || (nspec_peptides > 0) # prepare for merging similar equivalence classes
            nb = length(pep_info.common)
            if ((n1 < rspec_peptides * (n1 + nb)) || (n2 < rspec_peptides * (n2 + nb))) &&
               ((n1 < nspec_peptides) || (n2 < nspec_peptides))
                push!(sgXsg, minmax(sg1, sg2))
            end
        end
    end

    # build protein groups: sharing groups not contained in the others
    if !isempty(sgXsg)
        verbose && @info("Merging peptide-sharing groups that share more than " *
                         @sprintf("%.1f%%", (1-rspec_peptides)*100) * " decisive peptide(s) and " *
                         "have less than $nspec_peptides specific decisive peptide(s)")
        # build a graph of sharing groups to merge
        sg_graph = Graph(length(sg2supers));
        for (sg1, sg2) in sgXsg
            if isnothing(sg2supers[sg1]) && isnothing(sg2supers[sg2])
                # add an edge between "proper" sharing-groups
                add_edge!(sg_graph, sg1, sg2)
            end
        end
        sg_components = connected_components(sg_graph)
        nmerges = 0
        nmerged = 0
        for mergedsgs in sg_components
            (length(mergedsgs) == 1) && continue # skip trivial components
            nmerges += 1
            nmerged += length(mergedsgs)
            # add new sharing group
            nextsg = length(sg2subs) + 1
            push!(sg2subs, Set(mergedsgs))
            push!(sg2supers, nothing)
            for subsg in mergedsgs
                sups = isnothing(sg2supers[subsg]) ? (sg2supers[subsg] = Set{Int}()) : sg2supers[subsg]
                push!(sups, nextsg)
            end
        end
        verbose && @info("$(nmerged) sharing group(s) were merged into $(nmerges) new sharing groups")
    end

    # finally assemble protein groups
    verbose && @info("Assembling protein groups...")
    pep2npg = Dict{P, Int}()
    id2npg = Dict{Int, Int}()
    pgs = Vector{Protgroup{T,P}}()
    sg2pg = Dict{Int,Int}()

    function append_sg!(ids, peps, sg)
        if sg <= length(sg2ids)
            sgids = sg2ids[sg]
            union!(ids, sgids)
            for id in sgids
                union!(peps, id2allpeps[id])
            end
        else # merged sg
            @assert !isnothing(sg2subs[sg])
        end
        if !isnothing(sg2subs[sg])
            for subsg in sg2subs[sg]
                append_sg!(ids, peps, subsg)
            end
        end
        return ids
    end

    for (sg, supersgs) in enumerate(sg2supers)
        isnothing(supersgs) || continue # skip subsgs
        peps = Set{P}()
        ids = Set{Int}()
        append_sg!(ids, peps, sg)
        for pep in peps
            pep2npg[pep] = get(pep2npg, pep, 0) + 1
        end
        for id in ids
            id2npg[id] = get(id2npg, id, 0) + 1
        end
        @inbounds prots = Set(allprots[id] for id in ids)
        push!(pgs, Protgroup{T,P}(copy(prots), prots, copy(peps), peps))
        sg2pg[sg] = length(pgs)
    end
    verbose && @info("$(length(pgs)) protein groups assembled, referring to $(length(pep2npg)) peptide(s) of $(length(id2npg)) proteins")
    verbose && @info("Correct group-specific peptides and proteins")
    for pg in pgs
        filter!(pep -> pep2npg[pep] == 1, pg.spec_peps) # leave only ones that belong to one group
        filter!(prot -> id2npg[prot2id[prot]] == 1, pg.major_prots) # leave only proteins that belong to one group
    end
    nnospec = sum(pg -> length(pg.spec_peps) < minspec_peptides, pgs)
    if nnospec > 0
        verbose && @info("Removing $(nnospec) protein groups with less than $(minspec_peptides) specific peptides")
        filter!(pg -> length(pg.spec_peps) >= minspec_peptides, pgs)
    end
    return pgs
end

function dataframe(protgroups::AbstractVector{<:Protgroup};
            protein_ranks::Union{Dict, Nothing}=nothing,
            protgroup_col::Symbol = :protregroup_id,
            protein_col::Symbol = :protgroup_id,
            subobject::Symbol = :peptide,
            proteins_info::Union{AbstractDataFrame, Nothing}=nothing,
            extra_cols::Union{AbstractVector, Nothing}=[:organism, :gene_name, :protein_name, :is_contaminant])
    protsortby = isnothing(protein_ranks) ? identity :
        ac -> (get(protein_ranks, ac, 100), ac)

    res = DataFrame(
        [eachindex(protgroups) .- 1,
         [join(sort!(collect(prg.major_prots), by=protsortby), ';') for prg in protgroups],
         [join(sort!(collect(prg.prots), by=protsortby), ';') for prg in protgroups],
         [join(sort!(collect(prg.spec_peps)), ';') for prg in protgroups],
         [join(sort!(collect(prg.peps)), ';') for prg in protgroups]],
        [protgroup_col,
         Symbol("majority_", protein_col, "s"),
         Symbol(protein_col, "s"),
         Symbol("spec_", subobject, "_ids"), Symbol(subobject, "_ids")]
    )
    if !isnothing(proteins_info) && !isnothing(extra_cols)
        ac2row = Dict(ac => i for (i, ac) in enumerate(proteins_info[!, protein_col]))
        prg2rows = [begin
            rows = Vector{Int}()
            for ac in prg.major_prots
                ix = get(ac2row, ac, 0)
                (ix > 0) && push!(rows, ix)
            end
            rows
        end for prg in protgroups]
        for col in extra_cols
            # FIXME configurable aggregating of columns
            T = nonmissingtype(eltype(proteins_info[!, col]))
            if T <: AbstractString
                col_agg = [begin
                    items = sort!(unique(skipmissing(proteins_info[rows, col])))
                    isempty(items) ? missing : join(items, ';')
                end for rows in prg2rows]
            elseif T === Bool # FIXME assuming any is a right aggregation
                col_agg = [any(r -> proteins_info[r, col], rows) for rows in prg2rows]
            else
                @warn "Don't know how to aggregate column $col of type $T"
                continue
            end
            res[!, col] = col_agg
        end
    end
    return res
end

end
