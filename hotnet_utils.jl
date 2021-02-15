module HotnetUtils

using LinearAlgebra, LightGraphs, DataFrames, CategoricalArrays, CSV, HierarchicalHotNet, StatsBase

const HHN = HierarchicalHotNet
const FrameUtils = Main.FrameUtils
const GraphML = Main.GraphML
const FA = Main.ForceAtlas3

function import_reactomefi(reactomefi_path::AbstractString;
                           edgetypes_path::AbstractString=joinpath(dirname(reactomefi_path), "reactome_fi_directions.txt"),
                           verbose::Bool=false)
    # read ReactomeFI network data frame
    verbose && @info("Reading ReactomeFI network from $reactomefi_path")
    fi_df = CSV.read(reactomefi_path, DataFrame, header=true, delim='\t')
    rename!(x -> Symbol(lowercase(string(x))), fi_df)
    fi_df.direction = categorical(fi_df.direction)
    # generate the list of all genes (source & target) and convert to categorical
    fi_genes = sort!(convert(Vector{String}, unique!(vcat(fi_df.gene1, fi_df.gene2))))
    verbose && @info("$(nrow(fi_df)) interaction(s) of $(length(fi_genes)) gene(s) read")
    fi_gene2index = Dict(gene => i for (i, gene) in enumerate(fi_genes))
    fi_df.gene1 = levels!(categorical(fi_df.gene1), fi_genes)
    fi_df.gene2 = FrameUtils.matchcategorical(fi_df.gene2, fi_df.gene1)

    # read edge type info
    fi_edgetypes_df = CSV.read(edgetypes_path, DataFrame, header=true, delim='\t')
    fi_edgetypes_df.direction = FrameUtils.matchcategorical(fi_edgetypes_df.direction, fi_df.direction)
    fi_edgetypes_df[!, :src_to_dest] .= fi_edgetypes_df.src_to_dest .== 1
    fi_edgetypes_df[!, :dest_to_src] .= fi_edgetypes_df.dest_to_src .== 1
    categorical!(fi_edgetypes_df, :source)
    fi_edgetypes_df.target = FrameUtils.matchcategorical(fi_edgetypes_df.target, fi_edgetypes_df.source)
    fi_df = leftjoin(fi_df, fi_edgetypes_df, on=:direction)

    verbose && @info("Converting to directed graph...")
    fi_diedges_fwd_df = select(fi_df, [:gene1, :gene2, :score, :src_to_dest, :direction, :target])
    rename!(fi_diedges_fwd_df, :src_to_dest => :is_valid, :target => :diedge_type)
    fi_diedges_fwd_df[!, :is_reverse] .= false
    fi_diedges_rev_df = select(fi_df, [:gene1, :gene2, :score, :dest_to_src, :direction, :source])
    rename!(fi_diedges_rev_df, :dest_to_src => :is_valid, :gene1 => :gene2, :gene2 => :gene1, :source => :diedge_type)
    fi_diedges_rev_df[!, :is_reverse] .= true
    fi_diedges_df = vcat(fi_diedges_fwd_df, fi_diedges_rev_df)
    filter!(r -> r.is_valid, fi_diedges_df)
    verbose && @info("$(nrow(fi_diedges_df)) diedge(s)")

    return fi_diedges_df
end

# calculate average probability to stay in the original direct neighborhood of a vertex
# or travel further
function walkweight_stats(digraph::AbstractGraph;
                          restart_probs::AbstractVector{Float64}=0.05:0.05:0.95)
    walkmtxs = HHN.random_walk_matrix.(Ref(digraph), restart_probs)
    vtx_totprobs = vec.(sum.(walkmtxs, dims=1)) - diag.(walkmtxs)
    vtx_neiprobs = HHN.neighborhood_weights.(walkmtxs, Ref(digraph))
    vtx_extprobs = vtx_totprobs .- vtx_neiprobs

    stats_df = vcat(DataFrame(restart_prob = restart_probs,
                              vertices = "neighbors",
                              trans_probs = sum.(vtx_neiprobs)./nv(digraph)),
                    DataFrame(restart_prob = restart_probs,
                              vertices = "exterior",
                              trans_probs = sum.(vtx_extprobs)./nv(digraph)))
    categorical!(stats_df)
    return stats_df
end

scctree_matrix!(treematrix::AbstractMatrix{W},
                walkmatrix::AbstractMatrix{W},
                vertex_weights::AbstractVector{W}) where W =
    mul!(treematrix, walkmatrix, Diagonal(vertex_weights))

scctree_matrix(walkmatrix::AbstractMatrix, vertex_weights::AbstractVector) =
    scctree_matrix!(similar(walkmatrix), walkmatrix, vertex_weights)

function vertex2label_flows(vertex_flows::AbstractVector{Tuple{Int, Int}},
                            vertices_labels::AbstractVector;
                            buffer::AbstractVector{Tuple{Int, String}})
    flows = isnothing(buffer) ? Vector{String}() : empty!(buffer)
    for (v, n) in vertex_flows
        push!(flows, (n, vertices_labels[v]))
    end
    sort!(flows)
    return join(string.(last.(flows), '(', first.(flows), ')'), ' ')
end

function vertex2label_flows(vertex_flows::AbstractVector{Int},
                            vertices_labels::AbstractVector;
                            buffer::AbstractVector{String})
    flows = isnothing(buffer) ? Vector{String}() : empty!(buffer)
    for v in vertex_flows
        push!(flows, vertices_labels[v])
    end
    sort!(flows)
    return join(flows, ' ')
end

fix_reactomefi_iactiontype(type::AbstractString) =
    replace(replace(string(type), ">" => "&gt;"), "<" => "&lt;")

fix_reactomefi_iactiontype!(types::CategoricalVector) =
    recode!(types,
            (levels(types) .=>
                fix_reactomefi_iactiontype.(levels(types)))...)

function flowgraphml(tree::HHN.SCCTree, threshold::Number;
    walkmatrix::AbstractMatrix,
    vertices_weights::AbstractVector,
    vertices_labels::AbstractVector,
    vertices_info::Union{AbstractDataFrame, Nothing}=nothing,
    source_threshold::Number=1.5,
    sinks::AbstractVector,
    component_groups::Bool=true,
    flow_edges::Bool=false,
    verbose::Bool=false,
    layout_cut_threshold::Union{Number, Nothing}=1.5*threshold,
    step_threshold::Number=1.5*threshold,
    step_sinks::Union{AbstractVector, AbstractSet, Nothing}=isnothing(sinks) ? nothing : Set(sinks),
    extra_node_attrs::Union{AbstractVector, Nothing}=nothing,
    edge_pvalue_max::Union{Number, Nothing}=nothing,
    mincompsize::Integer=5,
    flow_metrics::Union{AbstractDataFrame, Nothing}=nothing,
    kwargs...
)
    nflows = nothing
    if !isnothing(flow_metrics)
        thresh_metrics = filter(r -> r.threshold == threshold, flow_metrics)
        if nrow(thresh_metrics) > 0
            @assert nrow(thresh_metrics) == 1
            nflows = thresh_metrics.nflows[1]
        end
    end
    if isempty(sinks) || (!isnothing(nflows) && nflows == 0) # sinks/flows are empty, so we fake the sinks to display big SCCs
        verbose && @info("No sinks, exporting all SCCs with size >=$(mincompsize)")
        used_sinks = eachindex(vertices_weights) # sink is everything
        used_mincompsize = mincompsize
    else
        used_sinks = sinks
        used_mincompsize = nothing
    end
    graph_def = HHN.export_flowgraph(tree, threshold, walkmatrix,
                findall(>=(source_threshold), vertices_weights), used_sinks;
                flow_edges=flow_edges, verbose=verbose,
                step_threshold=step_threshold, step_sinks=step_sinks,
                exported_sinks=sinks,
                mincompsize=used_mincompsize,
                kwargs...)

    vertices_df = graph_def.vertices
    vertices_df.node = string.("v", vertices_df.vertex)
    vertices_df.component_node = string.("g", vertices_df.component)
    vertices_df.vertex_type = ifelse.(vertices_df.is_source,
                                      ifelse.(vertices_df.is_sink, "source/sink", "source"),
                                      ifelse.(vertices_df.is_sink, "sink",
                                      ifelse.(vertices_df.is_steptrace, "step", "gene")))
    vertices_df.flows_to_genes = missings(String, nrow(vertices_df))
    vertices_df.flows_from_genes = missings(String, nrow(vertices_df))
    vertices_df.loops_through_genes = missings(String, nrow(vertices_df))

    flowsbuf = Vector{Tuple{Int, String}}()
    loopsbuf = Vector{String}()
    for r in eachrow(vertices_df)
        if !ismissing(r.flows_to)
            r.flows_to_genes = vertex2label_flows(r.flows_to, vertices_labels, buffer=flowsbuf)
        end
        if !ismissing(r.flows_from)
            r.flows_from_genes = vertex2label_flows(r.flows_from, vertices_labels, buffer=flowsbuf)
        end
        if !ismissing(r.loops_through)
            r.loops_through_genes = vertex2label_flows(r.loops_through, vertices_labels, buffer=loopsbuf)
        end
    end
    if !isnothing(vertices_info)
        vertices_df = leftjoin(vertices_df, vertices_info, on=:vertex)
    end

    # FA3 layout-specific attributes
    # FIXME remove it after the attributes could be appended to GraphML
    vertices_df[!, :layout_size] .= 0.0
    vertices_df[!, :layout_mass] .= 1.0
    vertices_df[!, :layout_x] .= 0.0
    vertices_df[!, :layout_y] .= 0.0
    if !isnothing(layout_cut_threshold)
        layout_conncomps = HHN.cut(tree, layout_cut_threshold)
        layout_conncomp_ixs = fill(0, HHN.nelems(layout_conncomps))
        @inbounds for (i, comp) in enumerate(layout_conncomps)
            layout_conncomp_ixs[comp] .= i
        end
        vertices_df.layout_component = layout_conncomp_ixs[vertices_df.vertex]
        if (nrow(vertices_df) > 0) && component_groups
            allowmissing!(vertices_df)
            compnodes_df = repeat(vertices_df[1:1, :], nrow(graph_def.components))
            compnodes_df[!, :] .= missing
            compnodes_df.node = string.("g", 1:nrow(compnodes_df))
            compnodes_df[!, :vertex_type] = "connected_component"
            append!(vertices_df, compnodes_df)
        end
    end

    edges_df = filter(r -> (r.source != r.target) &&
            (flow_edges || r.has_walk || r.has_original || r.has_original_rev),
            graph_def.edges)
    edges_df.source_node = string.("v", edges_df.source)
    edges_df.target_node = string.("v", edges_df.target)
    edges_df.edge_type =
        ifelse.(edges_df.has_original .| edges_df.has_original_rev,
                ifelse.(edges_df.has_trace, "trace", "interaction"),
        ifelse.(edges_df.has_walk, "walk", ifelse.(edges_df.has_flow, "flow", missing)))
    edges_df.mlog10_prob_perm_walkweight_greater = .-log10.(min.(coalesce.(edges_df.prob_perm_walkweight_greater, 1.0),
                                                                 coalesce.(edges_df.prob_perm_walkweight_greater_rev, 1.0),
                                                                 1.0))
    edges_df[!, :layout_weight] .= 1.0 # FIXME remove when it's possible to add GraphML attrs
    for fpcol in [:flowpaths, :flowpaths_rev]
        hasproperty(edges_df, fpcol) || continue
        # replace the lists of vertices along the paths with "→x→y→|→z→" string
        # the source and the target are skipped to avoid redundancy
        pathsbuf = IOBuffer()
        newflowpaths = missings(String, nrow(edges_df))
        for (i, paths) in enumerate(edges_df[!, fpcol])
            ismissing(paths) && continue
            for path in paths
                isempty(path) && continue # skip direct paths
                (position(pathsbuf) > 0) && print(pathsbuf, '|')
                @inbounds for v in path
                    print(pathsbuf, '→') # FIXME show real edge type
                    v ∈ eachindex(vertices_labels) ? print(pathsbuf, vertices_labels[v]) : print(pathsbuf, '?')
                end
                print(pathsbuf, '→')
            end
            (position(pathsbuf) > 0) && (newflowpaths[i] = String(take!(pathsbuf)))
        end
        edges_df[!, fpcol] = newflowpaths
    end
    if !isnothing(edge_pvalue_max)
        filter!(r -> r.prob_perm_walkweight_greater <= edge_pvalue_max, edges_df)
    end

    node_attrs = [:gene_name, :protein_description,
                  :vertex_type, :weight, :walkweight,
                  :is_hit, :prob_perm_walkweight_greater,
                  :flows_to_genes, :flows_from_genes, :loops_through_genes,
                  :layout_component, :layout_size, :layout_mass, :layout_x, :layout_y]
    isnothing(extra_node_attrs) || append!(node_attrs, extra_node_attrs)

    edge_attrs = [:edge_type, :interaction_type, :target_type, :source_type,
                  :walkweight, :walkweight_rev,
                  :walkpermweight_median, :walkpermweight_median_rev,
                  :walkpermweight_mad, :walkpermweight_mad_rev,
                  :walkpermweight_mean, :walkpermweight_mean_rev,
                  :walkpermweight_std, :walkpermweight_std_rev,
                  :prob_perm_walkweight_greater, :prob_perm_walkweight_greater_rev,
                  :has_original, :has_original_rev, :has_flow, :has_walk, :flowpaths, :flowpaths_rev,
                  :mlog10_prob_perm_walkweight_greater, :layout_weight]
    graph = GraphML.import_graph(vertices_df, edges_df,
                                 node_col=:node, parent_col=component_groups ? :component_node : nothing,
                                 source_col=:source_node, target_col=:target_node,
                                 node_attrs=intersect(node_attrs, propertynames(vertices_df)),
                                 edge_attrs=intersect(edge_attrs, propertynames(edges_df)),
                                 verbose=verbose)
    return graph
end

function layout_flowgraph!(graph::GraphML.Graph; scale::Number=40.0, progressbar::Bool=true)
    vertices_df = graph.node_data
    edges_df = graph.edge_data

    # update vertex layout attributes
    vertices_df.layout_size = ifelse.((vertices_df.vertex_type .== "gene") .|
                                      (vertices_df.vertex_type .== "step"), 0.2, 0.4)*(1/nrow(edges_df))
    vertices_df[!, :layout_mass] .= 1.0
    vertices_df[!, :layout_x] .= 0.0
    vertices_df[!, :layout_y] .= 0.0
    vertices_df[!, :anyoriginal_edge] .= false
    v2r = Dict(r.vertex => i for (i, r) in enumerate(eachrow(vertices_df)))

    # check which vertices are connected by original nodes and which not
    # and put the node degree information into adjacent edges
    sinkset = Set(vertices_df.vertex[vertices_df.is_sink])
    edges_df.has_sink = (edges_df.source .∈ Ref(sinkset)) .| (edges_df.target .∈ Ref(sinkset))
    vertices_df[!, :ninedges] .= 0
    vertices_df[!, :noutedges] .= 0
    for inedges_df in groupby(edges_df, :target)
        trg_row = v2r[inedges_df.target[1]]
        vertices_df[trg_row, :ninedges] = nrow(inedges_df)
        if any(r -> r.has_original || r.has_original_rev, eachrow(inedges_df))
            vertices_df[trg_row, :anyoriginal_edge] = true
        end
    end
    for outedges_df in groupby(edges_df, :source)
        src_row = v2r[outedges_df.source[1]]
        vertices_df[src_row, :noutedges] = nrow(outedges_df)
        if any(r -> r.has_original || r.has_original_rev, eachrow(outedges_df))
            vertices_df[src_row, :anyoriginal_edge] = true
        end
    end
    edges_df[!, :ends_have_original] .= false
    edges_df[!, :source_nedges] .= 0
    edges_df[!, :target_nedges] .= 0
    for i in 1:nrow(edges_df)
        src_row = v2r[edges_df.source[i]]
        trg_row = v2r[edges_df.target[i]]
        if vertices_df.anyoriginal_edge[src_row] &&
           vertices_df.anyoriginal_edge[trg_row]
           edges_df[i, :ends_have_original] = true
        end
        edges_df[i, :source_nedges] = vertices_df.ninedges[src_row] + vertices_df.noutedges[src_row]
        edges_df[i, :target_nedges] = vertices_df.ninedges[trg_row] + vertices_df.noutedges[trg_row]
    end
    # calculate weights depending on which vertices it connects
    weights = max.(coalesce.(edges_df.walkweight, 0.0), coalesce.(edges_df.walkweight_rev, 0.0))
    edges_df.layout_weight .= clamp.(weights * inv(mean(filter(>(0), weights))), 0.01, 5.0).^(0.25) .*
        ifelse.(edges_df.has_sink, ifelse.(edges_df.has_original .| edges_df.has_original_rev, 10.0,
                                           10.0 ./ sqrt.(min.(edges_df.source_nedges, edges_df.target_nedges))),
        ifelse.(edges_df.has_original .| edges_df.has_original_rev, 3.0,
        ifelse.(edges_df.ends_have_original,
                0.5 ./ min.(edges_df.source_nedges, edges_df.target_nedges).^2, # the vertices are connected by original edges, so weight this edge down
                1.0 ./ sqrt.(min.(edges_df.source_nedges, edges_df.target_nedges)))))

    # start FA3 layout
    fa_graph = FA.Graph(edges_df,
                        filter(r -> !ismissing(r.vertex), vertices_df),
                        src_node_col=:source, dest_node_col=:target, weight_col=:layout_weight,
                        node_col=:vertex, size_col=:layout_size, mass_col=:layout_mass)
    if hasproperty(vertices_df, :layout_component)
        fa_node_dislikes = Matrix{Float64}(undef, (length(fa_graph.nodes), length(fa_graph.nodes)))
        for i in CartesianIndices(fa_node_dislikes)
            fa_node_dislikes[i] = ifelse(vertices_df.layout_component[i[1]] == vertices_df.layout_component[i[2]], 0.0, 1.0)
        end
    else
        fa_node_dislikes = nothing
    end
    #=
    fa_node_likes = FA.socioaffinity(SimpleGraph(fa_graph), p=(0.5, 0.0), q=0.0)
    for (i, e) in enumerate(fa_graph.edges)
        elike = fa_node_likes[e.dest, e.src]
        if elike > 0
            fa_graph.edges[i] = FA.Edge(e.src, e.dest, e.weight / elike)
        end
    end
    =#

    FA.layout!(fa_graph, FA.ForceAtlas3Settings(fa_graph,
               outboundAttractionDistribution=false,
               attractionStrength=3.0, attractionEdgeWeightInfluence=0.1, jitterTolerance=0.1,
               repulsionStrength=isnothing(fa_node_dislikes) ? 4.0 : 2 .* (1.0 .+ fa_node_dislikes),
               repulsionFalloff=2.25,
               repulsionNodeModel=:Point,
               gravity=0.1, gravityFalloff=1.0, gravityShape=:Rod,
               gravityRodCorners=((-4.0, 0.0), (4.0, 0.0)), gravityRodCenterWeight=0.1),
               nsteps=1000, progressbar=progressbar)
    FA.layout!(fa_graph, FA.ForceAtlas3Settings(fa_graph,
               outboundAttractionDistribution=false,
               attractionStrength=3.0, attractionEdgeWeightInfluence=1.0, jitterTolerance=0.1,
               repulsionStrength = isnothing(fa_node_dislikes) ? 5.0 : 4 .+ 2 * fa_node_dislikes,
               repulsionFalloff = 2.25,
               repulsionNodeModel=:Circle,
               gravity=0.75, gravityFalloff=2, gravityShape=:Rod,
               gravityRodCorners=((-4.0, 0.0), (4.0, 0.0)), gravityRodCenterWeight=0.1),
               nsteps=5000, progressbar=progressbar)

    x_ix = findfirst(attr -> attr.name == "layout_x", graph.node_attrs)
    y_ix = findfirst(attr -> attr.name == "layout_y", graph.node_attrs)
    graph.node_attrs[x_ix].values[eachindex(fa_graph.nodes)] .= FA.extract_layout(fa_graph)[1] .* scale
    graph.node_attrs[y_ix].values[eachindex(fa_graph.nodes)] .= FA.extract_layout(fa_graph)[2] .* scale

    return graph
end

end
