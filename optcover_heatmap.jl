module OptCoverHeatmap

using PlotlyJS, Clustering, Distances, OptEnrichedSetCover, TextWrap, DataFrames
using Printf: @sprintf

using ..FrameUtils: frame2array

"""
    axis_order(mtx::AbstractMatrix, dimaxis::AbstractDataFrame;
               dims::Integer=1, order::Any=:hclu)

Order the `dims` axis of `mtx` matrix using `order` method.

The following `order` methods are supported:
  * `nothing`: keep the order "as is"
  * `:hclu` (the default): use hierarchical clustering with Ward
    linkage method and cosine distance between `mtx` columns or rows
    (see [`clustering`](@ref Clustering.hclust))
  * explicit order specified by a vector of indices
  * function that takes the `dimaxis` metadata data frame and returns
    the vector of element indices
"""
function axis_order(mtx::AbstractMatrix, dimaxis::AbstractDataFrame; dims::Integer=1, order::Any=:hclu)
    if isnothing(order)
        return axes(mtx, dims)
    elseif order == :hclu
        dist = pairwise(CosineDist(), mtx, dims=dims)
        dist[isnan.(dist)] .= 100.0
        hclu = hclust(dist, linkage=:ward)
        return hclu.order
    elseif order isa AbstractVector
        @assert length(order) == size(mtx, dims)
        return order
    elseif order isa Function
        return order(dimaxis)# collect(1:size(scores_mtx, 2))
    else
        throw(ArgumentError("Don't know how to handle $(dims == 2 ? "columns" : "rows") order: $order"))
    end
end

"""
    heatmap_matrices(covers_df;
                     experiment_axis_title="Cluster", term_axis_title="Term",
                     elements_label="genes",
                     experiment_cols = [:cluster, :cluster_genes, :ncluster],
                     term_cols = [:term_id, :term_collection, :term_genes, :term_name, :term_descr, :set_relevance],
                     process_experiment_axis::Function = df -> (df, df[experiment_cols[1]], df[experiment_cols[1]]),
                     process_term_axis::Function = df -> (df, df[term_cols[1]], df[term_cols[1]]),
                     term_order=:hclu, experiment_order=:hclu)

Prepares `covers_df` (output of `OptEnrichedSetCover.jl`) for display as an interactive heatmap.
"""
function heatmap_matrices(covers_df;
    experiment_axis_title="Cluster", term_axis_title="Term", elements_label="genes",
    experiment_cols = [:cluster, :cluster_genes, :ncluster],
    term_cols = [:term_id, :term_collection, :term_genes, :term_name, :term_descr, :set_relevance],
    process_experiment_axis::Function = df -> (df, df[experiment_cols[1]], df[experiment_cols[1]]),
    process_term_axis::Function = df -> (df, df[term_cols[1]], df[term_cols[1]]),
    term_order=:hclu, experiment_order=:hclu)

    weight_col = "set_overlap_log10pvalue" ∈ names(covers_df) ?
            :set_overlap_log10pvalue : :set_original_weight

    if "nmasked" ∈ names(covers_df)
        hasmasks = true
        covers_df.ntotal = covers_df.nmasked .+ covers_df.nunmasked
        extra_term_cols = [:ntotal]
        nmasked_mtx, _ = frame2array(covers_df, [term_cols, experiment_cols],
                                     data_col=:nmasked, default=0)
        isect_mtx, _ = frame2array(covers_df, [term_cols, experiment_cols],
                                   data_col=:intersect_genes, default="")
    else
        hasmasks = false
        extra_term_cols = Symbol[]
        nmasked_mtx = nothing
        isect_mtx = nothing
    end
    scores_mtx, (terms_axis, experiments_axis) = frame2array(covers_df, [[term_cols; extra_term_cols], experiment_cols],
                                                       data_col=weight_col, default=NaN)
    isnothing(nmasked_mtx) || (scores_mtx[nmasked_mtx .== 0] .= NaN)
    weights_mtx = frame2array(covers_df, [term_cols, experiment_cols],
                              data_col=:set_weight, default=NaN)[1]
    # prepare axis metadata, breakpoint labels and tooltips
    experiments_df, experiment_axis_labels, experiment_tips = process_experiment_axis(experiments_axis)
    experiments_df.axis_label = experiment_axis_labels
    experiments_df.axis_tip = experiment_tips
    terms_df, term_axis_labels, term_tips = process_term_axis(terms_axis)
    terms_df.axis_label = term_axis_labels
    terms_df.axis_tip = term_tips
    # generate tooltips for each matrix element
    tips_mtx = Matrix{String}(undef, size(scores_mtx))
    for i in 1:size(scores_mtx, 1), j in 1:size(scores_mtx, 2)
        tips_mtx[i, j] = string(experiment_axis_title, ": ", experiments_df.axis_tip[j],
                    " (", elements_label, "=", experiments_df[j, experiment_cols[end]], ")<br>\n",
                term_tips[i], ":<br>\n",
                    " ", hasmasks ? string(elements_label, "=", terms_df[i, :ntotal], ", ") : "",
                    "relevance=", @sprintf("%.3f", terms_df[i, :set_relevance]), "<br>\n",
                !ismissing(scores_mtx[i, j]) ? @sprintf("<b>P-value</b>=%.4e", exp10(scores_mtx[i, j])) :
                    "<b>P-value</b> missing", "<br>\n",
                hasmasks ? string(
                "<b>Common ", elements_label, " (", nmasked_mtx[i, j], ")</b>: ",
                    !ismissing(isect_mtx[i, j]) ? replace(wrap(isect_mtx[i, j], width=50,
                                 subsequent_indent="    "), r"\n" => "<br>\n") : "<none>", "<br>\n"
                ) : "",
                @sprintf("cover weight=%.3f", weights_mtx[i, j]))
    end

    hclust_mtx = ifelse.(isnan.(coalesce.(scores_mtx, NaN)), 0.0, (-scores_mtx).^0.5)
    term_ixs = axis_order(hclust_mtx, terms_df, dims=1, order=term_order)
    experiment_ixs = axis_order(hclust_mtx, experiments_df, dims=2, order=experiment_order)

    tips_mtx = ifelse.(isnan.(coalesce.(scores_mtx, NaN)), missing, tips_mtx)[term_ixs, experiment_ixs] # before scores_mtx is modified!
    scores_mtx = scores_mtx[term_ixs, experiment_ixs]
    terms_df = terms_df[term_ixs, :]
    experiments_df = experiments_df[experiment_ixs, :]
    return scores_mtx, tips_mtx, terms_df, experiments_df
end

function oesc_heatmap(covers_df;
    colorscale::AbstractString = "Blackbody", reversescale::Bool=true,
    zmax=0.0, zmin=-10,
    margin_b = nothing, margin_l = nothing, margin_t = 5, margin_r = 5,
    cell_width = nothing, cell_height = nothing,
    paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
    gridcolor = "#888", gridwidth = 1, transpose::Bool = false,
    matrixkwargs...
)
    scores_mtx, tips_mtx, terms_df, experiments_df = heatmap_matrices(covers_df; matrixkwargs...)
    if transpose
        xaxis_df = terms_df
        yaxis_df = experiments_df
    else
        scores_mtx = scores_mtx' # plotly expets matrix in row-major order, so we have to transpose this one
        tips_mtx = permutedims(tips_mtx, (2, 1))
        xaxis_df = experiments_df
        yaxis_df = terms_df
    end
    if any(nonunique(select(xaxis_df, :axis_label)))
        @warn "Non unique X axis labels: $(join(xaxis_df.axis_label[nonunique(select(xaxis_df, :axis_label))], ", "))"
    end
    if any(nonunique(select(yaxis_df, :axis_label)))
        @warn "Non unique Y axis labels: $(join(yaxis_df.axis_label[nonunique(select(yaxis_df, :axis_label))], ", "))"
    end
    res = plot(heatmap(z=scores_mtx, text=tips_mtx, hoverinfo="text",
            x=xaxis_df.axis_label, y=yaxis_df.axis_label,
            colorscale=colorscale, reversescale=reversescale,
            zmax=zmax, zmin=zmin, zauto=false, zsmooth=false, connectgaps=false,
            hoverongaps=false, transpose=false,
            hoverlabel=attr(bgcolor="#DDD", font=attr(color="#000")),
            xtype="array", ytype="array", xgap=gridwidth, ygap=gridwidth),
            Layout(paper_bgcolor=paper_bgcolor,
                   plot_bgcolor=plot_bgcolor,
                   modebar_bgcolor="#FFF",
                   xaxis = attr(type="category", tickson="boundaries", gridcolor=gridcolor, gridwidth=gridwidth),
                   yaxis = attr(type="category", tickson="boundaries", gridcolor=gridcolor, gridwidth=gridwidth)),
            config=PlotConfig(scrollZoom = false))
    if margin_l !== nothing
        res.plot.layout[:margin_l] = margin_l
    else
        margin_l = 0
    end
    if margin_r !== nothing
        res.plot.layout[:margin_r] = margin_r
    else
        margin_r = 0
    end
    if margin_b !== nothing
        res.plot.layout[:margin_b] = margin_b
    else
        margin_b = 0
    end
    if margin_t !== nothing
        res.plot.layout[:margin_t] = margin_t
    else
        margin_t = 0
    end
    if cell_width !== nothing
        res.plot.layout[:width] = cell_width * nrow(xaxis_df) + margin_l + margin_r
    end
    if cell_height !== nothing
        res.plot.layout[:height] = cell_height * nrow(yaxis_df) + margin_t + margin_b
    end
    return res
end

const ChangeStyle = Dict(
    "0" => (label="0", color="gray"),
    "+" => (label="&#9650;", color="red"),
    "-" => (label="&#9660;", color="blue"),
)

stylize_change(str) = stylize_identifier(str, ChangeStyle)

function stylize_identifier(id::AbstractString, id_info::AbstractDict{<:AbstractString})
    if haskey(id_info, id)
        info = id_info[id]
        return "<span style=\"color: $(info.color);\">$(info.label)</span>"
    else
        @warn "No info found for $id"
        return id
    end
end

stylize_treatment(id::AbstractString, treatment_info::AbstractDict{<:AbstractString}) =
    stylize_identifier(id, treatment_info)

function stylize_contrast(contrast::AbstractString, condition_info::AbstractDict{<:AbstractString})
    contrast_match = match(r"^(?<lhs>.+)_vs_(?<rhs>[^+-]+)(?<sign>[+-])?$", contrast)
    if (isnothing(contrast_match))
        @warn "No format match for contrast '$(contrast)'"
        return contrast
    else
        res = stylize_treatment(contrast_match[:lhs], condition_info) * " <b>vs</b> " *
              stylize_treatment(contrast_match[:rhs], condition_info);
        if (!isnothing(contrast_match[:sign]))
            res *= stylize_change(contrast_match[:sign]);
        end
        return res
    end
end

stylize_condition(str) = foldl(replace, [
    r"SC35M([^ +]*)" => s"SC35M(\1)",
    "()" => "(WT)",
    "SC35M(WT)" => "<span style=\"color: red;\">SC35M(<b>WT</b>)</span>",
    "SC35M(delNS1)" => "<span style=\"color: orange;\">SC35M(<b>&#x394;NS1</b>)</span>",
    "Mock" => "<span style=\"color: gray;\">Mock</span>"],
    init = str)

stylize_pulsechase(str) = foldl(replace, [
    "_cum" => "",
    "_scaled" => "<sub>scaled</sub>",
    r"native"i => "<span style=\"color: #3fbc97; font-weight: bold;\">native</span>",
    r"pulse"i => "<span style=\"color: #806fb2; font-weight: bold;\">pulse</span>",
    r"chase"i => "<span style=\"color: #488dcb; font-weight: bold;\">chase</span>"],
    init = str)

function process_cluster_axis(clu_df)
    clu_df,
    string.("c", clu_df.cluster),
    string.("cluster #", clu_df.cluster)
end

function process_effect_axis(eff_df)
    eff_df,
    string.(coalesce.(eff_df.effect_label, eff_df.effect), "&nbsp;",
            stylize_change.(eff_df.change)),
    string.("Effect: ", coalesce.(eff_df.effect_label, eff_df.effect), "&nbsp;",
            stylize_change.(eff_df.change))
end

function process_term_axis(terms_df)
    terms_df[!, :term_name] .= coalesce.(terms_df.term_name, terms_df.term_descr,
                                         terms_df.term_id, "<missing>")
    if startswith(string(terms_df.term_collection[1]), "Hs_SigDB_")
        term_axis_labels = replace.(terms_df.term_name, "_" => " ")
    else
        term_axis_labels = replace.(wrap.(terms_df.term_name, width=60),
                                    Ref(r"\n" => "<br>\n"))
    end
    nonunique_labels_mask = nonunique(select(terms_df, [:term_name], copycols=false))
    if any(nonunique_labels_mask) # workaround duplicate axis labels
        duplicate_labels_mask = term_axis_labels .∈ Ref(Set(term_axis_labels[nonunique_labels_mask]))
        term_axis_labels[duplicate_labels_mask] .= term_axis_labels[duplicate_labels_mask] .* " (" .* terms_df.term_id[duplicate_labels_mask] .* ")"
    end
    # workaround line breaks obscure the next line
    term_axis_labels = replace(term_axis_labels, r"(<br>\n)(.+)$" => s"\1<sup>\2</sup>", count=1)
    term_tip_labels = [string("Term ", r.term_id, " (",
                     replace(wrap(r.term_name, width=60,
                                  subsequent_indent="    "),
                             r"\n" => "<br>\n"), ")")
                        for r in eachrow(terms_df)]
    return terms_df, term_axis_labels, term_tip_labels
end

end
