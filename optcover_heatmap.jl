module OptCoverHeatmap

using PlotlyJS, Clustering, Distances, OptEnrichedSetCover, TextWrap, DataFrames
using Printf: @sprintf

using ..FrameUtils: frame2array

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

function heatmap_matrices(covers_df;
    mask_axis_title="Cluster", term_axis_title="Term", elements_label="genes",
    mask_cols = [:cluster, :cluster_genes, :ncluster],
    term_cols = [:term_id, :term_collection, :term_genes, :term_name, :term_descr, :set_relevance],
    process_mask_axis::Function = df -> (df, df[mask_cols[1]], df[mask_cols[1]]),
    process_term_axis::Function = df -> (df, df[term_cols[1]], df[term_cols[1]]),
    term_order=:hclu, mask_order=:hclu)

    covers_df.ntotal = covers_df.nmasked .+ covers_df.nunmasked
    scores_mtx, (terms_axis, masks_axis) = frame2array(covers_df, [[term_cols; :ntotal], mask_cols],
                                                       data_col=:set_overlap_log10pvalue,
                                                       default=NaN)
    nmasked_mtx = frame2array(covers_df, [term_cols, mask_cols],
                              data_col=:nmasked, default=0)[1]
    scores_mtx[nmasked_mtx .== 0] .= NaN
    isect_mtx = frame2array(covers_df, [term_cols, mask_cols],
                            data_col=:intersect_genes, default="")[1]
    weights_mtx = frame2array(covers_df, [term_cols, mask_cols],
                              data_col=:set_weight, default=NaN)[1]
    masks_df, mask_axis_labels, mask_tips = process_mask_axis(masks_axis)
    masks_df.axis_label = mask_axis_labels
    masks_df.axis_tip = mask_tips
    terms_df, term_axis_labels, term_tips = process_term_axis(terms_axis)
    terms_df.axis_label = term_axis_labels
    terms_df.axis_tip = term_tips
    tips_mtx = Matrix{String}(undef, size(scores_mtx))
    for i in 1:size(scores_mtx, 1), j in 1:size(scores_mtx, 2)
        tips_mtx[i, j] = string(mask_axis_title, ": ", masks_df.axis_tip[j],
                    " (", elements_label, "=", masks_df[j, mask_cols[end]], ")<br>\n",
                term_tips[i], ":<br>\n",
                    " ", elements_label, "=", terms_df[i, :ntotal],
                    ", relevance=", @sprintf("%.3f", terms_df[i, :set_relevance]), "<br>\n",
                @sprintf("<b>P-value</b>=%.4e", exp10(scores_mtx[i, j])), "<br>\n",
                "<b>Common ", elements_label, " (", nmasked_mtx[i, j], ")</b>: ",
                    !ismissing(isect_mtx[i, j]) ? replace(wrap(isect_mtx[i, j], width=50,
                                 subsequent_indent="    "), r"\n" => "<br>\n") : "<none>", "<br>\n",
                @sprintf("cover weight=%.3f", weights_mtx[i, j]))
    end

    hclust_mtx = ifelse.(isnan.(scores_mtx), 0.0, (-scores_mtx).^0.5)
    term_ixs = axis_order(hclust_mtx, terms_df, dims=1, order=term_order)
    mask_ixs = axis_order(hclust_mtx, masks_df, dims=2, order=mask_order)

    tips_mtx = ifelse.(isnan.(scores_mtx), missing, tips_mtx)[term_ixs, mask_ixs] # before scores_mtx is modified!
    scores_mtx = scores_mtx[term_ixs, mask_ixs]
    terms_df = terms_df[term_ixs, :]
    masks_df = masks_df[mask_ixs, :]
    return scores_mtx, tips_mtx, terms_df, masks_df
end

function oesc_heatmap(covers_df;
    colorscale::AbstractString = "Blackbody", reversescale::Bool=true,
    zmax=0.0, zmin=-10,
    margin_b = nothing, margin_l = nothing, cell_width = nothing, cell_height = nothing,
    paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
    gridcolor = "#888", gridwidth = 1, transpose::Bool = false,
    matrixkwargs...
)
    scores_mtx, tips_mtx, terms_df, masks_df = heatmap_matrices(covers_df; matrixkwargs...)
    if transpose
        xaxis_df = terms_df
        yaxis_df = masks_df
    else
        scores_mtx = scores_mtx' # plotly expets matrix in row-major order, so we have to transpose this one
        tips_mtx = permutedims(tips_mtx, (2, 1))
        xaxis_df = masks_df
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
            hoverongaps=false,
            hoverlabel=attr(bgcolor="#DDD", font=attr(color="#000")),
            xtype="array", ytype="array", xgap=gridwidth, ygap=gridwidth),
            Layout(paper_bgcolor=paper_bgcolor,
                   plot_bgcolor=plot_bgcolor,
                   modebar_bgcolor="#FFF",
                   xaxis = attr(tickson="boundaries", gridcolor=gridcolor, gridwidth=gridwidth),
                   yaxis = attr(tickson="boundaries", gridcolor=gridcolor, gridwidth=gridwidth)))
    if margin_l !== nothing
        res.plot.layout[:margin_l] = margin_l
        if cell_width !== nothing
            res.plot.layout[:width] = margin_l + cell_width * nrow(xaxis_df) + 120
        end
    end
    if margin_b !== nothing
        res.plot.layout[:margin_b] = margin_b
        if cell_height !== nothing
            res.plot.layout[:height] = margin_b + cell_height * nrow(yaxis_df) + 50
        end
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
