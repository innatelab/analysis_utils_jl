module SubheatmapUtils

using DataFrames
using VegaLite

fixlabel(::Missing) = missing
fixlabel(s::AbstractString) = replace(s, r"</?[^>]+>" => "")

multilinetovector(::Missing) = missing
multilinetovector(s::AbstractString) =
    split(replace(s, r"<br>\n?|\n" => "\n"), '\n')

function subheatmap_frame(z::AbstractMatrix,
                          row_axis::AbstractDataFrame,
                          colXsub_axis::AbstractDataFrame;
                          row_label_col::Symbol = :axis_label,
                          col_label_col::Symbol = :axis_label,
                          col_sub_cols::AbstractVector{Symbol},
                          col_cols = Not(col_sub_cols),
                          extra_row_cols::Union{AbstractVector{Symbol}, Nothing} = nothing,
                          extra_col_cols::Union{AbstractVector{Symbol}, Nothing} = nothing,
                          tips::Union{AbstractMatrix,Nothing} = nothing)
    @assert nrow(row_axis) == size(z, 1)
    @assert nrow(colXsub_axis) == size(z, 2)
    row_axis = copy(row_axis)
    row_axis[!, row_label_col] = fixlabel.(row_axis[!, row_label_col])
    #row_axis[!, row_label_col] = multilinetovector.(row_axis[!, row_label_col])
    col_axis = select(colXsub_axis, col_cols) |> unique!
    col_axis.col_index = 1:nrow(col_axis)
    colXsub_axis = copy(colXsub_axis, copycols=false)
    colXsub_axis.orig_order = 1:nrow(colXsub_axis)
    colXsub_axis = innerjoin(colXsub_axis, col_axis, on=names(col_axis)[1:end-1])
    sort!(colXsub_axis, :orig_order)
    subkey_col = col_sub_cols[1]
    # HACK for unique(CatArray) returning Array{String}
    subkeys = sort!(unique!(select(colXsub_axis, [subkey_col]))[!, subkey_col])
    @show subkeys
    if colXsub_axis[!, subkey_col] isa CategoricalArray
        # HACK subkeys and subkey_col have different pools now
        colXsub_axis.sub_index = searchsortedfirst.(Ref(CategoricalArrays.refs(subkeys)),
                                                    CategoricalArrays.refs(colXsub_axis[!, subkey_col]))
    else
        colXsub_axis.sub_index = searchsortedfirst.(Ref(subkeys), colXsub_axis[!, subkey_col])
    end
    @assert nrow(colXsub_axis) <= length(subkeys)*nrow(col_axis) # check that col_axis + subkey are really a key for colXsub_axis
    #col_axis[!, col_label_col] = multilinetovector.(col_axis[!, col_label_col])
    col_axis[!, col_label_col] = fixlabel.(col_axis[!, col_label_col])
    res = DataFrame(value = copy(vec(z)),
            row_index = repeat(1:nrow(row_axis), outer=size(z, 2)),
            col_index = repeat(colXsub_axis.col_index, inner=size(z, 1)),
            sub_index = repeat(colXsub_axis.sub_index, inner=size(z, 1)))
    res.row_label = row_axis[res.row_index, row_label_col]
    res.col_label = col_axis[res.col_index, col_label_col]
    res.sub_label = subkeys[res.sub_index]
    if !isnothing(tips)
        @assert size(tips) == size(z)
        res.tooltip = multilinetovector.(vec(tips))
    end
    if !isnothing(extra_row_cols)
        for col in extra_row_cols
            res[!, "row_$col"] = repeat(row_axis[!, col], outer=size(z, 2))
        end
    end
    if !isnothing(extra_col_cols)
        for col in extra_col_cols
            res[!, "col_$col"] = repeat(colXsub_axis[!, col], inner=size(z, 1))
        end
    end
    return res, row_axis, col_axis
end

function DiagonalTriangles(; scale::Number=1, offset::Number=0.1)
    hm = (1-offset)*scale
    vm = scale
    l = (2-offset)*scale
    return ["M$(-hm) $(vm)v$(-l)h$(l)z", "M$(-hm) $(vm)h$(l)v$(-l)z"]
end

const HotColorScheme = ["rgb(0,0,0)", "rgb(230,0,0)", "rgb(255,210,0)", "rgb(255,255,255)"]

function vegalite_subheatmap(subheatmap_df::AbstractDataFrame;
                             cell_width=20, cell_height=20, gap_width=1, grid_width=0.25,
                             labelLimit_l=200, labelLimit_b=200,
                             subshapes::AbstractArray=DiagonalTriangles(),
                             value_domain=(-10, 0), colorscheme=HotColorScheme,
                             xaxis_label::AbstractString = "column",
                             xaxis_tick_col=:axis_label,
                             yaxis_label::AbstractString = "row",
                             yaxis_tick_col=:axis_label,
                             subaxis_label::AbstractString = "sub",
                             coloraxis_label::AbstractString = "value"
)
    subheatmap_plot = subheatmap_df |> @vlplot(
        config={view={stroke=nothing}},
        height={step=cell_height+gap_width},
        width={step=cell_width+gap_width},
        mark={:point, filled=true},
        row={"rowblock_index:n"}, column={"colblock_index:n"},
        x={"col_label:n",
            sort=:col_index,
            axis = {
                title=xaxis_label,
                tickBand = "extent",
                grid = true,
                gridWidth = grid_width,
                labelLimit = labelLimit_b,
                labelExpr = "split(replace(datum.label, /<br>[\\n]?|[\\n]/, '\\n'), '\\n')",
                labelAlign = "left", labelBaseline = "middle", labelAngle = 45,
            },
            #scale = {domain = cols_df[!, xaxis_tick_col]}
        },
        y={"row_label:n",
            sort=:row_index,
            axis = {
                labelLimit = labelLimit_l,
                labelExpr = "split(replace(datum.label, /<br>[\\n]?|[\\n]/, '\\n'), '\\n')",
                labelOffset = "(1-length(split(replace(datum.label, /<br>[\\n]?|[\\n]/, '\\n'), '\\n')))*5",
                title = yaxis_label,
                titleX = -labelLimit_l-5,
                tickBand = "extent",
                grid = true,
                gridWidth = grid_width,
            },
            #scale = {domain = rows_df[!, yaxis_tick_col]}
        },
        shape={
            "sub_label:n", header={title=subaxis_label},
            sort={field="sub_index", order="descending"},
            scale={range=subshapes}, #domain=sublabels,
        },
        color={
            "value:q",
            scale={domain=value_domain, scheme=colorscheme, interpolate=:rgb},
            legend={title=coloraxis_label, titleOrient="right"},
        },
        tooltip={
            "tooltip:n"
        },
        resolve={scale={x="independent"}},
        opacity={ value=1 },
        size={ value=400 }
    )
    # HACK remove properties that are not present in the data
    hasproperty(subheatmap_df, :rowblock_index) || delete!(subheatmap_plot.encoding.params, "row")
    hasproperty(subheatmap_df, :collock_index) || delete!(subheatmap_plot.encoding.params, "column")
    hasproperty(subheatmap_df, :tooltip) || delete!(subheatmap_plot.encoding.params, "tooltip")
    return subheatmap_plot
end

end
