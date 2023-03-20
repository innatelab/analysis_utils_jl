module PlotlyUtils

using Printf: @sprintf
using PlotlyBase, PlotlyJS, Colors, DataFrames

"""
    variable_rows_layout(col2nrow::AbstractVector{<:Integer}) -> PlotlyJS.Layout

Creates plotly layout with regularly spaced ``N`` columns (`N = length(col2nrow)`)
and `col2nrow[i]` regularly spaced rows in each columns.
"""
function variable_rows_layout(col2nrow::AbstractVector{<:Integer})
    w, dx = PlotlyJS.subplot_size(1, length(col2nrow), false)[[1, 3]]
    out = Layout()
    x = 0.0
    subplot = 1
    for (col, nr) in enumerate(col2nrow)
        xdom = [x, x + w]::Vector{Float64}

        if nr > 0
            h, dy = PlotlyJS.subplot_size(nr, length(col2nrow), false)[[2, 4]]
            y = 1.0 # start from top
            for row in 1:nr
                out["xaxis$subplot"] = Dict([:domain=>copy(xdom), :anchor=>"y$subplot"])
                out["yaxis$subplot"] = Dict([:domain=>[y - h, y], :anchor=>"x$subplot"])
                #@show row col subplot sub_code x y w h
                subplot += 1
                y -= (dy + h)
            end
        end
        x += dx + w
    end
    return out
end

include(joinpath(@__DIR__, "color_utils.jl"))
include(joinpath(@__DIR__, "colorscales.jl"))
include(joinpath(@__DIR__, "data.jl"))
include(joinpath(@__DIR__, "tooltips.jl"))
include(joinpath(@__DIR__, "scatter.jl"))
include(joinpath(@__DIR__, "trace_registry.jl"))
include(joinpath(@__DIR__, "series.jl"))
include(joinpath(@__DIR__, "comparison_bracket.jl"))
include(joinpath(@__DIR__, "pulse_plots.jl"))
include(joinpath(@__DIR__, "hhotnet_plots.jl"))
include(joinpath(@__DIR__, "msinstrument_plots.jl"))

end
