module OptCoverPlots

using DataFrames, OptEnrichedSetCover, PlotlyJS
using OptEnrichedSetCover: MultiobjCoverProblemResult
using Printf: @sprintf

function plot_paretofront(cover::MultiobjCoverProblemResult;
                          plot_unfolded::Bool=false)
    fitfront = cover.folded_scores
    rawfitfront = cover.raw_scores
    best_folded_score = fitfront[best_index(cover)]
    zmin, zmax = length(first(fitfront)) == 3 ? extrema(getindex.(fitfront, 3)) : (0.0, 0.0)
    if zmin == zmax # 2d mode
        plots = [scatter(;x=getindex.(fitfront, 1), y=getindex.(fitfront, 2),
                         name="Pareto front",
                         text=[@sprintf("fitness=(%.3f, %.3f)<br>\nraw=(%.3f, %.3f, %.3f, %.3f)",
                                        fitfront[i]..., rawfitfront[i]...) for i in eachindex(fitfront)],
                         mode="lines+markers", hoverinfo="text", marker_size=5)]
        if plot_unfolded
            pushfirst!(plots,
                scatter(;x=getindex.(rawfitfront, 1), y=getindex.(rawfitfront, 2),
                        name="Raw Pareto front",
                        text=[@sprintf("fitness=(%.3f, %.3f)<br>\nraw=(%.3f, %.3f, %.3f, %.3f)",
                                       fitfront[i]..., rawfitfront[i]...) for i in eachindex(fitfront)],
                        mode="lines+markers", hoverinfo="text", marker_size=3))
        end
        push!(plots,
            scatter(;x=[best_folded_score[1]], y=[best_folded_score[2]],
                    name="Best",
                    text=[@sprintf("Best fitness=(%.3f, %.3f)<br>\nraw=(%.3f, %.3f, %.3f, %.3f)",
                                   best_folded_score..., rawfitfront[best_index(cover)]...)],
                    mode="markers", marker_size=8))
    else # 3d mode
        plots = [scatter3d(;x=getindex.(fitfront, 1), y=getindex.(fitfront, 2), z=getindex.(fitfront, 3),
                        name="Pareto front",
                        text=[@sprintf("fitness=(%.3f, %.3f, %.3f)<br>\nraw=(%.3f, %.3f, %.3f, %.3f)",
                                       fitfront[i]..., rawfitfront[i]...) for i in eachindex(fitfront)],
                        mode="markers", hoverinfo="text", marker_size=2)]
        if plot_unfolded
            pushfirst!(plots,
                scatter3d(;x=getindex.(rawfitfront, 1), y=getindex.(rawfitfront, 2), z=getindex.(rawfitfront, 3),
                          name="Raw Pareto front",
                          text=[@sprintf("fitness=(%.3f, %.3f, %.3f)<br>\nraw=(%.3f, %.3f, %.3f, %.3f)",
                                         fitfront[i]..., rawfitfront[i]...) for i in eachindex(fitfront)],
                          mode="markers", hoverinfo="text", marker_size=1.5)
            )
        end
        push!(plots,
            scatter3d(;x=[best_folded_score[1]], y=[best_folded_score[2]], z=[best_folded_score[3]],
                      name="Best",
                      text=[@sprintf("Best fitness=(%.3f, %.3f, %.3f)<br>\nraw=(%.3f, %.3f, %.3f, %.3f)",
                                     best_folded_score...,rawfitfront[best_index(cover)]...)],
                      mode="markers", marker_size=5))
    end
    return plot(plots, Layout(hovermode="closest"))
end

end
