function permstats_plot(realstats_df::AbstractDataFrame,
                        permstats_df::AbstractDataFrame,
                        deltas_df::AbstractDataFrame;
                        threshold_range::Union{Nothing, NTuple{2, <:Number}}=nothing,
                        yaxis_metric::Symbol = :maxcomponent_size)
    traces = Vector{GenericTrace}()
    permstats_filtered_df = sort!(isnothing(threshold_range) ? permstats_df :
                                  filter(r -> (threshold_range[1] <= coalesce(r.threshold_binmid, NaN)) &&
                                              (coalesce(r.threshold_binmid, NaN) <= threshold_range[2]),
                                         permstats_df), :threshold_bin)
    realstats_filtered_df = sort!(isnothing(threshold_range) ? realstats_df :
                                  filter(r -> threshold_range[1] <= coalesce(r.threshold, NaN) <= threshold_range[2],
                                         realstats_df), :threshold)
    isempty(realstats_filtered_df) && return nothing
    realstats_extrema = extrema(realstats_filtered_df[!, yaxis_metric])
    permstats_extrema = (minimum(permstats_filtered_df[!, Symbol(yaxis_metric, "_25")]),
                         maximum(permstats_filtered_df[!, Symbol(yaxis_metric, "_75")]))
    metric_range = (0.5*min(realstats_extrema[1], permstats_extrema[1]),
                    2.0*max(realstats_extrema[2], permstats_extrema[2]))
    permstats_filtered_df[!, Symbol(yaxis_metric, "_025_clamped")] = max.(permstats_filtered_df[!, Symbol(yaxis_metric, "_025")], metric_range[1])
    permstats_filtered_df[!, Symbol(yaxis_metric, "_975_clamped")] = min.(permstats_filtered_df[!, Symbol(yaxis_metric, "_975")], metric_range[2])

    append_CI_band_traces!(
        traces, "permuted", permstats_filtered_df,
        :threshold_binmid, Symbol(yaxis_metric, "_50"),
        (Symbol(yaxis_metric, "_25"), Symbol(yaxis_metric, "_75")),
        true, true) do trace, col
        trace[:legendgroup] = "permuted"
        trace[:showlegend] = col == Symbol(yaxis_metric, "_50")
        trace[:marker_color] = "gray"
        return trace
    end
    PlotlyUtils.append_CI_band_traces!(
        traces, "permuted", permstats_filtered_df,
        :threshold_binmid, nothing,
        (Symbol(yaxis_metric, "_025_clamped"), Symbol(yaxis_metric, "_975_clamped")),
        false, false) do trace, col
        trace[:legendgroup] = "permuted"
        trace[:showlegend] = false
        trace[:marker_color] = "gray"
        trace[:line_dash] = "dot"
        return trace
    end

    push!(traces,
        scatter(realstats_filtered_df,
                x=:threshold, y=yaxis_metric, name="real",
                marker_color="firebrick", legendgroup="real"))

    if deltas_df !== nothing
        for delta_df in groupby(deltas_df, [:stat, :value_type])
            tr_df = DataFrame(source = ["real", "permuted"],
                              stat = delta_df.stat[1],
                              value_type = delta_df.value_type[1],
                              threshold = delta_df.threshold[1],
                              metric = [delta_df.value[1], delta_df.permuted_value[1]],
                    )
            tr_df[!, :hover] .= @sprintf("<b>%s</b><br>threshold = %.5f<br>%s_%s(Δ) = %.4f<br>real = %.4f<br>permuted = %.4f",
                                         yaxis_metric, tr_df.threshold[1],
                                         delta_df.stat[1], delta_df.value_type[1],
                                         delta_df.delta[1],
                                         coalesce(tr_df.metric[1], NaN), coalesce(tr_df.metric[2], NaN))
            if any(ismissing, tr_df.metric)
                @warn "No $yaxis_metric found for $(delta_df.stat[1])(Δ)"
            elseif !isnothing(threshold_range) &&
                   !(threshold_range[1] <= tr_df.threshold[1] <= threshold_range[2])
                @warn "$yaxis_metric $(delta_df.stat[1])(Δ) threshold=$(tr_df.threshold[1]) is outside of displayed range"
            else
                push!(traces,
                      scatter(tr_df, x=:threshold, y=:metric, hovertext=:hover,
                              marker_size=2, marker_color=delta_df.stat[1] == "min" ? "blue" : "red",
                              showlegend=false, name=string(delta_df.stat[1], "(Δ)")))
            end
        end
    end

    plot(traces)
end
