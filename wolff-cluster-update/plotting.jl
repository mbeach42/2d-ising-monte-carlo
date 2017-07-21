module IsingPlots

using PyPlot

export allplots, plotE, plotM, plotCv, plotX

function allplots(T_range, L, E, M, Cv, X)
    fig = figure("all_plot",figsize=(10,4))
    subplot(221)
    xlabel("Temperature")
    ylabel("Energy")
    plot(T_range, E, ".", label = "L = 1")
    legend()
    subplot(222)
    xlabel("Temperature")
    ylabel("Magnetization")
    plot(T_range, abs(M), ".", color = "red")
    subplot(223)
    xlabel("Temperature")
    ylabel("Cv")
    plot(T_range, Cv, "." , color = "orange")
    subplot(224)
    xlabel("Temperature")
    ylabel("Sus")
    plot(T_range, X, ".", color = "green");
end

end
