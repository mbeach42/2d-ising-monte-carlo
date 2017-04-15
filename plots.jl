using PyPlot

function all_plots(T_range, E, M, Cv, X)
    fig = figure("all_plot",figsize=(10,4))
    subplot(221)
    xlabel("Temperature")
    ylabel("Energy")
    plot(T_range, E, ".")
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

function plot_E(T_range, E)
    figure("E_plot", figsize=(10,4))
    xlabel("Temperature")
    ylabel("Energy")
    plot(T_range, E, ".", color = "green");
end

function plot_M(T_range, M)
    figure("M_plot", figsize=(10,4))
    xlabel("Temperature")
    ylabel("Magnetization")
    plot(T_range, abs(M), ".", color = "green");
end

function plot_Cv(T_range, X)
    figure("Cv_plot", figsize=(10,4))
    xlabel("Temperature")
    ylabel(L"$C_v$")
    plot(T_range, X, ".", color = "green");
end

function plot_X(T_range, X)
    figure("X_plot", figsize=(10,4))
    xlabel("Temperature")
    ylabel(L"$\chi$")
    plot(T_range, X, ".", color = "orange");
end
