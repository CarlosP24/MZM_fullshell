using JLD2, CSV, CairoMakie, ColorSchemes, Colors, Makie.GeometryBasics

## Functions 
function schematics(ax)
    poly!(ax, Circle(Point2f(0,0), 10f0), color = colorant"#B9C8F7", strokewidth = 2)
    poly!(ax, Circle(Point2f(0,0), 8f0), color = colorant"#FBF5CD", strokewidth = 2)
    poly!(ax, Circle(Point2f(0,0), 5f0), color = :white, strokewidth = 2)
    poly!(ax, Circle(Point2f(-8, 8), 0.5f0), color = :white, strokewidth = 2)
    poly!(ax, Circle(Point2f(-8, 8), 0.1f0), color = :black, )
    arrows!([0, -5.3, -7.7, -8.3, -9.7], [0, 0, 0, 0, 0], [7.5, -2.2, 2.2, -1.2, 1.2], [0, 0, 0, 0, 0], color = :black, linewidth = 2, arrowsize = 5
    )
    
    text!(ax, .8, 0.5; text = L"R")
    text!(ax, -7.4, 0.5; text = L"w") 
    text!(ax, -9.8, 0.5; text = L"d") 
    text!(ax, -10, 7.2; text = L"z") 

    hidespines!(ax)
    hidedecorations!(ax)
    return
end


function plot_LDOS(ax, dict; colorrange = (1e-5, 1e-2))
    data = sum(values(dict["mjdict"]))
    xs, ys = dict["xs"], real.(dict["ys"])
    return  heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = colorrange, rasterize = 5, lowclip = :black)
end

## Figure

function FigTubular(files, zranges, ws)
    fig = Figure(resolution = (1100, 2/3 * 650), font = "CMU Serif Roman")

    xlabel = L"\Phi/\Phi_0"
    ylabel = L"$\omega$ (meV)"
    label = L"$$LDOS (arb. units)"
    yticks = [-0.2, 0, 0.2]
    xlabelpadding = -7

    for (row, rowfiles) in enumerate(files), (col, file) in enumerate(rowfiles)
        (row == 1) & (col == 1) ? continue : 1
            
        ax = CairoMakie.Axis(fig[row, col]; xlabel, ylabel, yticks, xlabelpadding)
        dict = load(file)

        plot_LDOS(ax, dict; colorrange = zranges[row][col])

        w = ws[row][col]
        text = w == 0 ? L"$w \rightarrow 0$" : L"$w = %$(w)$ nm"
        text!(ax, 0.5, 0.2; text = text, color = :white, align = (:center, :center))
        
        if (row == 2) & (col == 1)
            text!(ax, 1.25, -0.26; text = "Skewed \nCdGM \nanalogs", align = (:center, :bottom), color = :white, fontsize = 14)
            arrows!(ax, [1], [-0.2], [-0.1], [0.04], color = :white)
            arrows!(ax, [1], [-0.15], [-0.1], [0.07], color = :white)
            text!(ax, 1.3, 0; text = "Shifted \ngap", align = (:center, :center), color = :white, fontsize = 14)
        end
        
        col == 3 && Colorbar(fig[row, 4], colormap = :thermal, limits = (0, 1), label = label, ticklabelsvisible = true, ticks = [0,1], labelpadding = -10,  width = 10, ticksize = 2, ticklabelpad = 5)

        row != 2 && hidexdecorations!(ax, ticks = false) 

        row != 1 && col != 1 &&  hideydecorations!(ax, ticks = false)
        row == 1 && col != 2 && hideydecorations!(ax, ticks = false)
    
    end

    ax = CairoMakie.Axis(fig[1, 1]; aspect = DataAspect())
    schematics(ax)

    style = (; font = "CMU Serif Bold", fontsize = 20)
    Label(fig[1, 1, TopLeft()], "a"; padding = (-40, 0, -15, 0), style...)
    Label(fig[1, 2, TopLeft()], "b"; padding = (35, 0, -15, 0), style...)
    Label(fig[1, 3, TopLeft()], "c"; padding = (-20, 0, -15, 0), style...)
    Label(fig[2, 1, TopLeft()], "d"; padding = (-40, 0, -20, 0), style...)
    Label(fig[2, 2, TopLeft()], "e"; padding = (35, 0, -20, 0), style...)
    Label(fig[2, 3, TopLeft()], "f"; padding = (-20, 0, -20, 0), style...)


    colgap!(fig.layout, 1, -35)
    colgap!(fig.layout, 2, 15)
    colgap!(fig.layout, 3, 5)
    rowgap!(fig.layout, 1, 5)

    return fig

end

## Generate Figure

p = "Data/Fig_TCM_LDOS/"
try run(`mkdir Figures`) catch end

files = [
    [ 0, 68, 60],
    [50, 40, 30],
]

inpath(file) = "$(p)LDOS_$(file).jld2"

files = map(row -> inpath.(row),files)

zranges = [
    [(1e-5, 1e-2), (1.17e-3, 4e-2), (1e-3, 4e-2)],
    [(7e-4, 3e-2), (6e-4, 2e-2), (4e-4, 2e-2)],
]

ws = [
    [ 0,  0, 10],
    [20, 30, 40],
]

fig = FigTubular(files, zranges, ws)
save("Figures/FigTCM_LDOS.pdf", fig)
fig

