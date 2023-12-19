using JLD2, CairoMakie, Images, ImageMagick

## Functions

rast = 5

function find_border_num(data, xs, ys)
    mesh = collect(Iterators.product(xs, ys))

    edges = canny(data, (0.1, 0.2))
    border = findall( x -> x == 1, edges)

    return border, mesh[border]
end


function plot_PD(ax, dict; colorrange = (2, 4), dlim = 0.001, colormap = :haline, points = nothing,  kw...)

    xs, ys = dict["xs"], dict["ys"] 
    data = dict["PD"] 
    data_masked = getindex.(data, 1) .* (getindex.(data, 2) .> dlim)
    replace!(x -> x == 0 ? NaN : x, data_masked)

    cmap = cgrad(colormap, )

    data_plot = log10.(data_masked .* 5)


    hmap = heatmap!(ax, xs, ys, data_plot; colormap = cmap, colorrange,  rasterize = rast, kw...)

    if points !== nothing
        for point in points
            scatter!(ax, point[1]...; color = point[2])
        end
    end

    return hmap
end

function cuty(dict, ylims; isPf = false)
    s = size(dict["PD"])    
    xs = get(dict, "xs", zeros(s[1]))
    ys = get(dict, "ys", zeros(s[1]))
    PD = get(dict, "PD", zeros((s...)))
    PD0 = get(dict, "PD0", zeros((s...)))

    y1 = findmin(abs.(ys .- ylims[1]))[2]
    y2 = findmin(abs.(ys .- ylims[2]))[2]


    if isPf
        sy = size(ys)
        y1PD = (y1 * s[1] / sy[1] |> ceil |> Int) - 2
        y2PD = (y2 * s[2] / sy[1] |> ceil |> Int) + 2
        y1 -= 1
        y2 += 1
    else 
        y1PD = y1
        y2PD = y2
    end

    return Dict("xs" => xs, "ys" => ys[y1:y2], "PD" => PD[:, y1PD:y2PD], "PD0" => PD0[:, y1PD:y2PD])
end

## Figure

function Fig_SCM_PD(files; colorrange = (2, 4), αlims = (-200, 200))
    fig = Figure(resolution = (0.9 * 600, 4/3 * 650), font = "CMU Serif Roman")
    
    ylabel = L"$ \langle \alpha \rangle$ (meVnm)"
    ticks = ([colorrange[1], colorrange[2]], [L"10^%$(colorrange[1])", L"10^%$(colorrange[2])"])

    # Dome v α
    dictD = load(files["dome"])
    dict = cuty(dictD, αlims)

    xlabel = L"$U_{min}$ (meV)"
    xticks = [-1, -50, -100, -150]
    yticks = [ -200, -100, 0, 100, 200]
    ax = CairoMakie.Axis(fig[1, 1]; xlabel, ylabel, xticks, yticks, xreversed = true)

    hmap = plot_PD(ax, dict; colorrange)
    vlines!(ax, -30; color = :navyblue, linestyle = :dash)
    text!(ax, -100, -100; text = L"\mu  = 0", align = (:center, :center))

    hidexdecorations!(ax, ticks = false, grid = false)

   # Dome islands
   p = files["dome_islands"] 
   f = [[p * string(i) * "_" * string(j) * ".jld2" for j in 1:5] for i in 1:5] 

    rows = []
    Vs = []
    αs = []
    for file_row in f
        row_blocks = []
        for file in file_row
            dictf = load(file)
            push!(row_blocks, dictf["PD"])
            push!(Vs, dictf["xs"])
            push!(αs, dictf["ys"])
        end
        push!(rows, hcat(row_blocks...))
    end
    Vs = vcat(Vs...) |> unique!
    αs = vcat(αs...) |> unique!
    pd = vcat(rows...)

    dictb = cuty(Dict("xs" => Vs, "ys" => αs, "PD" => pd), αlims)

    border, μα_border = find_border_num(dictb["PD"], dictb["xs"], dictb["ys"])

    scatter!(ax, getindex.(μα_border, 1), getindex.(μα_border, 2); color = :red, markersize = 4)


    Colorbar(fig[1, 2] , hmap, label = L"$\xi_M$ (nm)", ticklabelsvisible = true,
    labelpadding = -20,  width = 10, ticksize = 2,
    ticklabelpad = 5, ticks = ticks, labelcolor = :black,)
    
    # Dome Pfaffian 
    dict_P = load(files["dome_Pfaffian"])
    ax = CairoMakie.Axis(fig[2, 1]; xlabel, ylabel, xticks, yticks, xreversed = true)

    PD = -dict_P["PD0"]

    dict = cuty(Dict("xs" => dictD["xs"], "ys" => dictD["ys"], "PD" => PD), αlims; isPf = true)
    #replace!(x -> x == -1 ? NaN : x, PD)


    hmap = heatmap!(ax, reverse(dict["xs"]), dict["ys"], dict["PD"]; colormap = [(:white, 0), :red], rasterize = rast)
    
    text!(ax, -83, 100; text = L"$Q = -1$", color = :white, align = (:center, :center), fontsize = 12)
    text!(ax, -60, 100; text = L"$Q = 1$", color = :black, align = (:center, :center), fontsize = 12)

    
    Colorbar(fig[2, 2]; colormap = cgrad([:white, :red,], categorical = true), ticks = ([0.25, 0.75], ["0", "1"]),   width = 10, ticksize = 2, ticklabelpad = 5, label = L"$N_M$", labelpadding = -7)
    
    # Rigid μ v α 
    dictR = load(files["rigid"])
    dict = cuty(dictR, αlims)
    xlabel = L"$\mu$ (meV)"
    ax = CairoMakie.Axis(fig[3, 1]; xlabel, ylabel, yticks)

    hmap = plot_PD(ax, dict; colorrange,)
    vlines!(ax, 0; color = :navyblue, linestyle = :dash)
    text!(ax, -15.5, 150; text = L"$U_{min} = -30$meV$", align = (:center, :center), fontsize = 12)

    hidexdecorations!(ax, ticks = false, grid = false)

    # Rigid islands
    dict = load(files["rigid_islands"])

    αs = dict["ys"]
    μs = dict["xs"]
    pd = dict["PD"]

    dictb = cuty(Dict("xs" => μs, "ys" => αs, "PD" => pd), αlims)

    border, μα_border = find_border_num(dictb["PD"], dictb["xs"], dictb["ys"])



    scatter!(ax, getindex.(μα_border, 1), getindex.(μα_border, 2); color = :red, markersize = 4)

    scatter!(ax, -11, 20; color = :orange, markersize = 8)  
    scatter!(ax, 2, 20; color = :white, markersize = 8)

    scatter!(ax, -8.5, 200; color = :purple)

    Colorbar(fig[3, 2] , hmap, label = L"$\xi_M$ (nm)", ticklabelsvisible = true,
    labelpadding = -20,  width = 10, ticksize = 2,
    ticklabelpad = 5, ticks = ticks, labelcolor = :black,)
    
    # Rigid Pfaffian 
    dict_P = load(files["rigid_Pfaffian"])

    dict = cuty(Dict("xs" => dictR["xs"], "ys" => dictR["ys"], "PD" => dict_P["PD0"]), αlims; isPf = true)

    ax = CairoMakie.Axis(fig[4, 1]; xlabel, ylabel, yticks)
    hmap = heatmap!(ax, dict["xs"], dict["ys"], -dict["PD"]; colormap = [(:white, 0), :red], rasterize = rast)
    text!(ax, 3, 100; text = L"$Q = -1$", color = :white, align = (:center, :center), fontsize = 12)
    text!(ax, -3, 100; text = L"$Q = 1$", color = :black, align = (:center, :center), fontsize = 12)



    

    Colorbar(fig[4, 2]; colormap = cgrad([:white, :red,], categorical = true), ticks = ([0.25, 0.75], ["0", "1"]),   width = 10, ticksize = 2, ticklabelpad = 5, label = L"$N_M$", labelpadding = -7)

    style = (font = "CMU Serif Bold", fontsize = 20)
    Label(fig[1, 1, TopLeft()], "a";  padding = (-40, 0, -22, 0), style...)

    Label(fig[2, 1, TopLeft()], "b";  padding = (-40, 0, -22, 0), style...)

    Label(fig[3, 1, TopLeft()], "c";  padding = (-40, 0, -22, 0), style...)

    Label(fig[4, 1, TopLeft()], "d";  padding = (-40, 0, -22, 0), style...)

    colgap!(fig.layout, 1, 5) 
    rowgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 3, 5)
    return fig
end

## Generate Figure

p = "Data/Fig_SCM_PD/"
try run(`mkdir Figures`) catch end

files = Dict("dome" => "$(p)SCM_dome_α.jld2", 
"dome_Pfaffian" => "$(p)PD_SCM_dome_n=1.jld2",
"rigid" => "$(p)SCM_rigid_μα.jld2",
"rigid_Pfaffian" => "$(p)PD_SCM_rigid_n=1.jld2",
"dome_islands" => "$(p)SCM_dome/SCM_dome_μα_",
"rigid_islands" => "$(p)SCM_dome/SCM_rigid.jld2")

fig = Fig_SCM_PD(files)
save("Figures/FigSCM_PD.pdf", fig)
fig
