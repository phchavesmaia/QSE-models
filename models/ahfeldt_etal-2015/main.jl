# *********************
# **** Load Files  **** 
# *********************
using  Plots, LaTeXStrings, FixedEffectModels, GeoStats, GeoIO, CSV, DataFrames, Statistics, LinearAlgebra
import CairoMakie as Mke

function load_dir(dir::String)
    files = readdir(dir)
    for file in files
        if endswith(file, ".jl")
            include(joinpath(dir, file))
        end
    end
end

try 
    cd("D:/Dropbox/learn-julia/qse/models/ahfeldt_etal-2015/")
catch
    cd("C:/Users/pedro.maia/Dropbox/learn-julia/qse/models/ahfeldt_etal-2015/")
end
load_dir("functions")