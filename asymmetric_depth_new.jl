
import Pkg
using Pkg
Pkg.activate("asymdepth_env")
# Pkg.add(PackageSpec(name="NautyGraphs", version="0.6.3"))
Pkg.add(["NautyGraphs", "GraphIO", "Graphs", "SQLite", "DataFrames", "DBInterface", "Combinatorics", "Transducers", "TimerOutputs", "BenchmarkTools", "Glob", "Plots", "JET", "LoggingExtras", "FLoops"])

using SparseArrays, LinearAlgebra
using NautyGraphs, Graphs, SHA
using NautyGraphs: DenseNautyGraph
using TimerOutputs, BenchmarkTools, Glob
using GraphIO.Graph6
using Combinatorics
using SQLite, Serialization
using Logging, LoggingExtras
using FLoops

function canonical_hash_groupsize(g::DenseNautyGraph, options::NautyGraphs.NautyOptions=NautyGraphs.default_options(g))
    canong, _, _, statistics = NautyGraphs._densenauty(g, options)
    return statistics.grpsize1 * 10^statistics.grpsize2, hash_sha(canong)
end

function hash_sha(x...)
    io = IOBuffer()
    write(io, (htol(x) for x in x)...)
    return reinterpret(UInt128, sha256(take!(io)))[1]
end

@inline function check_duplicates(arr::AbstractVector{UInt128})
    @inbounds for i in 1:length(arr)-1
        if arr[i] == arr[i + 1]
            return true
        end
    end
    return false
end


function asymmetric_depth(g6_string::String)
    g6 = Graph6._g6StringToGraph(g6_string)
    A = Graphs.adjacency_matrix(g6)
    g = NautyGraph(A)
    i = 1
    canonical_table = UInt128[]
    n = g.n_vertices
    grpsize::Int64, _, _ = NautyGraphs.nauty(g, false)
    if grpsize > 1
        return 0
    end

    @inbounds for i in i:n-1
        empty!(canonical_table)
        @inbounds for comb in combinations(1:n, n - i)
            subg, _ = Graphs.induced_subgraph(g, comb)
            grpsize, canon_hash_sub = canonical_hash_groupsize(subg)
            if grpsize > 1
                println("comb: ", comb)
                println("canon_hash_sub: ", canon_hash_sub)
                return i
            end
            push!(canonical_table, canon_hash_sub)
        end
        sort!(canonical_table)
        if check_duplicates(canonical_table)
            return i
        end
    end
    return n - 1
end
