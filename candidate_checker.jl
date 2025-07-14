
# %%
import Pkg
using Pkg
Pkg.activate("test_env")
Pkg.add(PackageSpec(name="NautyGraphs", version="0.4.0"))
Pkg.add(["GraphIO", "Graphs", "SQLite", "DataFrames", "DBInterface", "Combinatorics", "TimerOutputs", "BenchmarkTools", "Glob","Plots", "ThreadsX", "JET"])

using SparseArrays, LinearAlgebra
using NautyGraphs, Graphs
using NautyGraphs: WordType
using NautyGraphs: DenseNautyGraph
using NautyGraphs: AbstractNautyGraph
using Graphs, TimerOutputs, BenchmarkTools, Glob
using GraphIO.Graph6
using NautyGraphs
using Combinatorics
using Base.Threads
using SQLite
using Serialization

function asymmetric_depth_checkpoint(g::NautyGraph, heuristics::Bool, verbose::Bool=false, max_depth::Int=typemax(Int), start_depth::Int=1, checkpoint_file::String="checkpoint.jls", result_file::String="final_result.jls")
    # Attempt to load checkpoint state if it exists.
    if isfile(checkpoint_file)
        # state = deserialize(checkpoint_file)
        i = state["i"]
        canonical_table = state["canonical_table"]
        n_g = state["n_g"]
        if verbose
            println("Resuming from checkpoint: i = ", i)
        end
    else
        i = start_depth
        canonical_table = Set{UInt64}()
        n_g = 0
    end

    n = g.n_vertices
    grpsize::Int64, _, _ = NautyGraphs.nauty(g, false)
    if grpsize > 1
        # Save result before returning.
        # serialize(result_file, 0)
        return 0
    end

    @inbounds for i in i:n-1
        if verbose
            println("i: ", i)
        end
        
        it = 0
        @inbounds for comb in combinations(1:n, n - i)
            it += 1
            if it < n_g
                continue
            end 

            subg, _ = Graphs.induced_subgraph(g, comb)
            grpsize, _, _, canon_hash_sub::UInt64 = NautyGraphs._nautyhash(subg)
            if grpsize > 1 || (canon_hash_sub in canonical_table)
                # Finished computation; save result i.
                println("comb: ", comb)
                println("canon_hash_sub: ", canon_hash_sub)
                # serialize(result_file, i)
                # if isfile(checkpoint_file)
                #     rm(checkpoint_file)
                # end
                return i
            end
            push!(canonical_table, canon_hash_sub)
            n_g += 1
        end



        # Save progress after each outer-loop iteration.
        empty!(canonical_table)

        if i == max_depth
            return max_depth
        end

    end

    # Remove the checkpoint if computation finished.
    # if isfile(checkpoint_file)
    #     rm(checkpoint_file)
    # end
    serialize(result_file, n - 1)
    return n - 1
end

# asymmetric conference graphs on 25 vertices
# g6_lines = """X~zfCqTc{YPR`jUQidaeNRLXIrIMphoxKsVXKixPZCnD[fBHuQl"""
# g6_lines =  """X~zfCqTc{YPT`fUQidaeNRKxItIMpholosZFKjXHZGnDZDYHwuF"""
g6_lines =  """~?@_s`aaOoI?L?iE??o_EO?ZA??@?G?_O?GK?AG??PG?@D??AJC????_??_C??O?O?B??_?_?@?@G?@?@O??_?k???A?????OG???@@????AB????AG????@C_????PO????AJG??????C?????G@?????O?_????O?G????E?@????@O?@????E??A????F??????A???????@?_??????OO??????AC???????GW???????PO???????P_???????G[???????A?????????OG???????@@????????AC????????AE????????@D?????????P_????????AFO??????????C?????????GC?????????OA?????????O?G????????E?@????????@O?A????????E??A????????F??????????A???????????@?_??????????OO??????????AC???????????GW???????????PO???????????P_???????????G[???????????A?????????????OG???????????@@????????????AC????????????AE????????????@D?????????????P_????????????AF?????????????G??????????????OG?????????????OO?????????????GO?????????????AE??????????????PO?????????????@E??????????????AF"""

# using Serialization
# using ThreadsX

# # Files to store progress and collected results.
# progress_file = "progress.jls"
# results_file = "all_results.jls"

# # Load previous progress if available (last processed index)
# last_idx = isfile(progress_file) ? deserialize(progress_file) : 0

# # Load previous results if available
# results = isfile(results_file) ? deserialize(results_file) : Dict{Int, Int}()

# g6_lines = split(g6_lines, "\n")

# for (idx, g6_str) in enumerate(g6_lines)
#     println("Processing graph ", idx, ": ", g6_str)
#     if isempty(strip(g6_str))
#         continue
#     end
#     # Skip already processed graphs based on progress.
#     if idx <= last_idx
#         println("Skipping graph ", idx, " (already processed)")
#         continue
#     end

#     println("Processing graph ", idx, ": ", g6_str)
#     g = Graph6._g6StringToGraph(g6_str)
#     A = Graphs.adjacency_matrix(g)
#     ng = NautyGraph(A)
    
#     # Use unique checkpoint and result files per graph.
#     ckpt_file = "experiments/candidate_checkpoints/checkpoint_$(idx).jls"
#     res_file  = "experiments/candidate_checkpoints/final_result_$(idx).jls"

#     print(typeof(ckpt_file))

#     res = asymmetric_depth__checkpoint(ng, false, true, 19, 17, ckpt_file, res_file)
    
#     results[idx] = res
#     # Update progress file and results file after processing each graph.
#     serialize(progress_file, idx)
#     serialize(results_file, results)
    
#     println("Graph ", idx, " result: ", res)
# end

# println("All results saved in ", results_file)

using ThreadsX
using Serialization
using Base.Threads: ReentrantLock

progress_file_set = "progress_set.jls"
results_file = "all_results.jls"

processed_indices = isfile(progress_file_set) ? deserialize(progress_file_set) : Set{Int}()
results = isfile(results_file) ? deserialize(results_file) : Dict{Int, Int}()

g6_lines = split(g6_lines, "\n")

# Prepare list of indices that haven't been processed and are non-empty.
indices_to_process = [ idx for (idx, line) in enumerate(g6_lines)
    if !(idx in processed_indices)
]

println("Indices to process: ", indices_to_process)

const progress_lock = ReentrantLock()

results_list = ThreadsX.map(indices_to_process) do idx
    println("Processing graph ", idx, ": ", g6_lines[idx])
    
    g = Graph6._g6StringToGraph(g6_lines[idx])
    A = Graphs.adjacency_matrix(g)
    ng = NautyGraph(A)
    
    # Use unique checkpoint and result files per graph.
    ckpt_file = "projects/partial_symmetries_/experiments/2025_experiment_smallest_graph_satisfying_bound_asymmetric_depth/checkpoint_$(idx).jls"
    res_file  = "projects/partial_symmetries_/experiments/2025_experiment_smallest_graph_satisfying_bound_asymmetric_depth/final_result_$(idx).jls"
    
    res = asymmetric_depth__checkpoint(ng, false, true, 6, 1, ckpt_file, res_file)
    
    # Save progress and result immediately.
    # lock(progress_lock) do
    #     processed_indices = isfile(progress_file_set) ? deserialize(progress_file_set) : Set{Int}()
    #     push!(processed_indices, idx)
    #     serialize(progress_file_set, processed_indices)
        
    #     results = isfile(results_file) ? deserialize(results_file) : Dict{Int, Int}()
    #     results[idx] = res
    #     serialize(results_file, results)
    # end
    # (idx, res)
end

println("All results saved in ", results_file)
println("Results: ", results_list)


# sort -n -m hash_batch_*.txt -o sorted_hashes.txt

##
