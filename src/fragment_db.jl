using BioSequences
import JSON

type Fragment 
    seq::DNASequence
    description::String
end


type FragDB
    l::Int
    max_rep_target_dist::Int
    reps_to_fragments::Dict{DNASequence, Array{Fragment}}
end

function create_fragment_db(l::Int, max_rep_target_dist::Int) :: FragDB
    db = FragDB(l, 
                max_rep_target_dist, 
                Dict{DNAKmer{k}, Array{DNASequence}}(), 
                Dict{DNASequence, Array{DNASequence}}())
    return db
end

function hamming_dist(a::DNASequence, b::DNASequence, thresh::Int) :: Bool
    d = 0
    for i in 1:length(a)
        d += a[i] != b[i] ? 1 : 0 
        if d > thresh
            return false
        end
    end
    return d <= thresh
end

function add_fragment_to_db!(db::FragDB, frag::DNASequence, fraginfo::String)
    fragObj = Fragment(frag, fraginfo)
    matched = false
    for rep in keys(db.reps_to_fragments)
        if hamming_dist(frag, rep, db.max_rep_target_dist)
            push!(db.reps_to_fragments[rep], fragObj)
            matched = true
            break
        end
    end

    if !matched
        db.reps_to_fragments[frag] = [fragObj]
    end
end

function write_db_to_file(db::FragDB, fname::String)
    j = Dict{String, Any}[]
    for rep in keys(db.reps_to_fragments)
        seq = convert(String, rep)
        n = length(db.reps_to_fragments[rep])
        obj = Dict("rep" => seq,
                   "seqs" => String[],
                   "descriptions" => String[])
        for frag in db.reps_to_fragments[rep]
            push!(obj["descriptions"], frag.description)
            push!(obj["seqs"], convert(String, frag.seq))
        end
        push!(j, obj)
    end
    obj = Dict("l" => db.l,
               "max_rep_target_dist" => db.max_rep_target_dist,
               "clusters" => j)  
    iof = open(fname, "w")
    JSON.print(iof, obj)
    close(iof)
end

function read_db_from_file(db::FragDB, fname::String)

end