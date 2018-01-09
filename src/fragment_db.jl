
type Fragment 
    seq::DNASequence
    description::String
end


type FragDB
    k::Int
    l::Int
    max_rep_target_dist::Int
    kmers_to_reps::Dict{DNAKmer, Array{DNASequence}}
    reps_to_fragments::Dict{DNASequence, Array{Fragment}}
end

function create_fragment_db(l::Int, max_rep_target_dist::Int)
    k = floor((l - max_rep_target_dist) / (max_rep_target_dist + 1))
    db = FragDB(k, 
                l, 
                max_rep_target_dist, 
                Dict{DNAKmer{k}, Array{DNASequence}}(), 
                Dict{DNASequence, Array{DNASequence}}())
    return db
end

function hamming_dist(a::DNASequence, b::DNASequence)
    d = 0
    for i in len(a)
        d += a[i] != b[i] ? 1 : 0 
    end
    return d
end

function add_fragment_to_db!(db::FragDB, frag::DNASequence, fraginfo::String)
    
    matched = false
    kmers = map(x -> x[2], each(DNAKmer{k}, seq))
    for kmer in kmers
        for rep in get(db.kmers_to_reps, kmer, [])
            d = hamming_dist(frag, rep)
            if d <= db.max_rep_target_dist
                push!(db.reps_to_fragments[rep], frag)
                matched = true
                break
            end
        end
        if matched
            break
        end
    end

    if !matched
        for kmer in kmers
            try
                push!(db.kmers_to_reps[kmer], frag)
            catch error
                if isa(error, KeyError)
                    db.kmers_to_reps[kmer] = [frag]
                end
            end
        end
        db.reps_to_fragments[frag] = Fragment(frag, fraginfo)
    end
end

function write_db_to_file(db::FragDB, fname::String)
    j = Array{Dict{String, Array{String}}}()
    for rep in keys(db.reps_to_fragments)
        seq = convert(String, rep)
        obj = Dict('rep' => rep,
                   'seqs' => [],
                   'descriptions' => [])
        for frag in db.reps_to_fragments[rep]
            push!(obj, frag.description)
            push!(obj, convert(String, frag.seq))
        end
        push!(j, obj)
    end
    obj = Dict('k' => db.k,
               'l' => db.l,
               'max_rep_target_dist' => db.max_rep_target_dist,
               'clusters' => j)  
    iof = open(fname, 'w')
    json.print(iof, obj)
    close(iof)
end

function read_db_from_file(db::FragDB, fname::String)

end