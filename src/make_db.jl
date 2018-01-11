using ArgParse
using BioSequences
import GZip

include("fragment_db.jl")

function parse_cl()
    s = ArgParseSettings()
    @add_arg_table s begin
        "-l"
            arg_type = Int
        "-d"
            arg_type = Int
        "out"
        "fastas"
            nargs = '*'
        
    end
    return parse_args(ARGS, s)
end


make_subseqs(record::FASTA.Record, k::Int) = Channel(ctype=DNASequence) do c
    seq = FASTA.sequence(record)
    for i in 1:length(seq) - k +1
        subseq = seq[i:i + k - 1]
        push!(c, subseq)
    end
end


function parse_fasta(input_fasta, db::FragDB)
    reader = FASTA.Reader(GZip.open(input_fasta, "r"))
    for record in reader
        seq_info = string(FASTA.identifier(record), " ", input_fasta)
        for subseq in make_subseqs(record, db.l)
            add_fragment_to_db!(db, subseq, seq_info)
        end
    end
    close(reader)
end


begin
    args = parse_cl()
    db = create_fragment_db(args["l"], args["d"])
    for fasta in args["fastas"]
        parse_fasta(fasta, db)
    end
    write_db_to_file(db, args["out"])
end
