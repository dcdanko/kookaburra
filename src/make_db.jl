using ArgParse


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


@task function make_subseqs(record, k)
    seq = FASTA.sequence(record)
    for i in 1:len(seq) - k + 1
        subseq = seq[i:i+k]
        produce(subseq)
    end

end

function parse_fasta(input_fasta, db::FragDB)
    reader = FASTA.Reader(open(input_fasta, "r"))
    for record in reader
        seq_info = FASTA.identifier(record) + " " + FASTA.description(record)
        for subseq in make_subseqs(record, db.l)
            add_fragment_to_db!(db, subseq, seq_info)
        end
    end
    close(reader)
end

begin
    args = parse_cl()
    db = create_fragment_db(args.l, args.d)
    for fasta in args.fastas
        parse_fasta(fasta, db)
    end
    write_db_to_file(db, args.out)
end
