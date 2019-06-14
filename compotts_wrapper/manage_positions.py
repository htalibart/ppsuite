def get_real_aln_df(res_aln_file, compotts_objects):
    df_res = pd.read_csv(res_aln_file)
    c_names = list(df_res.columns)
    df_dict = {}
    for k in range(2):
        df_dict[c_names[k]] = compotts_objects[k].get_real_positions(df_res[c_names[k]].tolist())
    df = pd.DataFrame(df_dict)
    return df


def get_seqs_aligned(res_aln_file, compotts_objects):
    df = get_real_aln_df(res_aln_file, compotts_objects)
    positions = [df['pos_ref'].tolist(), df['pos_2'].tolist()]
    seqs_aligned = ["",""]

    for k in range(2):
        seqs_aligned[k]+='-'*(max(positions[0][0], positions[1][0])-1)+compotts_objects[k].real_seq[positions[k][0]]

    old = [positions[0][0], positions[1][0]]
    for pos in range(1,len(positions[0])):
        gap_length = max([positions[k][pos]-old[k] for k in range(2)])-1
        for k in range(2):
            seqs_aligned[k]+='-'*gap_length+compotts_objects[k].real_seq[positions[k][pos]]
        old = [positions[0][pos],positions[1][pos]]
    gap_length = max([len(compotts_objects[k].real_seq)-old[k] for k in range(2)])
    for k in range(2):
        seqs_aligned[k]+='-'*gap_length
    return seqs_aligned



def get_seq_trimmed_for_ref(res_aln_file, ref_object, query_object):
    df = pd.read_csv(res_aln_file)
    aligned_positions = df['pos_ref'].tolist()
    trimmed_seq=""
    for pos in range(ref_object.mrf.ncol):
        if pos in aligned_positions:
            ind = aligned_positions.index(pos)
            pos_2 = df['pos_2'].tolist()[ind]
            trimmed_seq+=query_object.get_real_letter_at_trimmed_pos(pos_2)
        else:
            trimmed_seq+='-'
    return trimmed_seq


def get_msas_aligned(res_aln_file, train_msa_files, output_msa):
    print("merging "+train_msa_files[0]+" and "+train_msa_files[1])
    df = pd.read_csv(res_aln_file)
    aligned_positions = [df['pos_ref'].tolist(), df['pos_2'].tolist()]
    aligns = [SeqIO.parse(f, "fasta") for f in train_msa_files]
    records = []
    for k in range(2):
        align = aligns[k]
        for record in align:
            new_record = record
            new_seq=""
            for pos in aligned_positions[k]:
               new_seq+=record[pos]
            new_record.seq=Seq(new_seq)
            records.append(new_record)
    new_alignment = MultipleSeqAlignment(records)
    with open(output_msa, 'w') as f:
        AlignIO.write(new_alignment, f, "fasta")


def get_pos_aligned_at_pos(res_aln_file, pos):
    df = pd.read_csv(res_aln_file)
    return int(df[df.pos_ref==pos]['pos_2'])



def get_aligned_v_scores(res_aln_file, v_scores):
    aligned_v_scores = np.zeros(len(pd.read_csv(res_aln_file)['pos_ref'].tolist()))
    pos=0
    for i,k in zip(pd.read_csv(res_aln_file)['pos_ref'].tolist(), pd.read_csv(res_aln_file)['pos_2'].tolist()):
        aligned_v_scores[pos] = v_scores[i][k]
        pos+=1
    return aligned_v_scores
