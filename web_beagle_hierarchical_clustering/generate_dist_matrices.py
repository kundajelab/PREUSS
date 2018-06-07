import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="accepts output from Web Beagle; generates pairwise distance matrix from RNA sequences")
    parser.add_argument("--web_beagle_output")
    parser.add_argument("--o")
    return parser.parse_args()

def main():
    args=parse_args()
    
    #read in the web beagle output file 
    data=open(args.web_beagle_output,'r').read().strip().split('\n')

    #generate dictionaries to keep track of seq_identity, str_identity, str_similarity
    seq_id_dict=dict()
    str_id_dict=dict()
    str_sim_dict=dict() 
    seqs=set([]) 
    for line in data:
        if line.startswith('>'):
            tokens=line.split('|')
            seq1=tokens[0].strip('>')
            seq2=tokens[1]
            seqs.add(seq1)
            seqs.add(seq2) 
            seq_identity=tokens[3].split(':')[1]
            str_identity=tokens[4].split(':')[1]
            str_similarity=tokens[5].split(':')[1]

            #populate seq_id_dict 
            if seq1 not in seq_id_dict:
                seq_id_dict[seq1]=dict()
                seq_id_dict[seq1][seq1]='100'
            if seq2 not in seq_id_dict:
                seq_id_dict[seq2]=dict()
                seq_id_dict[seq2][seq2]='100'
            seq_id_dict[seq1][seq2]=seq_identity
            seq_id_dict[seq2][seq1]=seq_identity

            #populate str_id_dict
            if seq1 not in str_id_dict:
                str_id_dict[seq1]=dict()
                str_id_dict[seq1][seq1]='100'
            if seq2 not in str_id_dict:
                str_id_dict[seq2]=dict()
                str_id_dict[seq2][seq2]='100'
            str_id_dict[seq1][seq2]=str_identity
            str_id_dict[seq2][seq1]=str_identity

            #populate str_sim_dict
            if seq1 not in str_sim_dict:
                str_sim_dict[seq1]=dict()
                str_sim_dict[seq1][seq1]='100'
            if seq2 not in str_sim_dict:
                str_sim_dict[seq2]=dict()
                str_sim_dict[seq2][seq2]='100'
            str_sim_dict[seq1][seq2]=str_similarity
            str_sim_dict[seq2][seq1]=str_similarity
    
    outf_seq_id=open(args.o+".seq_id.txt",'w')
    outf_str_id=open(args.o+".str_id.txt",'w')
    outf_str_sim=open(args.o+'.str_sim.txt','w')
    seqs=list(seqs)
    #write header
    header='\t'+'\t'.join(seqs)
    outf_seq_id.write(header+'\n')
    outf_str_id.write(header+'\n')
    outf_str_sim.write(header+'\n')
    for seqname in seqs:
        outf_seq_id.write(seqname)
        outf_str_id.write(seqname)
        outf_str_sim.write(seqname) 
        for seqname2 in seqs:
            try:
                cur_seq_identity=seq_id_dict[seqname][seqname2]
            except:
                cur_seq_identity='NA'
            try:
                cur_str_identity=str_id_dict[seqname][seqname2]
            except:
                cur_str_identity='NA'
            try:
                cur_str_similarity=str_sim_dict[seqname][seqname2]
            except:
                cur_str_similarity='NA'
            outf_seq_id.write('\t'+cur_seq_identity)
            outf_str_id.write('\t'+cur_str_identity)
            outf_str_sim.write('\t'+cur_str_similarity)
        outf_seq_id.write('\n')
        outf_str_id.write('\n')
        outf_str_sim.write('\n')
        
                       

if __name__=="__main__":
    main()
    
