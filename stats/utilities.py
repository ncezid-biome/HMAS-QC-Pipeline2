def revcomp(myseq):
    rc = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G', 'U' : 'A', 'Y' : 'R', 'R' : 'Y', 'K':'M', 'M':'K','B':'V',\
            'D':'H', 'H':'D', 'V':'B', 'N':'N'}
    seq = [rc[n] if n in rc else n for n in myseq] # allow non-IUPAC code to stay as is
    return("".join(list(reversed(seq))))


class Primers:
    'Parses mothur oligos file and extracts the primer names. Stole from A. Jo Williams-Newkirk'
    
    def __init__(self, fname):
        self.fname = fname
        self.pseqs = dict()
        Primers.reader(self, self.fname, self.pseqs)
        self.pnames = self.pseqs.keys()
    
    def reader(self, fname, pseqs):
        with open(fname, 'r') as infile:
            for line in infile:
                if line.startswith("primer"):
                    tmp = line.split('\t')
                    pseqs[tmp[3].strip('\n')] = [tmp[1], revcomp(tmp[2])]


def create_fasta_dict(fasta):
    '''
    this method reads a fasta file and convert it into a dictionary, with seq_ID being the key and actual sequence
    being the value

    Note: if there are duplicate seq ID, only the first seq_ID/sequence will be saved. 

    Parameters
    ----------
    fasta: String name of the fasta file

    Returns a dictionary

    '''
    seq_dict = {}
    with open(fasta, 'r') as f:
        for ind, row in enumerate(f.readlines(), start=1):
            if ind%2 == 1:
                last_seqID = row.strip().split()[0][1:] #removes the '>'
            else:
                last_seq = row.strip()
                if last_seqID in seq_dict:
                    print (f"warning ! has a duplicate seq ID {last_seqID}")
                else:
                    seq_dict[last_seqID] = last_seq
    return seq_dict
