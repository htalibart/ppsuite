import os
from Bio import SeqIO
from Bio import ExPASy


class SwissprotAccessor:
    """
    Accessor to Swissprot records using the website (default) or local .dat files
    Author : F. Coste
    """
    def __init__(self, sprot_filepath = None):
        """
        @param sprot_filepath:  path of swissprot flat file. 
                                If None, accessor will send requests to the web server 
        """
        self.sprot_filepath = sprot_filepath
        if sprot_filepath is None:
            # use web server, no need of index
            self.index = None
        else:
            index_filepath = os.path.splitext(sprot_filepath)[0] + ".idx"
            if (not os.path.isfile(index_filepath)):
                # need to create index file
                if (not os.path.isfile(self.sprot_filepath)):
                    errmsg = f"""
                    Swissprot data file missing: {sprot_filepath} 
                    Download and uncompress file uniprot_sprot.dat.gz 
                    from ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/
                    """
                    raise OSError(errmsg)
                print(f"Indexing {self.sprot_filepath}...")
                self.index = SeqIO.index_db(index_filepath, sprot_filepath, "swiss")
            else:
                # use existing index file
                self.index = SeqIO.index_db(index_filepath)
                
        
    def get_record(self, accession_id):
        if self.index:
            return self.index[accession_id]
        else:
            with ExPASy.get_sprot_raw(accession_id) as sp_handle:
                return SeqIO.read(sp_handle, "swiss")
        
    def iter_records(self, accession_ids):
        for acc in accession_ids:
            try:
                yield self.get_record(acc)
            except:
                print(f"Warning: Failed to retrieve {acc}! Are you sure of this accession id ???")
