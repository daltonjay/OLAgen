import pandas as pd

class GlobalState:
    def __init__(self):
        self._global_muts = None
        self._target_names = None
        self._global_AA_seqs = None
        self._global_SeqIO_seqs = None
        self._global_storage_df = pd.DataFrame(columns=['Target', 'SNP', 'VP_Probe', 'WT_Probe', 'CP_Probe'])
        self._primer_row_counter = -1

    @property
    def global_muts(self):
        return self._global_muts

    @global_muts.setter
    def global_muts(self, value):
        self._global_muts = value

    @property
    def target_names(self):
        return self._target_names

    @target_names.setter
    def target_names(self, value):
        self._target_names = value

    @property
    def global_AA_seqs(self):
        return self._global_AA_seqs

    @global_AA_seqs.setter
    def global_AA_seqs(self, value):
        self._global_AA_seqs = value

    @property
    def global_SeqIO_seqs(self):
        return self._global_SeqIO_seqs

    @global_SeqIO_seqs.setter
    def global_SeqIO_seqs(self, value):
        self._global_SeqIO_seqs = value

    @property
    def global_storage_df(self):
        return self._global_storage_df

    @global_storage_df.setter
    def global_storage_df(self, value):
        self._global_storage_df = value

    @property
    def primer_row_counter(self):
        return self._primer_row_counter

    @primer_row_counter.setter
    def primer_row_counter(self, value):
        self._primer_row_counter = value
