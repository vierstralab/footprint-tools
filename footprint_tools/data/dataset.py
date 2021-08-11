from footprint_tools.data.loader import data_loader
from footprint_tools.data.utils import numpy_collate, numpy_collate_concat

from genome_tools import genomic_interval
from footprint_tools import cutcounts   

from tqdm import tqdm

import pandas as pd

import logging
logger = logging.getLogger(__name__)

class base_dataset(object):
    def batch_iter(self, **kwargs):
        raise  NotImplementedError

    def load_all(self, **kwargs):
        """Load an entire dataset"""
        return [x for x in self.batch_iter(**kwargs)]

class dataset(base_dataset):
    """All datasets should subclass this class. All subclases should
    override `__len__` and `__getitem__` to support integer indexing"""

    def __getitem__(self, index):
        """Return one sample"""
        raise NotImplementedError
    
    def __len__(self):
        """Return number of all elements"""
        raise NotImplementedError

    def _batch_iterable(self, batch_size=1, num_workers=0, **kwargs):
        
        dl = data_loader(self, 
                        batch_size=batch_size,
                        collate_fn=numpy_collate,
                        num_workers=num_workers,
                        **kwargs)

        return dl

    def batch_iter(self, batch_size=1, num_workers=0, **kwargs):
        
        dl =  self._batch_iterable(batch_size=batch_size,
                                    num_workers=num_workers,
                                    **kwargs)
        
        return iter(dl)

    def load_all(self, batch_size=1, **kwargs):
        """Load all data"""
        return numpy_collate_concat([x for x in tqdm(self.batch_iter(batch_size, **kwargs))])
