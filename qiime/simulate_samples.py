
from commands import getstatusoutput
from numpy import arange, array
from numpy.random import shuffle, normal, exponential
from qiime.util import get_tmp_filename
import os

COMMAND_STR = \
"""#!/usr/bin/env Rscript
library('gtools');
prior_vals = %s
total_prior_knowledge = %s
prior_vals = total_prior_knowledge*(prior_vals/sum(prior_vals))
taxa = length(prior_vals)
samples = %s
sequencing_depth = round(%s)
d = matrix(0,nrow=taxa,ncol=samples)
for (i in 1:samples){
    pvs = rdirichlet(1,alpha=prior_vals)
    d[,i]=table(factor(sample(taxa,sequencing_depth,pvs,replace=TRUE),1:taxa))
}
write.table(d, file=%s, sep=',', quote=FALSE, col.names=FALSE, row.names=FALSE)
"""

def null_from_normal(samples, otus, mean, std, ints=True, clip_on=True,
    sparsity=None):
    """create null from a normal distribution"""
    data = normal(mean, std, samples*otus).reshape(otus, samples)
    if clip_on:
        data = data.clip(0) #replace entries less than 0 with 0
    if ints:
        data = data.round() #default is rounding to nearest int
    if sparsity:
        sparse_inds = arange(otus*samples)
        shuffle(sparse_inds)
        data.flatten()[sparse_inds[:int(otus*samples*sparsity)]] = 0.
        data.reshape(otus, samples)
    return data

def null_from_exponential(samples, otus, scale, ints=True):
    """create null from exponential distribution"""
    data = exponential(scale, samples*otus).reshape(otus, samples)
    if ints:
        data = data.round() #default is rounding to nearest int
    return data

def null_from_data(data, tpk, Rseed=None):
    """generates null from dirichlet model of data based on row sums
    Inputs:
     tpk - int, total prior knowledge to allow the rdirichlet code. higher 
      tpk will result in the simulated table more closely matching the 
      initial data.
     Rseed - int/None, whether or not to seed R at a given value.
    """
    prior_vals = data.sum(1)
    pvs_str =  'c(%s)' % ','.join(map(str,prior_vals))
    sams_str = data.shape[1] # num cols
    seq_depth = data.sum(0).mean()
    out_fp = get_tmp_filename()
    command_str = COMMAND_STR % (pvs_str, tpk, sams_str, seq_depth, out_fp)
    if Rseed!=None:
        command_str = command_str[:23]+'set.seed('+str(Rseed)+')\n'+\
            command_str[23:]
    command_file = get_tmp_filename()
    
    open(command_file, 'w').write(command_str)
    
    cmd_status, cmd_output = getstatusoutput('R --slave < ' + command_file)
    if cmd_status==32512:
        raise ValueError, 'Most likely you do not have R installed, ' +\
            'which is a dependency for QIIME'
    elif cmd_status==256:
        raise ValueError, 'Most likely you do not have gtools library ' +\
            'installed in R installed, which is a dependency for QIIME'
                
    lines = map(str.rstrip , open(out_fp).readlines())
    
    return array([map(float,line.split(',')) for line in lines])