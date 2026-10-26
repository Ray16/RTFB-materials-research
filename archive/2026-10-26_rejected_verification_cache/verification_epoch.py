"""Reuse large immutable-input verification within one CPU dispatch pass.

Small source/configuration files always receive a full hash check. Large binary
objects are checked again after any metadata change, including ctime, and after
each explicit epoch boundary. This helper never runs in scientific workers.
"""
import json
from pathlib import Path
import time


class VerificationEpoch:
    def __init__(self,verify,minimum_bytes=16*1024**2,report=None):
        self.verify_original=verify;self.minimum_bytes=minimum_bytes;self.report=report
        self.epoch=0;self.cache={};self.reset_counts()

    def reset_counts(self):
        self.counts=dict(full_checks=0,large_binary_full_checks=0,large_binary_bytes_hashed=0,
            large_binary_checks_reused=0,large_binary_bytes_not_rehashed=0,full_check_seconds=0.)

    @staticmethod
    def stamp(path):
        stat=path.stat()
        return (stat.st_dev,stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns)

    def begin(self):
        if self.report and self.epoch:
            with self.report.open('a') as stream:stream.write(json.dumps(dict(epoch=self.epoch,checked=time.time(),**self.counts))+'\n')
        self.epoch+=1;self.cache.clear();self.reset_counts()

    def verify(self,mapping):
        for original,expected in mapping.items():
            path=Path(original).resolve();before=self.stamp(path)
            large=path.suffix.lower() in ('.npy','.h5') and before[2]>=self.minimum_bytes
            key=(str(path),expected)
            if large and self.cache.get(key)==before:
                self.counts['large_binary_checks_reused']+=1
                self.counts['large_binary_bytes_not_rehashed']+=before[2];continue
            started=time.perf_counter();self.verify_original({original:expected})
            self.counts['full_check_seconds']+=time.perf_counter()-started
            after=self.stamp(path)
            if before!=after:raise RuntimeError('Verification input changed while hashing: '+str(path))
            self.counts['full_checks']+=1
            if large:
                self.cache[key]=after;self.counts['large_binary_full_checks']+=1
                self.counts['large_binary_bytes_hashed']+=after[2]
