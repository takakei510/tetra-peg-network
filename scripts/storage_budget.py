"""Bounded writes for retained analysis, exports and ZIPs; no external packages."""
import csv
import gzip
import io
import json
import shutil
from pathlib import Path
from benchmark_resources import directory_bytes


class StorageLimit(ValueError):
    pass


class WriteBudget:
    """One process, one allocation. Counts writes conservatively, including ZIP rewrites."""
    def __init__(self, limit):
        self.limit = int(limit)
        self.written = 0

    def open(self, path):
        return LimitedFile(path, self)

    def json(self, path, value):
        with self.open(path) as f:
            f.write((json.dumps(value, ensure_ascii=False, indent=2)+'\n').encode())

    def csv(self, path, rows, compressed=False):
        if not rows:
            raise ValueError('No rows')
        fields = list(dict.fromkeys(k for row in rows for k in row))
        with self.open(path) as raw:
            if compressed:
                with gzip.GzipFile(filename='', mode='wb',fileobj=raw,mtime=0,compresslevel=1) as gz:
                    with io.TextIOWrapper(gz,encoding='utf-8-sig',newline='') as text:
                        writer=csv.DictWriter(text,fieldnames=fields)
                        writer.writeheader();writer.writerows(rows)
            else:
                # detach so the outer context owns the file.
                text=io.TextIOWrapper(raw,encoding='utf-8-sig',newline='')
                try:
                    writer=csv.DictWriter(text,fieldnames=fields)
                    writer.writeheader();writer.writerows(rows);text.flush()
                finally:
                    text.detach()


class LimitedFile(io.BufferedIOBase):
    def __init__(self, path, budget):
        self.path = Path(path)
        self.stream = self.path.open('xb')
        self.budget = budget

    def writable(self):return True
    def seekable(self):return True
    def tell(self):return self.stream.tell()
    def seek(self,*args):return self.stream.seek(*args)
    def write(self, data):
        if self.budget.written+len(data)>self.budget.limit:
            raise StorageLimit('Storage write allocation reached; completed samples preserved. Increase --max-storage-gib and rerun.')
        count=self.stream.write(data)
        self.budget.written+=count
        return count
    def flush(self):
        if not self.stream.closed:self.stream.flush()
    def close(self):
        if not self.closed:
            try:super().close()
            finally:self.stream.close()


def export_budget(root, maximum, min_free):
    remaining=min(maximum-directory_bytes(root),shutil.disk_usage(root).free-min_free)
    # Space was reserved before scheduling; fail before creating an export if none remains.
    if remaining<1024*1024:
        raise StorageLimit('Insufficient export space; original completed samples preserved')
    return WriteBudget(remaining)
