import hashlib


def _hash_string(s):
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


class HashSeq(str):
    def __new__(cls, v):
        # print('validate!!', v)
        instance = super().__new__(cls, v.upper())
        return instance

    @property
    def hash_value(self):
        h = _hash_string(self)
        return h


class HashSeqList(list):
    def append(self, o, /):
        if type(o) is str:
            super().append(HashSeq(o))
        elif type(o) is HashSeq:
            super().append(o)
        else:
            raise ValueError("bad type")

    @property
    def hash_value(self):
        h_list = [x.hash_value for x in self]
        hash_seq = "_".join(sorted(h_list))
        return _hash_string(hash_seq)
