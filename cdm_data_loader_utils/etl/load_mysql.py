import logging
import json
import polars as pl
from sqlalchemy import text
from sqlalchemy.orm import Session

logger = logging.getLogger(__name__)


LOADER_CONFIGURATION = {
    'contigset': [
        ['hash_contigset'],
        {
            'hash_contigset': (None, str),
            'n_contigs': (None, int),
            'size': (None, int),
        }
    ],
    'contig': [
        ['hash_contig'],
        {
            'hash_contig': (None, str),
            'length': (None, int),
            'gc_content': (None, float),
            'base_count': (None, dict)
        }
    ],
    'contigset_x_contig': [
        ['hash_contigset_x_contig'],
        {
            'hash_contigset_x_contig': (None, str),
            'hash_contigset': (None, str),
            'hash_contig': (None, str),
            'contig_index': (None, int)
        }
    ],
    'protocol': [
        ['protocol_id'],
        {
            'protocol_id': (None, str),
            'parent_protocol_id': (None, str),
            'name': (None, str),
            'version': (None, str),
            'identifier': (None, str),
            'inputs': (None, str),
            'outputs': (None, str),
        }
    ],
    'feature': [
        ['feature_id'],
        {
            'feature_id': (None, str),
            'hash_contigset_x_contig': (None, str),
            'source': (None, str),
            'protocol_id': (None, str),
            'type': (None, str),
            'start': (None, int),
            'end': (None, int),
            'strand': (None, str),
            'cds_phase': (None, int),
        }
    ],
    'protein': [
        ['hash_protein_sequence'],
        {
            'hash_protein_sequence': (None, str),
            'length': (None, int),
            'sequence': (None, str)
        }
    ],
    'feature_x_protein': [
        ['feature_id', 'hash_protein_sequence'],
        {
            'feature_id': (None, str),
            'hash_protein_sequence': (None, str),
            'has_stop_codon': (None, bool)
        }
    ]
}


class LoadSQL:

    def __init__(self, engine, table_name: str, columns_pk: list, columns_config: dict):
        self.engine = engine
        self.batch_size = 100000
        self.table_name = table_name
        self.columns_pk = columns_pk
        self.columns_config = columns_config

    @staticmethod
    def build_from_config(engine, table_name, config):
        columns_pk = config[0]
        columns_config = config[1]
        return LoadSQL(engine, table_name, columns_pk, columns_config)

    def load(self, filename, pos=0):
        fn_fetch_keys = self.fetch_keys(self.table_name, self.columns_pk)
        self._load(filename, self.table_name, fn_fetch_keys, self.columns_pk, self.columns_config, pos=pos)

    def fetch_keys(self, table_name, pk):
        def fn_fetch_keys():
            table_ids = set()
            with self.engine.connect() as con:
                select_columns = ', '.join(pk)
                query = f"SELECT {select_columns} FROM {table_name};"
                rs = con.execute(text(query))
                for o in rs:
                    pk_tuple = tuple([o[i] for i in range(len(pk))])
                    table_ids.add(pk_tuple)
            return table_ids

        return fn_fetch_keys

    def load_contigset(self, filename, pos=0):
        table_name = 'contigset'
        pk = ['hash_contigset']
        columns = {
            'hash_contigset': (None, str),
            'n_contigs': (None, int),
            'size': (None, int),
        }
        fn_fetch_keys = self.fetch_keys(table_name, pk)
        self._load(filename, table_name, fn_fetch_keys, pk, columns, pos)

    def load_contig(self, filename, pos=0):
        table_name = 'contig'
        pk = ['hash_contig']
        columns = {
            'hash_contig': (None, str),
            'length': (None, int),
            'gc_content': (None, float),
            'base_count': (None, dict)
        }
        fn_fetch_keys = self.fetch_keys(table_name, pk)
        self._load(filename, table_name, fn_fetch_keys, pk, columns, pos)

    def load_contigset_x_contig(self, filename, pos=0):
        table_name = 'contigset_x_contig'
        pk = ['hash_contigset_x_contig']
        columns = {
            'hash_contigset_x_contig': (None, str),
            'hash_contigset': (None, str),
            'hash_contig': (None, str),
            'contig_index': (None, int)
        }
        fn_fetch_keys = self.fetch_keys(table_name, pk)
        self._load(filename, table_name, fn_fetch_keys, pk, columns, pos)

    @staticmethod
    def process_type(v, t):
        if v is None:
            return 'NULL'
        elif t is str:
            return f"'{v}'"
        elif t is int:
            return int(v)
        elif t is float:
            return float(v)
        elif t is bool:
            if v:
                return 'TRUE'
            else:
                return 'FALSE'
        elif t is dict:
            json_str = json.dumps({k: v for k, v in v.items() if v is not None})
            return f"'{json_str}'"
        raise ValueError(f'invalid type {t}')

    def process_rows(self, table_ids: set, rows: list, columns_pk: list, columns: dict):

        insert_columns = list(columns)

        sql_values = []
        for row in rows:
            logger.debug(row)

            pk = tuple([row[x] for x in columns_pk])

            logger.debug(pk)

            if pk not in table_ids:
                table_ids.add(pk)

                sql_row = []
                for column_name in insert_columns:
                    df_id, column_type = columns[column_name]
                    v = row[column_name] if df_id is None else row[df_id]
                    sql_row.append(str(self.process_type(v, column_type)))

                sql_str = ', '.join(sql_row)

                logger.debug(sql_str)
                sql_values.append(f'({sql_str})')

        return insert_columns, sql_values

    def _load(self, filename, table_name, fn_fetch_keys, columns_pk, columns, pos=0):
        lazy_frame = pl.scan_parquet(filename)
        df_size = lazy_frame.select(pl.len()).collect()['len'][0]
        table_ids = fn_fetch_keys()

        print(df_size, len(table_ids))

        while pos < df_size:
            with Session(self.engine) as session:
                rows = lazy_frame.slice(pos, self.batch_size).collect(streaming=False).rows(named=True)
                insert_columns, sql_values = self.process_rows(table_ids, rows, columns_pk, columns)

                sql_columns = ', '.join(insert_columns)
                if len(sql_values) > 0:
                    sql_insert = f"""
                    INSERT INTO {table_name} 
                      ({sql_columns}) 
                    VALUES 
                    """
                    sql_insert += ',\n'.join(sql_values)
                    sql_insert += ';'
                    print('payload len:', len(sql_insert))
                    session.execute(text(sql_insert))
                    session.commit()
            pos += self.batch_size
            print(pos)
