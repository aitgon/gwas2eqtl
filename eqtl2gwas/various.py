#%%
import sqlite3

import numpy
from sqlalchemy import insert


def sql_insert_ignore(table, conn, keys, data_iter):
    """SQL insert or ignore for method args of pandas to_sql"""
    sqlite3.register_adapter(numpy.int64, lambda val: int(val))
    sqlite3.register_adapter(numpy.int32, lambda val: int(val))
    conn.execute(insert(table.table).prefix_with("OR IGNORE"), [dict(zip(keys, row)) for row in data_iter])