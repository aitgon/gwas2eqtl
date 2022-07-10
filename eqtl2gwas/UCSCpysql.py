import pandas
from sqlalchemy import create_engine

class UCSCpysql:
    """This class easily get UCSC tables for annotation based on SQL queries"""

    host = 'genome-euro-mysql.soe.ucsc.edu'
    user = 'genome'

    def __init__(self, sql, db):

        self.sql = sql
        self.database = db
        conn_dic = {'host': UCSCpysql.host, 'user': UCSCpysql.user, 'db': db}
        self.engine = create_engine('mysql+pymysql://{user}@{host}/{db}'.format(**conn_dic))

    def get_annotation_df(self):
        """Returns a DF with the given annotation"""

        with self.engine.connect() as con:
            return pandas.read_sql(sql=self.sql, con=con)
        return None
