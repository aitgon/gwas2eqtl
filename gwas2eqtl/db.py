from sqlalchemy import Column, Integer, String, SmallInteger, Float, UniqueConstraint
from sqlalchemy.orm import declarative_base

# declarative base class
Base = declarative_base()


class tophits(Base):
   """scripts/tophits2db2.py"""
   __tablename__ = "tophits"
   __table_args__ = (UniqueConstraint('chrom', 'pos', 'gwas_id', name='_chrom_pos_uc'),)

   id = Column('id', String(63), primary_key=True)
   chrom = Column('chrom', SmallInteger, nullable=False)
   pos = Column('pos', Integer, nullable=False)
   rsid= Column('rsid', Integer, nullable=False)
   nea = Column('nea', String(255), nullable=False)
   ea = Column('ea', String(255), nullable=False)
   pval = Column('pval', String(32), nullable=False)  # comma sep list of pvals
   beta = Column('beta', String(32), nullable=False)  # comma sep list of betas
   se = Column('se', Float, nullable=False)
   eaf = Column('eaf', Float, nullable=True)
   n = Column('n', Integer, nullable=False)
   gwas_id = Column('gwas_id', String(63), primary_key=True)
   pos19 = Column('pos19', Integer, nullable=False)
