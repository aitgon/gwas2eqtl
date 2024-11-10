from sqlalchemy import Column, Integer, String, Float, UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

# Declarative meta class
class Coloc(Base):

    __tablename__ = 'coloc'
    __table_args__ = (
            UniqueConstraint('rsid', 'egene', 'eqtl_identifier', 'gwas_identifier', 'leadeqtl_pos', 'leadeqtl_egene'),
            )

    id = Column(Integer, primary_key=True, autoincrement=True)
    chrom = Column(Integer)
    pos = Column(Integer)
    rsid = Column(String)
    ref = Column(String)
    alt = Column(String)
    egene_symbol = Column(String)
    egene = Column(String)
    eqtl_beta = Column(Float)
    eqtl_pvalue = Column(Float)
    eqtl_identifier = Column(String)
    gwas_beta = Column(Float)
    gwas_pvalue = Column(Float)
    gwas_identifier = Column(String)
    gwas_trait_name = Column(String)
    pp_h4 = Column(Float)
    PP_H4_abf = Column('PP.H4.abf', Float)
    nsnps = Column(Integer)
    PP_H3_abf = Column('PP.H3.abf',Float)
    PP_H2_abf = Column('PP.H2.abf', Float)
    PP_H1_abf = Column('PP.H1.abf', Float)
    PP_H0_abf = Column('PP.H0.abf', Float)
    leadeqtl_pos = Column(Integer)
    leadeqtl_egene = Column(String)
