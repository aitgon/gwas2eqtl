from sqlalchemy import Column, Integer, String, SmallInteger, Float, UniqueConstraint, ForeignKey
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
   pval = Column('pval', String(31), nullable=False)  # comma sep list of pvals
   beta = Column('beta', String(31), nullable=False)  # comma sep list of betas
   se = Column('se', Float, nullable=False)
   eaf = Column('eaf', Float, nullable=True)
   n = Column('n', Integer, nullable=False)
   gwas_id = Column('gwas_id', String(63), primary_key=True)
   pos19 = Column('pos19', Integer, nullable=False)
   Column('coloc_lead_pos', Integer, nullable=False)
   Column('coloc_variant_id', String(50), nullable=False)
   Column('coloc_region', String(50), nullable=False)



class coloc(Base):
   """scripts/insrt_coloc.py"""
   __tablename__ = "coloc"
   __table_args__ = (UniqueConstraint('chrom', 'pos', 'alt', 'eqtl_gene_id', 'gwas_id', 'eqtl_id', name='_coloc_uc'),)

   id = Column('id', String(127), primary_key=True)
   chrom = Column('chrom', SmallInteger, nullable=False)
   pos = Column('pos', Integer, nullable=False)
   rsid= Column('rsid', Integer, nullable=False)
   ref = Column('ref', String(127), nullable=False)
   alt = Column('alt', String(127), nullable=False)
   pval = Column('gwas_pval', Float, nullable=False)
   beta = Column('gwas_beta', Float, nullable=False)
   eqtl_gene_id = Column('eqtl_gene_id', String(15), nullable=False)
   gwas_id = Column('gwas_id', String(127), nullable=False)
   eqtl_pval = Column('eqtl_pval', Float, nullable=False)
   eqtl_beta = Column('eqtl_beta', Float, nullable=False)
   eqtl_id = Column('eqtl_id', String(127), nullable=False)
   pp_h4_abf = Column('pp_h4_abf', Float, nullable=False)
   pp_h3_abf = Column('pp_h3_abf', Float, nullable=False)
   pp_h2_abf = Column('pp_h2_abf', Float, nullable=False)
   pp_h1_abf = Column('pp_h1_abf', Float, nullable=False)
   pp_h0_abf = Column('pp_h0_abf', Float, nullable=False)
   snp_pp_h4 = Column('snp_pp_h4', Float, nullable=False)
   coloc_variant_id = Column('coloc_variant_id', String(63), nullable=False)
   coloc_region = Column('coloc_region', String(63), nullable=False)
   nsnps = Column('nsnps', SmallInteger, nullable=False)
