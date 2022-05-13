from gwastoeqtl.BaseColoc import Base, Coloc
from pandas import ExcelWriter
from sqlalchemy import create_engine
import csv
import os
import pandas
import pathlib
import sys

#%%############################################################################
#
# some Parameters
#
###############################################################################
from gwastoeqtl.Logger import Logger
from gwastoeqtl.UCSCpysql import UCSCpysql

eqtl_identifier_tsv_path = "config/eqtl_Kasela_2017_CD8.ods"
gwas_identifier_ods_path = "config/gwas_ieu-a1162.ods"
coloc_dir = "results/coloc/1:121000000-122000000/5e-08/500000"
tsv_path = os.path.join(coloc_dir, "cat", "coloc.tsv")
ods_path = os.path.join(coloc_dir, "cat", "coloc.ods")
db_path = os.path.join("results", "db.sqlite")


#%%
help_cmd_str = "todo"
try:
    eqtl_identifier_tsv_path = sys.argv[1]
    gwas_identifier_ods_path = sys.argv[2]
    coloc_dir = sys.argv[3]
    tsv_path = sys.argv[4]
    ods_path = sys.argv[5]
    db_path = sys.argv[6]
    if len(sys.argv) > 7:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

###############################################################################
#
# select immune cell types
#
###############################################################################

#%%
eqtl_df = pandas.read_csv(eqtl_identifier_tsv_path, sep="\t")
eqtl_df = eqtl_df.loc[eqtl_df['ftp_path'].str.contains('/ge/|/microarray/', regex=True, na=False), ]
eqtl_identifier_lst = (eqtl_df['ftp_path'].str.replace('.all.tsv.gz', '', regex=True)).str.split('/', expand=True)[10].tolist()

#%%
gwas_identifier_df = pandas.read_excel(gwas_identifier_ods_path, header=0)
gwas_identifier_lst = gwas_identifier_df["identifier"].tolist()

#%%############################################################################
#
# Concatenate all colocalisation studies H3
#
###############################################################################
coloc_df = pandas.DataFrame()
Logger.info("len(eqtl_identifier_lst): " + str(len(eqtl_identifier_lst)))
Logger.info("len(gwas_identifier_lst): " + str(len(gwas_identifier_lst)))
coloc_line_lst = []
for eqtl_identifier in eqtl_identifier_lst:
    eqtl_dir_path = os.path.join(os.getcwd(), coloc_dir, eqtl_identifier)
    for gwas_identifier in gwas_identifier_lst:
        coloc_tsv_path = os.path.join(eqtl_dir_path, "{}.tsv".format(gwas_identifier))
        if os.path.isfile(coloc_tsv_path) and pathlib.Path(coloc_tsv_path).stat().st_size > 0:
            # csv reader is much faster than pandas
            with open(coloc_tsv_path) as fin:
                csvreader = csv.reader(fin, delimiter="\t")
                columns = next(csvreader)
                for line in csvreader:
                    coloc_line_lst.append(line)

coloc_df = pandas.DataFrame.from_records(coloc_line_lst, columns=columns)
Logger.info("coloc_df.shape[0]: " + str(coloc_df.shape[0]))
columns = ['chrom', 'pos', 'rsid', 'ref', 'alt', 'egene_symbol', 'egene', 'eqtl_beta', 'eqtl_pvalue', 'eqtl_identifier',
     'gwas_beta', 'gwas_pvalue', 'gwas_identifier', 'gwas_trait_name', 'pp_h4', 'PP.H4.abf', 'nsnps', 'PP.H3.abf',
     'PP.H2.abf', 'PP.H1.abf', 'PP.H0.abf', 'leadeqtl_pos', 'leadeqtl_egene', ]
if coloc_df.shape[0] > 0:
    # %%
    Logger.info("Annotate genes")
    database = 'hg38'
    sql = 'select distinct knownAttrs.geneId, kgXref.geneSymbol from knownAttrs, kgXref where knownAttrs.kgID=kgXref.kgID'
    df = UCSCpysql(sql=sql, db=database).get_annotation_df()
    df.rename({'geneId': 'egene', 'geneSymbol': 'egene_symbol'}, axis=1, inplace=True)
    df['egene'] = df['egene'].str.split('.').str[0]
    coloc_df = coloc_df.merge(df, on='egene', how='left')
    coloc_df = coloc_df[columns]
    coloc_df.drop_duplicates(inplace=True)
    coloc_df.sort_values(by=['PP.H4.abf', 'pp_h4', 'chrom', 'pos'], ascending=[False, False, True, True], inplace=True)

    #%% TSV
    Logger.info("Write as TSV")
    pathlib.Path(os.path.dirname(ods_path)).mkdir(parents=True, exist_ok=True)
    coloc_df.to_csv(tsv_path, sep="\t", index=False, header=True)

    #%% ODS
    Logger.info("Write as ODS")
    pathlib.Path(os.path.dirname(ods_path)).mkdir(parents=True, exist_ok=True)
    coloc_df["PP.H4.abf"] = pandas.to_numeric(coloc_df["PP.H4.abf"])
    coloc_cutoff_df = coloc_df.loc[coloc_df["PP.H4.abf"] >= 0.8, ]
    # import pdb; pdb.set_trace()
    with ExcelWriter(ods_path) as fout:
        coloc_cutoff_df.to_excel(fout, index=False, engine="odfpy")

    #%% SQLITE
    Logger.info("Insert into SQLITE")
    pathlib.Path(os.path.dirname(db_path)).mkdir(exist_ok=True, parents=True)
    engine_str = 'sqlite:///{db_path}'.format(db_path=db_path)
    engine = create_engine(engine_str, echo=False)
    Base.metadata.create_all(engine)
    coloc_cutoff_dic_lst = coloc_cutoff_df.to_dict('records')
    with engine.begin() as conn:
        conn.execute(Coloc.__table__.insert().prefix_with("OR IGNORE"), coloc_cutoff_dic_lst)
else:
    print("Warning: No colocalisations!!")
    with ExcelWriter(ods_path) as writer:
        pandas.DataFrame(columns=columns).to_excel(writer, index=False)
