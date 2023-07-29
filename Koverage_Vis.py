import os
import yaml
import pandas as pd
import datapane as dp
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from datetime import datetime

OUTWD = os.getcwd()
SAMCOV = OUTWD + '/results/sample_coverage.tsv'
ALLFILE = OUTWD + '/results/all_coverage.tsv'
CONFYAML = OUTWD + '/koverage.config.yaml'

## Data frame creation ##

DF = pd.read_csv(SAMCOV, sep='\t')
ADF = pd.read_csv(ALLFILE, sep='\t')

## sequence name extraction from results Data frame

seqnames = []
for i in DF['Sample'].unique():
    seqnames.append(i)

###############
# Title
###############

DATE = f"{datetime.now():%d-%m-%Y}"
READSLIST = ", ".join(f"{item}" for item in seqnames)

#Get ref from config.yaml
with open(CONFYAML, 'r') as stream:
    try:
        yaml_content = yaml.safe_load(stream)
        REFVALUE = yaml_content["args"]["ref"]
    except yaml.YAMLError as exc:
        pass


wonk = """# Koverage Report""" + "\n" + "### " + "Report Date: " + DATE + "\n" +"### " + "Reference Sequence: " + REFVALUE + "\n" + "### " + "Reads Processed: " + READSLIST + " ###"
TITLE = dp.Group(dp.Text("""![](https://raw.githubusercontent.com/beardymcjohnface/Koverage/main/koverage.png)"""), dp.Text(wonk), columns=2)

sample_cov_desc =  """| Column | Description |
| --------- | ---------------------------- |
| **Sample** | Sample name derived from read file name |
| **Contig** | Contig ID from assembly FASTA |
| **Count** | Raw mapped read count |
| **RPM** | Reads per million |
| **RPKM** | Reads per kilobase million |
| **RPK** | Reads per kilobase |
| **TPM** | Transcripts per million |
| **Mean** | Estimated mean read depth |
| **Median** | Estimated median read depth |
| **Hitrate** | Estimated fraction of contig with depth > 0 |
| **Variance** | Estimated read depth variance |
<a href="https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/" target="_blank"> Link for further explanation of RPKM and TPM </a>
"""

###########################
# Sample Coverage
###########################

graphs = []

def qualgraph(seqname, df):
    df = df[(df['Sample'] == seqname)]
    df = df.sort_values('Count', ascending=False)
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(go.Bar(x=df['Contig'], y=df['Count'], name="Count"), secondary_y=False)
    fig.add_trace(go.Scatter(x=df['Contig'], y=df['Mean'], name="Mean Depth"), secondary_y=True)
    fig.update_xaxes(title_text=f"{REFVALUE} Contig Number")
    fig.update_yaxes(title_text="Count", secondary_y=False)
    fig.update_yaxes(title_text="Mean", secondary_y=True)
    fig.update_layout(title=seqname)
    graphs.append(dp.Group(dp.Plot(fig), dp.DataTable(df), label=seqname))
    
##regex to shorten contig names - matches megahit and spades assemblies.##

DF['Contig'] = DF['Contig'].str.extract(r'([A-Za-z0-9]+_[A-Za-z0-9]+)')

for plink in seqnames:
    qualgraph(plink, DF)

## changes layout of results (tabbed or not) depending on 1 or more sequence results

if len(seqnames) > 1:
    sample_coverage = dp.Blocks(
    dp.Text("## Sample Coverage"),
    dp.Select(
        blocks=[*graphs]
        )
    )
else:
    sample_coverage = dp.Blocks(
        dp.Text("## Sample Coverage"),
        blocks=[*graphs]
        )
    
SAMPCOV = dp.Group(dp.Text(sample_cov_desc), sample_coverage, label = "Sample Coverage")
VIS = [SAMPCOV]

#############################
# All Coverage
############################

all_cov_desc =  """Column | Description
----- | ----------------------
**Contig** | Contig ID from assembly FASTA
**Count** | Raw mapped read count
**RPM** | Reads per million
**RPKM** | Reads per kilobase million
**RPK** | Reads per kilobase
**TPM** | Transcripts per million
<a href="https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/" target="_blank"> Link for further explanation of RPKM and TPM </a>
"""

##regex to shorten contig names - matches megahit and spades assemblies.##

ADF['Contig'] = ADF['Contig'].str.extract(r'([A-Za-z0-9]+_[A-Za-z0-9]+)')
ADF = ADF.sort_values('Count', ascending=False)

click = []
def butt(val):
    plonk = dict(label=val, method="update", args=[{"y": [ADF[val]]}])
    click.append(plonk)
head = list(ADF)
for i in head:
    butt(i)
fig = go.Figure()
fig.add_trace(go.Bar(x=ADF['Contig'], y=ADF['Count'], name="Count"))
fig.update_xaxes(title_text=f"{REFVALUE} Contig Number")
fig.update_layout(title_text='All Coverage', autosize=True,
    updatemenus=[
        dict(
            type="buttons",
            bgcolor='mediumspringgreen',
            bordercolor='black',
            xanchor="left",
            yanchor="top",
            direction = "left",
            pad={"r": 10, "t": 10},
            x=0.11,
            y=1.4,
            showactive=True,
            buttons=click)])

ALLCOV = dp.Group(dp.Text(all_cov_desc),dp.Plot(fig), dp.DataTable(ADF), label="All Coverage")

VIS.append(ALLCOV)

#################
# Report building
#################

report = dp.Blocks(TITLE, dp.Select(type=dp.SelectType.TABS,blocks=VIS))

dp.save_report(report, path=OUTWD + "/Koverage_Report.html")
