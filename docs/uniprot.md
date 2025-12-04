# Mapping of UniProt to CDM schema

* UniProt schema in XSD format: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot.xsd
* UniProt XML / dat / FASTA file downloads: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/

--

UniProt entries do not have a unique ID as the same protein may have several different accessions.

Therefore, each entry added to the CDM will have a unique ID generated for it
use python function `uuid.uuid4()` (referred to henceforth as `CDM:00000000-0000-0000-12345678`)

The UniProt accessions will be stored in the `identifier` table so that searching by any of the accessions will retrieve the appropriate CDM entry.


TODO:
- add in DataSource table with provenance info for UniProt download - date accessed, where accessed from, etc.

### CDM table `datasource`
```
name: UniProt archaea
source: UniProt
url: https://uniprot.org/downloads/path/to/file
accessed: <today's date>
version: 115
```

Opening stanza of each UniProt entry:
```xml
<entry dataset="Swiss-Prot" created="2025-02-05" modified="2025-04-09" version="115" xmlns="https://uniprot.org/uniprot">
```
TODO: add created / modified / version to CDM schema

### CDM table `entity`
```
entity_id: CDM:00000000-0000-0000-1234567
entity_type: protein
data_source:
created: 2025-02-05
updated: 2025-04-09

```

# UniProt schema top level elements:

## UniProt schema element `accession`
```xml
<xs:element name="accession" type="xs:string" maxOccurs="unbounded"/>
```

### Examples
```xml
    <accession>Q8U0V2</accession>
    <accession>A0A5C0XS91</accession>
```

### CDM table `identifier`

https://kbase.github.io/cdm-schema/Identifier/

```
entity_id: CDM:00000000-0000-0000-1234567
identifier: UniProt:Q8U0V2
source: UniProt
description: UniProt accession

entity_id: CDM:00000000-0000-0000-1234567
identifier: UniProt:A0A5C0XS91
source: UniProt
description: UniProt accession
```

## UniProt schema element `name`
```xml
<xs:element name="name" type="xs:string" maxOccurs="unbounded"/>
```

### Examples
```xml
<name>HIS2_PYRIL</name>
```

### CDM table `name`

https://kbase.github.io/cdm-schema/Name/

```
entity_id: CDM:00000000-0000-0000-1234567
name: HIS2_PYRIL
description: UniProt protein name
source: UniProt
```

## UniProt schema element `protein`
```xml
<xs:element name="protein" type="proteinType"/>
```

### Examples
```xml
    <protein>
        <recommendedName>
            <fullName evidence="1">AMP phosphorylase</fullName>
            <shortName evidence="1">AMPpase</shortName>
            <ecNumber evidence="1">2.4.2.57</ecNumber>
        </recommendedName>
        <alternativeName>
            <fullName evidence="1">Nucleoside monophosphate phosphorylase</fullName>
            <shortName evidence="1">NMP phosphorylase</shortName>
        </alternativeName>
    </protein>
```

### CDM table `name`

https://kbase.github.io/cdm-schema/Name/

Ignore EC numbers (capture them elsewhere)
```
entity_id: CDM:00000000-0000-0000-1234567
name: AMP phosphorylase
description: UniProt recommended full name
source: UniProt

entity_id: CDM:00000000-0000-0000-1234567
name: AMPpase
description: UniProt recommended short name
source: UniProt

entity_id: CDM:00000000-0000-0000-1234567
name: Nucleoside monophosphate phosphorylase
description: UniProt alternative full name
source: UniProt

entity_id: CDM:00000000-0000-0000-1234567
name: NMP phosphorylase
description: UniProt alternative short name
source: UniProt
```

## UniProt schema element `gene`
```xml
<xs:element name="gene" type="geneType" minOccurs="0" maxOccurs="unbounded"/>
```
- ignore

### Examples
```xml
    <gene>
        <name type="primary">gvpN11</name>
        <name evidence="17" type="synonym">gvpN</name>
        <name evidence="15" type="synonym">p-gvpN</name>
        <name type="ordered locus">VNG_5033G</name>
    </gene>
    <gene>
        <name evidence="25" type="primary">gvpN1</name>
        <name evidence="25" type="ordered locus">VNG_6032G</name>
    </gene>
```

## UniProt schema element `organism`
```xml
<xs:element name="organism" type="organismType"/>
```

We only need taxon ID

### Examples
```xml
    <organism>
        <dbReference type="NCBI Taxonomy" id="273063"/>
```


CDM Association table:
```
subject: CDM:00000000-0000-0000-1234567
object: NCBITaxon:273063
```

TODO:
- is there a 'has taxon' predicate?

## UniProt schema element `organismHost`
```xml
<xs:element name="organismHost" type="organismType" minOccurs="0" maxOccurs="unbounded"/>
```
- ignore


## UniProt schema element `geneLocation`
```xml
<xs:element name="geneLocation" type="geneLocationType" minOccurs="0" maxOccurs="unbounded"/>
```
- describes non-nuclear gene locations (organelles and plasmids)
- ignore

## UniProt schema element `reference`
```xml
<xs:element name="reference" type="referenceType" maxOccurs="unbounded"/>
```

### Examples
```xml
    <reference key="1">
        <citation type="journal article" date="2016" name="Genome Announc." volume="4" first="0" last="0">
            <title>Complete genome sequence of the hyperthermophilic and piezophilic archaeon Thermococcus barophilus Ch5, capable of growth at the expense of hydrogenogenesis from carbon monoxide and formate.</title>
            <authorList>
                <person name="Oger P."/>
                <person name="Sokolova T.G."/>
                <person name="Kozhevnikova D.A."/>
                <person name="Taranov E.A."/>
                <person name="Vannier P."/>
                <person name="Lee H.S."/>
                <person name="Kwon K.K."/>
                <person name="Kang S.G."/>
                <person name="Lee J.H."/>
                <person name="Bonch-Osmolovskaya E.A."/>
                <person name="Lebedinsky A.V."/>
            </authorList>
            <dbReference type="PubMed" id="26769929"/>
            <dbReference type="DOI" id="10.1128/genomea.01534-15"/>
        </citation>
        <scope>NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA]</scope>
        <source>
            <strain>Ch5</strain>
        </source>
    </reference>
    <reference key="2">
        <citation type="submission" date="1998-08" db="EMBL/GenBank/DDBJ databases">
            <authorList>
                <person name="Aravalli R.N."/>
            </authorList>
        </citation>
        <scope>NUCLEOTIDE SEQUENCE [GENOMIC DNA]</scope>
        <source>
            <strain>ATCC 35092 / DSM 1617 / JCM 11322 / P2</strain>
        </source>
    </reference>
```

TODO:
- AJ will remap

## UniProt schema element `comment`
```xml
<xs:element name="comment" type="commentType" nillable="true" minOccurs="0" maxOccurs="unbounded"/>
```

Pull out catalytic activity, cofactor, anything else with an ontology ref in the dbReference section.

### Examples
```xml
    <comment type="catalytic activity">
        <reaction evidence="2">
            <text>a primary alcohol + NAD(+) = an aldehyde + NADH + H(+)</text>
            <dbReference type="Rhea" id="RHEA:10736"/>
            <dbReference type="ChEBI" id="CHEBI:15378"/>
            <dbReference type="ChEBI" id="CHEBI:15734"/>
            <dbReference type="ChEBI" id="CHEBI:17478"/>
            <dbReference type="ChEBI" id="CHEBI:57540"/>
            <dbReference type="ChEBI" id="CHEBI:57945"/>
            <dbReference type="EC" id="1.1.1.1"/>
        </reaction>
    </comment>

    <comment type="cofactor">
        <cofactor evidence="1">
            <name>Mg(2+)</name>
            <dbReference type="ChEBI" id="CHEBI:18420"/>
        </cofactor>
        <text evidence="1">Binds 1 Mg(2+) ion per subunit.</text>
    </comment>

```

### CDM table `association`

search for "Rhea" or "EC" in dbReference

RHEA:10736
```
subject: CDM:00000000-0000-0000-1234567
predicate: TODO
object: RHEA:10736
evidence_type: <get from evidence stanzas>
publications:
supporting_objects:
annotation_date:
protocol_id:
primary_knowledge_source:
aggregator:
```

EC:1.1.1.1
```
subject: CDM:00000000-0000-0000-1234567
predicate: ??
object: EC:1.1.1.1
evidence_type: <get from evidence stanzas>
```




TODO:
- evidence: retrieve from `evidence` section
- predicate: depends on what is being captured. Align with GO annotations.


#### List of comment types

```xml
<xs:enumeration value="catalytic activity"/>
<xs:enumeration value="subcellular location"/>
<xs:enumeration value="allergen"/>
<xs:enumeration value="alternative products"/>
<xs:enumeration value="biotechnology"/>
<xs:enumeration value="biophysicochemical properties"/>
<xs:enumeration value="caution"/>
<xs:enumeration value="cofactor"/>
<xs:enumeration value="developmental stage"/>
<xs:enumeration value="disease"/>
<xs:enumeration value="domain"/>
<xs:enumeration value="disruption phenotype"/>
<xs:enumeration value="activity regulation"/>
<xs:enumeration value="function"/>
<xs:enumeration value="induction"/>
<xs:enumeration value="miscellaneous"/>
<xs:enumeration value="pathway"/>
<xs:enumeration value="pharmaceutical"/>
<xs:enumeration value="polymorphism"/>
<xs:enumeration value="PTM"/>
<xs:enumeration value="RNA editing"/>
<xs:enumeration value="similarity"/>
<xs:enumeration value="sequence caution"/>
<xs:enumeration value="subunit"/>
<xs:enumeration value="tissue specificity"/>
<xs:enumeration value="toxic dose"/>
<xs:enumeration value="online information"/>
<xs:enumeration value="mass spectrometry"/>
<xs:enumeration value="interaction"/>
```
TODO:
- which do we want to keep?



## UniProt schema element `dbReference`
```xml
<xs:element name="dbReference" type="dbReferenceType" minOccurs="0" maxOccurs="unbounded"/>
```

### Examples

EC:
```xml
    <dbReference type="EC" id="3.6.1.31" evidence="1"/>
```
dbReference["type"] + ":" + dbReference["id"]
EC:3.6.1.31

Pfam:
```xml
    <dbReference type="Pfam" id="PF01503">
        <property type="entry name" value="PRA-PH"/>
        <property type="match status" value="1"/>
    </dbReference>
```
Pfam:PF01503

PDB:
```xml
    <dbReference type="PDB" id="3CF4">
        <property type="method" value="X-ray"/>
        <property type="resolution" value="2.00 A"/>
        <property type="chains" value="G=1-170"/>
    </dbReference>
```
* PDB:3CF4

TODO: check whether the PDB ID in the dbReference field always appears in the <evidence> section, and that the PDB IDs in the evidence section always appear in the <dbReference> section

RefSeq:
```xml
    <dbReference type="RefSeq" id="WP_011762956.1">
        <property type="nucleotide sequence ID" value="NC_008701.1"/>
    </dbReference>
```
* NCBI protein seq: https://www.ncbi.nlm.nih.gov/protein/WP_011762956.1/


* NCBI nucleotide seq ID: https://www.ncbi.nlm.nih.gov/nuccore/NC_008701.1
(protein has sequence that is contained within this genome seq ID)
```
subject: CDM:....
predicate: is_contained_by ? is_part_of?
object: NC_008701.1
```


TODO:
- is there a relation ontology term for this genome => protein relationship?

see https://obofoundry.org/ontology/ro.html

EMBL:
```xml
    <dbReference type="EMBL" id="CP000504">
        <property type="protein sequence ID" value="ABL88381.1"/>
        <property type="molecule type" value="Genomic_DNA"/>
    </dbReference>
```
* EMBL:CP000504 https://www.ebi.ac.uk/ena/browser/view/CP000504


### CDM table `association`

```
subject: CDM:00000000-0000-0000-12345678
object: EMBL:CP000504
```
```
subject: CDM:00000000-0000-0000-12345678
object: EMBL:CP000504
```
```
subject: CDM:00000000-0000-0000-12345678
object: Pfam:PF01503
```

#### Initial approach

- add in subject/object associations, add in predicate later.

TODO:
- find all possible values for dbReferences
- distinguish between CDM `association`s and `identifier`s
- use ID mapping service to retrieve this data?


## UniProt schema element `proteinExistence`
```xml
<xs:element name="proteinExistence" type="proteinExistenceType"/>
```

### Examples
```xml
    <proteinExistence type="inferred from homology"/>
```

### CDM table `protein`

https://kbase.github.io/cdm-schema/Protein/

see `evidence_for_existence` enum: https://kbase.github.io/cdm-schema/EvidenceForExistence/

```
protein_id: CDM:00000000-0000-0000-1234567
evidence_for_existence: inferred_from_homology
```

## UniProt schema element `keyword`
```xml
<xs:element name="keyword" type="keywordType" minOccurs="0" maxOccurs="unbounded"/>
```
- ignore

## UniProt schema element `feature`
```xml
<xs:element name="feature" type="featureType" minOccurs="0" maxOccurs="unbounded"/>
```
- protein sequence features
- ignore for now
- F-UP

### Examples
```xml
    <feature type="chain" id="PRO_0000136455" description="Histidine biosynthesis bifunctional protein HisIE">
        <location>
            <begin position="1"/>
            <end position="209"/>
        </location>
    </feature>
    <feature type="region of interest" description="Phosphoribosyl-AMP cyclohydrolase">
        <location>
            <begin position="1"/>
            <end position="123"/>
        </location>
    </feature>
```


## UniProt schema element `evidence`
```xml
<xs:element name="evidence" type="evidenceType" minOccurs="0" maxOccurs="unbounded"/>
```

### CDM table `association`

correlate with evidence="x" in the entry's `dbReference` and other fields

### Examples

UniProt:
```xml
    <evidence type="ECO:0000255" key="1">
        <source>
            <dbReference type="HAMAP-Rule" id="MF_01021"/>
        </source>
    </evidence>
```
CDM:
```
evidence_type: ECO:0000255
supporting_objects: HAMAP-Rule:MF_01021
```

UniProt:
```xml
    <evidence type="ECO:0000269" key="2">
        <source>
            <dbReference type="PubMed" id="16042384"/>
        </source>
    </evidence>
```

CDM:
```
evidence_type: ECO:0000269
publications: PMID:16042384
```

UniProt:
```xml
    <evidence type="ECO:0000305" key="3"/>
```
CDM:
```
evidence_type: ECO:0000305
```

UniProt:
```xml
    <evidence type="ECO:0000305" key="4">
        <source>
            <dbReference type="PubMed" id="16042384"/>
        </source>
    </evidence>
```
CDM:
```
evidence_type: ECO:0000305
publications: PMID:16042384
```

UniProt:
```xml
    <evidence type="ECO:0007829" key="5">
        <source>
            <dbReference type="PDB" id="1ZPS"/>
        </source>
    </evidence>
```
CDM:
```
evidence_type: ECO:0007829
supporting_objects: PDB:1ZPS
```

TODO:
- check GAF parser and see what goes in the supporting objects field vs what in the publications field

## UniProt schema element `sequence`
```xml
<xs:element name="sequence" type="sequenceType"/>
```

### Examples
```xml
<sequence length="295" mass="32846" checksum="8D513D8710954F13" modified="2010-05-18" version="1">MADSPFPCVSVIVPVYNDPTGIRDTLTALTKQTYPTERVNILPIDNGSTDETRDVIRQFEREHENVTLVVEDEIQGSYAARNTGIEQATGEIFAFVDADMYMDESWLETAVDAMDEAAYVGCDIELVTNGEDTLPARFDAQTAFPIAQYIRQQQYAPTCGLLVSREVVDDVGPFDERLVSGGDSEFGSRVANAGYRQAFAPAATLYHPVRDSFSSLVKKELRVGRGLCQRQEYYADRFGRPGIPPRPSGVKSPDEESSGLDFERLLFGVLSVVMTGVRALGYYREYVRYVRGNTR</sequence>
```

CDM table `protein`

https://kbase.github.io/cdm-schema/Protein/

```
protein_id: CDM:00000000-0000-0000-1234567
length: 295
sequence: MADSPFPCVSVIVPVYNDPTGIRDTLTALTKQTYPTERVNILPIDNGSTDETRDVIRQFEREHENVTLVVEDEIQGSYAARNTGIEQATGEIFAFVDADMYMDESWLETAVDAMDEAAYVGCDIELVTNGEDTLPARFDAQTAFPIAQYIRQQQYAPTCGLLVSREVVDDVGPFDERLVSGGDSEFGSRVANAGYRQAFAPAATLYHPVRDSFSSLVKKELRVGRGLCQRQEYYADRFGRPGIPPRPSGVKSPDEESSGLDFERLLFGVLSVVMTGVRALGYYREYVRYVRGNTR
```

TODO:
- how is sequence checksum generated?
- what does version refer to?
- ask tech team about file storage -- sequences for proteins, features, scaffolds, etc. - which go into db and which just go as files?
