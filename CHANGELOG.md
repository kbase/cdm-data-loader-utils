# CDM Data Loaders Changelog

- [CDM Data Loaders Changelog](#cdm-data-loaders-changelog)
    - [v0.1.12](#v0112)
    - [v0.1.11](#v0111)
    - [v0.1.10](#v0110)
    - [v0.1.9](#v019)
    - [v0.1.8](#v018)
    - [v0.1.7](#v017)
    - [v0.1.6](#v016)
    - [v0.1.5](#v015)
    - [v0.1.4](#v014)
    - [v0.1.3](#v013)
    - [v0.1.2](#v012)
    - [v0.1.1](#v011)
    - [v0.1.0](#v010)


### v0.1.12

- Dedicated UniProt and UniRef XML import pipelines using custom XSD-aware XML-to-JSON converter.

### v0.1.11

- Preliminary XSV, XML, and JSONLines importers

### v0.1.10

- Add in the ability to configure loggers using a file
- Bugfix for UniProt pipeline

### v0.1.9

- Splitting out NCBI REST API pipeline into annotation and dataset components.
- NCBI FTP pipeline and accompanying notebooks added.

### v0.1.8

- Add rotating file log handler for easier debugging.

### v0.1.7

- Add in AllTheBacteria file download client.

### v0.1.6

- Make NCBI REST API client more resilient to errors and ensure existing imports are not lost.

### v0.1.5

- Add batch size parameter to the NCBI REST API interface.


### v0.1.4

- Add in NCBI REST API interface.


### v0.1.3

- Add in file batcher for use with file-based importers.


### v0.1.2

- Update XML File Splitter to use the latest version, which includes the `gzip` parameter.

### v0.1.1

- Add [XML File Splitter](https://github.com/ialarmedalien/xml_file_splitter) to the container.

### v0.1.0

- Initial release.
