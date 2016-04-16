# Schedulr
*scheduling as microservice*

## Description

***Schedulr*** is light weighted scheduler or job distribution system between single processes(`workr`). It is based on web techologies and provides microservice which is able to perform job distribution. `Schedulr` also provide frontend service for brief overview and job attaching. Another significant benefit of microservice archetecture is easy accesability(one can manage jobs via one's `ipython` notebook.

## Usage

Just run `Schedulr` in the following manner
```bash
    python schedulr.py
```
and then spawn as many `workrs` as needed
```bash
    python workr.py
```

## Credits 

&copy; [Daniel Bershatsky](mailto:daniel.bershatsky@skolkovotech.ru), 2015-2016