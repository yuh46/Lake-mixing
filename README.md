# Characterizing vertical mixing in lakes and reservoirs

A parsimonious 1D turbulent diffusion model with an artificial circulation term to simulate both turbulent diffusion and artificial mixing.

## Purpose

Purpose of this software package is to characterize vertical mixing and assess the viability of artificial circulation in lakes and reservoirs.


## Folder structure

Source code of the turbulent diffusion model is stored in the file *TurbulentDiffusionModel*, and the example input data file is in *Inputdata*.

## Prerequisites

TBD

## Example tutorial

### Model set-up

1. Parameter values listed under the *User defined parameters* are calibrated values for Jordan Lake study area. 
2. Updating path of the stored files

```make
data= xlsread('**UPDATE PATH**/Inputdata.xlsx','Input','B2:K19753','',@convertSpreadsheetExcelDates);
```

### Running the model

The model is then run and will generate figures to visualize results (e.g., water temperature at various depths, turbulent diffusion at variopus depths, effective diffusion (including artificial circulation) at various depths).

## References 

Han, Y.; Smithheart, J.W.; Smyth, R.L.; Aziz, T.N.; Obenour, D.R. (*in review*). **Assessing Vertical Mixing as a Potential Control on Cyanobacteria
Dominance in Shallow Turbid Reservoirs**. *Lake and Reservoir Management*.

## Acknowledgements

This study was funded by the National Science Foundation (Project CBET-1545624) and by NC Water Resources Research Institute (Project 16-01-U).


## License

TBD
