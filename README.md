# Lab2_bioinfomatyka_strukturalna
Program implementuje algorytm Needlemana–Wunscha do globalnego dopasowywania dwóch sekwencji aminokwasowych.

Algorytm wykorzystuje programowanie dynamiczne do znalezienia dopasowania o maksymalnej punktacji.

## Wymagania
- Python3 
- Numpy 
- Biopython 

Instalacja biblioteke: 

python -m pip install -r requirements.txt 

## Uruchomienie 
python Needleman_wunsch_alg.py examples/para_sekwencji.fasta 

Zapis wyników: 

python Needleman_wunsch_alg.py examples/para_sekwencji.fasta --output examples/output.txt 

Wyświetlenie macierzy: 

python Needleman_wunsch_alg.py examples/para_sekwencji.fasta --matrix 

## Dane 
Program przyjmuje plik FASTA zawierający dwie niepuste sekwencje 

Przykład:
>seq1 
AACG 

>target 
AATCG 

## System punktacji 

Zastosowano taki system punktacji: 
- Match: +1 
- Mismatch: -1 
- Gap: -2 

# Wyniki 

Dla sekwencji AACG i AATCG: 

Score: 2 
Liczba dopasowań: 1 
Znalezione dopasowanie:
AA-CG 
AATCG 

Wyniki są zapisywane domyślnie do pliku output.txt 

Implementacja wykorzystuje NumPy do przechowywania macierzy i Biopython do odczytywania plików FASTA.

