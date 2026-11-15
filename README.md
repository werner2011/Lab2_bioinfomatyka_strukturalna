# Laboratorium bioinformatyka strukturalna 
Punkt II - część projektowa
Program implementuje w Pythonie algorytm Needlemana–Wunscha, służący do **globalnego dopasowywania dwóch sekwencji biologicznych**. 

Algorytm wyznacza optymalną punktację za pomocą programowania dynamicznego, a następnie odtwarza **wszystkie optymalne dopasowania** metodą rekurencyjną *traceback*.

## 1. Wymagania i instalacja

- Python 3
- Numpy
- Biopython do odczytu plików FASTA (`Bio.SeqIO`)

`argparse` i `pathlib` są częścią standardowej biblioteki Pythona, więc nie wymagają instalacji.

### 1.1. Pobranie projektu


```bash
git clone https://github.com/werner2011/Lab2_bioinfomatyka_strukturalna.git
```

### 1.2. Utworzenie środowiska i instalacja bibliotek

```bash
python -m pip install -r requirements.txt
```

## 2. Uruchamianie 

```bash
python Needleman_wunsch_alg.py examples/para_sekwencji.fasta
```

Zapis do wybranego pliku .txt:

```bash
python Needleman_wunsch_alg.py examples/para_sekwencji.fasta --output output.txt
```

Wyświetlenie macierzy w terminalu

```bash
python Needleman_wunsch_alg.py examples/para_sekwencji.fasta --matrix
```

Informacje o wymaganych argumentach

```bash
python Needleman_wunsch_alg.py --help
```

## 3. Format danych wejściowych i wyjściowych

Program odczytuje **jeden plik FASTA z dokładnie dwoma rekordami**. Pierwszy to `Query`, a drugi `Target`. Przykład `examples/para_sekwencji.fasta`:

```fasta
>query
AACG
>target
AATCG
```

Program zapisuje punktację, liczbę optymalnych dopasowań i wszystkie otrzymane pary sekwencji. Dla wyżej wymienionego pliku oczekiwany wynik to:

```text
Query: AACG
Target: AATCG
Score: 2
Liczba dopasowan: 1

Dopasowanie 1:
AA-CG
AATCG
```

## 4. Metoda i system punktowania

| | Punktacja |
|---|---:|
| Match - identyczne symbole | +1 |
| Mismatch - różne symbole | -1 |
| Gap - każda pozycja luki | -2 |

Algorytm przebiega w trzech etapach:

1. **Inicjalizacja** utworzenie macierzy `(n+1) x (m+1)` i wypełnienie jej pierwszego wiersza oraz kolumny karami za luki.
2. **Programowanie dynamiczne** wypełnianie macierzy maksimum z trzech możliwości: przejścia z przekątnej, z góry i z lewej.
3. **Traceback** cofanie od prawego dolnego rogu do `(0, 0)`


## 5. Przykłady użycia


### 5.1. `para_sekwencji.fasta`

```text
AA-CG
AATCG
```

Cztery identyczne symbole i jedna luka: `4 × (+1) + (−2) = 2`.

```bash
python Needleman_wunsch_alg.py examples/para_sekwencji.fasta --output output.txt
```

### 5.2. `DNA_ex.fasta`

Sekwencje pochodzą ze struktur zdeponowanych w [RCSB PDB](https://www.rcsb.org/):

- [2DND](https://www.rcsb.org/structure/2DND): `CGCAAATTTGCG`
- [1BNA](https://www.rcsb.org/structure/1BNA): `CGCGAATTCGCG` 

Dwa niedopasowania i dziesięć identycznych pozycji dają **Score = 8**.

```bash
python Needleman_wunsch_alg.py examples/DNA_ex.fasta --output dna_output.txt
```

### 5.3. Ubikwityna `ubikwityna.fasta`

Porównanie sekwencji ludzkiej ubikwityny o długości 76 aminokwasów:

- [PDB 1UBQ](https://www.rcsb.org/structure/1UBQ)
- [PDB 2RSU](https://www.rcsb.org/structure/2RSU)

Otrzymano 75 identycznych pozycji i jedną substytucję, bez luk: **Score = 74**.

```bash
python Needleman_wunsch_alg.py examples/ubikwityna.fasta --output protein_output.txt
```

### 5.4. Wiele optymalnych dopasowań `more_align.fasta`


```fasta
>query
AGGTC
>target
AGCTTG
```

Program znajduje **4 różne dopasowania** z tą samą punktacją **−1**:

```text
AG-GTC    AGG-TC    AGGT-C    AGGTC-
AGCTTG    AGCTTG    AGCTTG    AGCTTG
```

```bash
python Needleman_wunsch_alg.py examples/more_align.fasta --output more_align_output.txt
```

### 5.5. Neuropeptydy `oxy_vas.fasta`

Porównano sekwencje aminokwasowe dwóch ludzkich neuropeptydów:

- **Oksytocyna** `CYIQNCPLG`.
- **Argininowa wazopresyna** `CYFQNCPRG`

Oba neuropeptydy mają 9 reszt i różnią się w pozycjach 3 oraz 8. Siedem identycznych pozycji i dwa różniące się aminokwasy dają **Score = 5**. 
Pomimo dużego podobieństwa sekwencji peptydy pełnią różne funkcje fizjologiczne.

```bash
python Needleman_wunsch_alg.py examples/oxy_vas.fasta --output neuropeptydy.txt
```


## 6. Źródła
- [RCSB Protein Data Bank](https://www.rcsb.org/) struktury [1BNA](https://www.rcsb.org/structure/1BNA), [2DND](https://www.rcsb.org/structure/2DND), [1UBQ](https://www.rcsb.org/structure/1UBQ), [2RSU](https://www.rcsb.org/structure/2RSU).
