# Laboratorium — algorytm Needlemana–Wunscha

Implementacja **algorytmu Needlemana–Wunscha** do globalnego dopasowywania par sekwencji biologicznych (nukleotydowych/aminokwasowych).

Algorytm wykorzystuje **programowanie dynamiczne** do wyznaczenia maksymalnej punktacji dopasowania. Następnie wykonuje *traceback*, odtwarzając **optymalne dopasowania**.

## 1. Wymagania

- Python 3
- NumPy
- Biopython

Instalacja wymaganych bibliotek:

```bash
python -m pip install -r requirements.txt
```


## 2. Uruchomienie

Uruchomienie z plikiem FASTA:

```bash
python Needleman_wunsch_alg.py examples/para_sekwencji.fasta
```

Zapis do pliku tekstowego:

```bash
python Needleman_wunsch_alg.py examples/para_sekwencji.fasta --output examples/output.txt
```

Wyświetlenie macierzy punktacji w terminalu:

```bash
python Needleman_wunsch_alg.py examples/para_sekwencji.fasta --matrix
```

## 3. Format danych wejściowych

Program oczekuje **jednego pliku FASTA zawierającego dwie sekwencje**. Pierwsza sekwencja jest traktowana jako `Query`, druga jako `Target`.

Przykład — `examples/para_sekwencji.fasta`:

```fasta
>query
AACG
>target
AATCG
```


## 4. System punktowania i zasada działania

| Zdarzenie | Punktacja |
|---|---:|
| Match — identyczne symbole | +1 |
| Mismatch — różne symbole | −1 |
| Gap — pojedyncza pozycja luki | −2 |

Algorytm składa się z trzech etapów:

1. **Inicjalizacja macierzy:** pierwszy wiersz i pierwsza kolumna otrzymują kolejne kary za luki.
2. **Programowanie dynamiczne:** każda komórka przyjmuje maksimum z trzech przejść: po przekątnej, z góry i z lewej.
3. **Traceback:** odtwarzane są ścieżki prowadzące do optymalnej punktacji, w tym wszystkie równorzędne rozwiązania.

Rekordy FASTA są odczytywane przy pomocy `Bio.SeqIO`.

## 5. Przykłady 

W katalogu `examples/` znajdują się trzy pliki FASTA 

### 5.1. Plik — `para_sekwencji.fasta`

Sekwencje: `AACG` (Query) i `AATCG` (Target).

```text
Query:  AA-CG
Target: AATCG
```

- **Wynik:** 2
- **Liczba optymalnych dopasowań:** 1

Uruchomienie:

```bash
python Needleman_wunsch_alg.py examples/para_sekwencji.fasta
```

### 5.2. Sekwencje DNA — `DNA_ex.fasta`

Porównano dwie sekwencje DNA, których struktury znajdują się w bazie PDB:

- **[2DND](https://www.rcsb.org/structure/2DND):** 
- **[1BNA](https://www.rcsb.org/structure/1BNA):** 

```text
2DND: CGCAAATTTGCG
1BNA: CGCGAATTCGCG
```

- **Wynik:** 8
- **Liczba optymalnych dopasowań:** 1

Uruchomienie:

```bash
python Needleman_wunsch_alg.py examples/DNA_ex.fasta --output DNA_output.txt
```

### 5.3. Sekwencje ubikwityny — `ubikwityna.fasta`

Porównano dwie sekwencje o długości **76 aminokwasów**:

- **[1UBQ](https://www.rcsb.org/structure/1UBQ):** ubikwityna typu dzikiego; struktura wyznaczona metodą krystalografii rentgenowskiej.
- **[2RSU](https://www.rcsb.org/structure/2RSU):** wariant ubikwityny **Q41N**; struktura wyznaczona metodą NMR.

Wariant Q41N zawiera substytucję **glutaminy (Q) na asparaginę (N) w pozycji 41**.

- **Wynik:** 74
- **Liczba optymalnych dopasowań:** 1
- **Interpretacja:** 75 identycznych aminokwasów, jedna substytucja, nie wykryto luk.

Uruchomienie:

```bash
python Needleman_wunsch_alg.py examples/ubikwityna.fasta --output ubikwityna_output.txt
```



## Źródła

- [RCSB Protein Data Bank](https://www.rcsb.org/) - rekordy struktur [1BNA](https://www.rcsb.org/structure/1BNA), [2DND](https://www.rcsb.org/structure/2DND), [1UBQ](https://www.rcsb.org/structure/1UBQ) i [2RSU](https://www.rcsb.org/structure/2RSU).
