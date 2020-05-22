# ImageProcessing

## Gruppo
Gruppo: strcpy(group, "name");

## Descrizione

### Parte 1

- `void ip_mat_free(ip_mat *a)`: Libera memoria se la matrice `a` è non nulla.
- `void ip_mat_init_random(ip_mat *t, float mean, float var)`: Genera una matrice di valori casuali, presi dalla matrice iniziale `t`, utilizzando la distribuzione normale di media `mean` e deviazione standard `std`.

### Parte 2
### Parte 3
### Extra
### Ausiliarie

## Funzioni extra
Ogni eventuale aggiunta extra aggiungetela in fondo a `ip_lib.c` sotto `/*** FUNZIONI EXTRA ***/` e in `ip_lib.h` fate lo stesso.

## Funzioni ausiliarie
In caso di funzioni ausiliarie aggiungere `/**** FUNZIONE AUSILIARIA ***/` sopra la descrizione della funzione in `ip_lib.c`.

## Usage

Per compilare tutto usare
```sh
make all
```
oppure 
```sh
make
```

Per compilare ed eseguire i test
```
make test
./test
```

Per compilare ed eseguire il programma
```sh
make main
./main <params>
```

Sono inclusi due script di test `runner_gauss.sh` e `runner_noise.sh`
Per eseguirli usare
```sh
make main
./runner_<gauss|noise>.sh
```